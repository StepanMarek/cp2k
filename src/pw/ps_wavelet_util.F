!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Performs a wavelet based solution of the Poisson equation.
!> \author Florian Schiffmann (09.2007,fschiff)
! **************************************************************************************************
MODULE ps_wavelet_util

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE ps_wavelet_base,                 ONLY: f_poissonsolver,&
                                              p_poissonsolver,&
                                              s_poissonsolver
   USE ps_wavelet_fft3d,                ONLY: fourier_dim
   USE pw_grid_types,                   ONLY: pw_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_util'

   ! *** Public data types ***

   PUBLIC :: PSolver, &
             P_FFT_dimensions, &
             S_FFT_dimensions, &
             F_FFT_dimensions

CONTAINS

! **************************************************************************************************
!> \brief Calculate the Poisson equation $\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)$
!>     from a given $\rho$, for different boundary conditions an for different data distributions.
!>     Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!> \param geocode Indicates the boundary conditions (BC) of the problem:
!>             'F' free BC, isolated systems.
!>                 The program calculates the solution as if the given density is
!>                 "alone" in R^3 space.
!>             'S' surface BC, isolated in y direction, periodic in xz plane
!>                 The given density is supposed to be periodic in the xz plane,
!>                 so the dimensions in these direction mus be compatible with the FFT
!>                 Beware of the fact that the isolated direction is y!
!>             'P' periodic BC.
!>                 The density is supposed to be periodic in all the three directions,
!>                 then all the dimensions must be compatible with the FFT.
!>                 No need for setting up the kernel.
!> \param iproc label of the process,from 0 to nproc-1
!> \param nproc number of processors
!> \param n01 global dimension in the three directions.
!> \param n02 global dimension in the three directions.
!> \param n03 global dimension in the three directions.
!> \param hx    grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param hy grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param hz grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param rhopot main input/output array.
!>                 On input, it represents the density values on the grid points
!>                 On output, it is the Hartree potential, namely the solution of the Poisson
!>                 equation PLUS (when ixc/=0) the XC potential PLUS (again for ixc/=0) the
!>                 pot_ion array. The output is non overlapping, in the sense that it does not
!>                 consider the points that are related to gradient and WB calculation
!> \param karray kernel of the poisson equation. It is provided in distributed case, with
!>                 dimensions that are related to the output of the PS_dim4allocation routine
!>                 it MUST be created by following the same geocode as the Poisson Solver.
!> \param pw_grid ...
!> \date February 2007
!> \author Luigi Genovese
!> \note The dimensions of the arrays must be compatible with geocode, nproc,
!>     ixc and iproc. Since the arguments of these routines are indicated with the *, it
!>     is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
! **************************************************************************************************
   SUBROUTINE PSolver(geocode, iproc, nproc, n01, n02, n03, hx, hy, hz, &
                      rhopot, karray, pw_grid)
      CHARACTER(len=1), INTENT(in)                       :: geocode
      INTEGER, INTENT(in)                                :: iproc, nproc, n01, n02, n03
      REAL(KIND=dp), INTENT(in)                          :: hx, hy, hz
      REAL(KIND=dp), DIMENSION(*), INTENT(inout)         :: rhopot
      REAL(KIND=dp), DIMENSION(*), INTENT(in)            :: karray
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      INTEGER                                            :: i1, i2, i3, iend, istart, j2, m1, m2, &
                                                            m3, md1, md2, md3, n1, n2, n3, nd1, &
                                                            nd2, nd3, nlim, nwb, nwbl, nwbr, nxc, &
                                                            nxcl, nxcr, nxt
      REAL(KIND=dp)                                      :: factor, hgrid, red_fact, scal
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: zf

!the order of the finite-difference gradient (fixed)
!calculate the dimensions wrt the geocode

      IF (geocode == 'P') THEN
         CALL P_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      ELSE IF (geocode == 'S') THEN
         CALL S_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      ELSE IF (geocode == 'F') THEN
         CALL F_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      ELSE
         CPABORT("PSolver: geometry code not admitted")
      END IF
      !array allocations
      ALLOCATE (zf(md1, md3, md2/nproc))

      !  CALL timing(iproc,'Exchangecorr  ','ON')
      !dimension for exchange-correlation (different in the global or distributed case)
      !let us calculate the dimension of the portion of the rhopot array to be passed
      !to the xc routine
      !this portion will depend on the need of calculating the gradient or not,
      !and whether the White-Bird correction must be inserted or not
      !(absent only in the LB ixc=13 case)

      !nxc is the effective part of the third dimension that is being processed
      !nxt is the dimension of the part of rhopot that must be passed to the gradient routine
      !nwb is the dimension of the part of rhopot in the wb-postprocessing routine
      !note: nxc <= nwb <= nxt
      !the dimension are related by the values of nwbl and nwbr
      !      nxc+nxcl+nxcr-2 = nwb
      !      nwb+nwbl+nwbr = nxt
      istart = iproc*(md2/nproc)
      iend = MIN((iproc + 1)*md2/nproc, m2)

      nxc = iend - istart
      nwbl = 0
      nwbr = 0
      nxcl = 1
      nxcr = 1

      nwb = nxcl + nxc + nxcr - 2
      nxt = nwbr + nwb + nwbl

      !calculate the actual limit of the array for the zero padded FFT
      IF (geocode == 'P') THEN
         nlim = n2
      ELSE IF (geocode == 'S') THEN
         nlim = n2
      ELSE IF (geocode == 'F') THEN
         nlim = n2/2
      END IF

      !!$  print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
      !!$  print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1

      IF (istart + 1 <= m2) THEN
         red_fact = 1._dp
         CALL scale_and_distribute(m1, m3, md1, md2, md3, nxc, rhopot, zf, nproc, red_fact)
      ELSE IF (istart + 1 <= nlim) THEN !this condition assures that we have perform good zero padding
         DO i2 = istart + 1, MIN(nlim, istart + md2/nproc)
            j2 = i2 - istart
            DO i3 = 1, md3
               DO i1 = 1, md1
                  zf(i1, i3, j2) = 0._dp
               END DO
            END DO
         END DO
      END IF

      !this routine builds the values for each process of the potential (zf), multiplying by scal
      IF (geocode == 'P') THEN
         !no powers of hgrid because they are incorporated in the plane wave treatment
         scal = 1._dp/REAL(n1*n2*n3, KIND=dp)
         CALL P_PoissonSolver(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, zf, &
                              scal, hx, hy, hz, pw_grid%para%group)
      ELSE IF (geocode == 'S') THEN
         !only one power of hgrid
         scal = hy/REAL(n1*n2*n3, KIND=dp)
         CALL S_PoissonSolver(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, karray, zf, &
                              scal, pw_grid%para%group)
      ELSE IF (geocode == 'F') THEN
         hgrid = MAX(hx, hy, hz)
         scal = hgrid**3/REAL(n1*n2*n3, KIND=dp)
         CALL F_PoissonSolver(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, karray, zf, &
                              scal, pw_grid%para%group)
         factor = 0.5_dp*hgrid**3
      END IF

      !  call timing(iproc,'PSolv_comput  ','ON')

      !the value of the shift depends on the distributed i/o or not
      IF (geocode == 'F') THEN
         red_fact = 1._dp
      ELSE
         red_fact = -fourpi
      END IF

      CALL scale_and_distribute(m1, m3, md1, md2, md3, nxc, zf, rhopot, nproc, red_fact)

      DEALLOCATE (zf)

   END SUBROUTINE PSolver

! **************************************************************************************************
!> \brief Calculate four sets of dimension needed for the calculation of the
!>     convolution for the periodic system
!> \param n01 original real dimensions (input)
!> \param n02 original real dimensions (input)
!> \param n03 original real dimensions (input)
!> \param m1 original real dimension, with m2 and m3 exchanged
!> \param m2 original real dimension, with m2 and m3 exchanged
!> \param m3 original real dimension, with m2 and m3 exchanged
!> \param n1 the first FFT dimensions, for the moment supposed to be even
!> \param n2 the first FFT dimensions, for the moment supposed to be even
!> \param n3 the first FFT dimensions, for the moment supposed to be even
!> \param md1 the n1,n2,n3 dimensions. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param md2 the n1,n2,n3 dimensions. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param md3 the n1,n2,n3 dimensions. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param nd1 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd2 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd3 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nproc ...
!> \date October 2006
!> \author Luigi Genovese
!> \note This four sets of dimensions are actually redundant (mi=n0i),
!>     due to the backward-compatibility
!>     with the other geometries of the Poisson Solver.
!>     The dimensions 2 and 3 are exchanged.
! **************************************************************************************************
   SUBROUTINE P_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      INTEGER, INTENT(in)                                :: n01, n02, n03
      INTEGER, INTENT(out)                               :: m1, m2, m3, n1, n2, n3, md1, md2, md3, &
                                                            nd1, nd2, nd3
      INTEGER, INTENT(in)                                :: nproc

      INTEGER                                            :: l1, l2, l3

!dimensions of the density in the real space

      m1 = n01
      m2 = n03
      m3 = n02

      ! real space grid dimension (suitable for number of processors)
      l1 = m1
      l2 = m2
      l3 = m3 !beware of the half dimension
      CALL fourier_dim(l1, n1)
      IF (n1 == m1) THEN
      ELSE
         PRINT *, 'the FFT in the x direction is not allowed'
         PRINT *, 'n01 dimension', n01
         CPABORT("")
      END IF
      l1 = l1 + 1
      CALL fourier_dim(l2, n2)
      IF (n2 == m2) THEN
      ELSE
         PRINT *, 'the FFT in the z direction is not allowed'
         PRINT *, 'n03 dimension', n03
         CPABORT("")
      END IF
      CALL fourier_dim(l3, n3)
      IF (n3 == m3) THEN
      ELSE
         PRINT *, 'the FFT in the y direction is not allowed'
         PRINT *, 'n02 dimension', n02
         CPABORT("")
      END IF

      !dimensions that contain the unpadded real space,
      ! compatible with the number of processes
      md1 = n1
      md2 = n2
      md3 = n3
      DO WHILE (nproc*(md2/nproc) .LT. n2)
         md2 = md2 + 1
      END DO

      !dimensions of the kernel, 1/8 of the total volume,
      !compatible with nproc
      nd1 = n1/2 + 1
      nd2 = n2/2 + 1
      nd3 = n3/2 + 1
      DO WHILE (MODULO(nd3, nproc) .NE. 0)
         nd3 = nd3 + 1
      END DO

   END SUBROUTINE P_FFT_dimensions

! **************************************************************************************************
!> \brief Calculate four sets of dimension needed for the calculation of the
!>     convolution for the surface system
!> \param n01 original real dimensions (input)
!> \param n02 original real dimensions (input)
!> \param n03 original real dimensions (input)
!> \param m1 original real dimension, with 2 and 3 exchanged
!> \param m2 original real dimension, with 2 and 3 exchanged
!> \param m3 original real dimension, with 2 and 3 exchanged
!> \param n1 the first FFT dimensions, for the moment supposed to be even
!> \param n2 the first FFT dimensions, for the moment supposed to be even
!> \param n3 the double of the first FFT even dimension greater than m3
!>           (improved for the HalFFT procedure)
!> \param md1 the n1,n2 dimensions.
!> \param md2 the n1,n2,n3 dimensions.
!> \param md3 the half of n3 dimension. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param nd1 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd2 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd3 fourier dimensions for which the kernel is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nproc ...
!> \date October 2006
!> \author Luigi Genovese
!> \note This four sets of dimensions are actually redundant (mi=n0i),
!>     due to the backward-compatibility
!>     with the Poisson Solver with other geometries.
!>     Dimensions n02 and n03 were exchanged
! **************************************************************************************************
   SUBROUTINE S_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      INTEGER, INTENT(in)                                :: n01, n02, n03
      INTEGER, INTENT(out)                               :: m1, m2, m3, n1, n2, n3, md1, md2, md3, &
                                                            nd1, nd2, nd3
      INTEGER, INTENT(in)                                :: nproc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'S_FFT_dimensions'

      INTEGER                                            :: handle, l1, l2, l3

!dimensions of the density in the real space

      CALL timeset(routineN, handle)
      m1 = n01
      m2 = n03
      m3 = n02

      ! real space grid dimension (suitable for number of processors)
      l1 = m1
      l2 = m2
      l3 = m3 !beware of the half dimension
      CALL fourier_dim(l1, n1)
      IF (n1 == m1) THEN
      ELSE
         PRINT *, 'the FFT in the x direction is not allowed'
         PRINT *, 'n01 dimension', n01
         CPABORT("")
      END IF
      l1 = l1 + 1
      CALL fourier_dim(l2, n2)
      IF (n2 == m2) THEN
      ELSE
         PRINT *, 'the FFT in the z direction is not allowed'
         PRINT *, 'n03 dimension', n03
         CPABORT("")
      END IF
      DO
         CALL fourier_dim(l3, n3)
         IF (MODULO(n3, 2) == 0) THEN
            EXIT
         END IF
         l3 = l3 + 1
      END DO
      n3 = 2*n3

      !dimensions that contain the unpadded real space,
      ! compatible with the number of processes
      md1 = n1
      md2 = n2
      md3 = n3/2
      DO WHILE (nproc*(md2/nproc) .LT. n2)
         md2 = md2 + 1
      END DO

      !dimensions of the kernel, 1/8 of the total volume,
      !compatible with nproc

      !these two dimensions are like that since they are even
      nd1 = n1/2 + 1
      nd2 = n2/2 + 1

      nd3 = n3/2 + 1
      DO WHILE (MODULO(nd3, nproc) .NE. 0)
         nd3 = nd3 + 1
      END DO
      CALL timestop(handle)

   END SUBROUTINE S_FFT_dimensions

! **************************************************************************************************
!> \brief Calculate four sets of dimension needed for the calculation of the
!>     zero-padded convolution
!> \param n01 original real dimensions (input)
!> \param n02 original real dimensions (input)
!> \param n03 original real dimensions (input)
!> \param m1 original real dimension with the dimension 2 and 3 exchanged
!> \param m2 original real dimension with the dimension 2 and 3 exchanged
!> \param m3 original real dimension with the dimension 2 and 3 exchanged
!> \param n1 ...
!> \param n2 ...
!> \param n3 the double of the first FFT even dimension greater than m3
!>           (improved for the HalFFT procedure)
!> \param md1 half of n1,n2,n3 dimension. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param md2 half of n1,n2,n3 dimension. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!>                 properly enlarged to be compatible with the FFT dimensions n_i.
!>                 md2 is further enlarged to be a multiple of nproc
!> \param nd1 fourier dimensions for which the kernel FFT is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd2 fourier dimensions for which the kernel FFT is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nd3 fourier dimensions for which the kernel FFT is injective,
!>                 formally 1/8 of the fourier grid. Here the dimension nd3 is
!>                 enlarged to be a multiple of nproc
!> \param nproc ...
!> \date February 2006
!> \author Luigi Genovese
!> \note The dimension m2 and m3 correspond to n03 and n02 respectively
!>     this is needed since the convolution routine manage arrays of dimension
!>     (md1,md3,md2/nproc)
! **************************************************************************************************
   SUBROUTINE F_FFT_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      INTEGER, INTENT(in)                                :: n01, n02, n03
      INTEGER, INTENT(out)                               :: m1, m2, m3, n1, n2, n3, md1, md2, md3, &
                                                            nd1, nd2, nd3
      INTEGER, INTENT(in)                                :: nproc

      INTEGER                                            :: l1, l2, l3

!dimensions of the density in the real space, inverted for convenience

      m1 = n01
      m2 = n03
      m3 = n02
      ! real space grid dimension (suitable for number of processors)
      l1 = 2*m1
      l2 = 2*m2
      l3 = m3 !beware of the half dimension
      DO
         CALL fourier_dim(l1, n1)
         IF (MODULO(n1, 2) == 0) THEN
            EXIT
         END IF
         l1 = l1 + 1
      END DO
      DO
         CALL fourier_dim(l2, n2)
         IF (MODULO(n2, 2) == 0) THEN
            EXIT
         END IF
         l2 = l2 + 1
      END DO
      DO
         CALL fourier_dim(l3, n3)
         IF (MODULO(n3, 2) == 0) THEN
            EXIT
         END IF
         l3 = l3 + 1
      END DO
      n3 = 2*n3

      !dimensions that contain the unpadded real space,
      ! compatible with the number of processes
      md1 = n1/2
      md2 = n2/2
      md3 = n3/2
      DO WHILE (nproc*(md2/nproc) .LT. n2/2)
         md2 = md2 + 1
      END DO

      !dimensions of the kernel, 1/8 of the total volume,
      !compatible with nproc
      nd1 = n1/2 + 1
      nd2 = n2/2 + 1
      nd3 = n3/2 + 1

      DO WHILE (MODULO(nd3, nproc) .NE. 0)
         nd3 = nd3 + 1
      END DO

   END SUBROUTINE F_FFT_dimensions

! **************************************************************************************************
!> \brief ...
!> \param m1 ...
!> \param m3 ...
!> \param md1 ...
!> \param md2 ...
!> \param md3 ...
!> \param nxc ...
!> \param rhopot ...
!> \param zf ...
!> \param nproc ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE scale_and_distribute(m1, m3, md1, md2, md3, nxc, &
                                   rhopot, zf, nproc, factor)

      !Arguments----------------------
      INTEGER, INTENT(in)                                :: m1, m3, md1, md2, md3, nxc, nproc
      REAL(KIND=dp), DIMENSION(md1, md3, md2/nproc), &
         INTENT(inout)                                   :: zf, rhopot
      REAL(KIND=dp), INTENT(in)                          :: factor

      CHARACTER(len=*), PARAMETER :: routineN = 'scale_and_distribute'

      INTEGER                                            :: handle, j1, j3, jp2

      CALL timeset(routineN, handle)

      IF (nxc .GE. 1) THEN
         DO jp2 = 1, nxc
            DO j3 = 1, m3
               DO j1 = 1, m1
                  zf(j1, j3, jp2) = factor*rhopot(j1, j3, jp2)
               END DO
               DO j1 = m1 + 1, md1
                  zf(j1, j3, jp2) = 0._dp
               END DO
            END DO
            DO j3 = m3 + 1, md3
               DO j1 = 1, md1
                  zf(j1, j3, jp2) = 0._dp
               END DO
            END DO
         END DO
         DO jp2 = nxc + 1, md2/nproc
            DO j3 = 1, md3
               DO j1 = 1, md1
                  zf(j1, j3, jp2) = 0._dp
               END DO
            END DO
         END DO
      ELSE
         zf = 0._dp
      END IF
      CALL timestop(handle)

   END SUBROUTINE scale_and_distribute
END MODULE ps_wavelet_util
