!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods for the SBE code (more model based, as opposed to sbe_utils, which are code based)
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
MODULE sbe_methods
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE kpoint_types,                    ONLY: kpoint_type, &
                                              get_kpoint_info

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_methods'

   PUBLIC :: integrate_rs_matrices

CONTAINS
! TODO : For GPU parallelism, probably need to write a kernel instead of a subroutine and
!        handle the kpoint indices on same level as the matrix indices
! **************************************************************************************************
!> \brief Integrate (sum) the RS-matrices with a given value of wavevector to get k-space matrices
!> \note Does not handle spin index - given arrays for a specific index only
!> \param k_array_real Matrix in the kspace - result, real part
!> \param k_array_real Matrix in the kspace - result, imaginary part
!> \param kpoint_index Specific k-point index
!> \param rs_array Rank local array of real-space matrices
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE integrate_rs_matrices(k_array_real, k_array_imag, kpoint_index, rs_array, kpoints, index_to_cell)
      REAL(kind=dp), DIMENSION(:,:), POINTER            :: k_array_real
      REAL(kind=dp), DIMENSION(:,:), POINTER            :: k_array_imag
      INTEGER                                           :: kpoint_index
      ! dimension 1 : image cell index, dimension 2 : rows, dimension 3 : cols
      REAL(kind=dp), DIMENSION(:,:,:), POINTER          :: rs_array
      TYPE(kpoint_type), POINTER                        :: kpoints
      INTEGER, DIMENSION(:, :), POINTER                 :: index_to_cell
      CHARACTER(len=*), PARAMETER                       :: routineN="integrate_rs_matrices"
      INTEGER                                           :: handle, &
                                                           i, j, k, &
                                                           n_images, &
                                                           row_size, &
                                                           col_size
      INTEGER, DIMENSION(2)                             :: mat_size
      REAL(kind=dp), DIMENSION(3)                       :: kp
      REAL(kind=dp), DIMENSION(:, :), POINTER           :: xkp
      REAL(kind=dp)                                     :: cosfactor, &
                                                           sinfactor

      CALL timeset(routineN, handle)
      ! Get matrix constants
      ! mat_size = SHAPE(k_array_real)
      ! row_size = mat_size(1)
      ! col_size = mat_size(2)
      row_size = SIZE(k_array_real, 1)
      col_size = SIZE(k_array_real, 2)
      n_images = SIZE(rs_array, 1)
      ! Get the k-point coordinates (in rec. lattice vector units)
      CALL get_kpoint_info(kpoints, xkp=xkp)
      DO i=1,3
         kp(i) = xkp(i, kpoint_index)
      END DO
      ! Now, iterate over realspace and build the sum
      ! TODO : OMP
      k_array_real(:,:) = 0.0_dp
      k_array_imag(:,:) = 0.0_dp
      DO i=1,n_images
         ! Construct the phase and factors
         cosfactor = COS(twopi*(kp(1) * index_to_cell(1,i) + &
                 kp(2) * index_to_cell(2,i) + kp(3) * index_to_cell(3,i)))
         sinfactor = SIN(twopi*(kp(1) * index_to_cell(1,i) + &
                 kp(2) * index_to_cell(2,i) + kp(3) * index_to_cell(3,i)))
         DO j=1,row_size
            DO k=1,col_size
               k_array_real(j,k) = k_array_real(j,k) + cosfactor * rs_array(i, j, k)
               k_array_imag(j,k) = k_array_imag(j,k) + sinfactor * rs_array(i, j, k)
            END DO
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE integrate_rs_matrices
END MODULE sbe_methods
