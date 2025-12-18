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
   USE kinds,                          ONLY: dp
   USE mathconstants,                  ONLY: twopi
   USE kpoint_types,                   ONLY: kpoint_type, &
                                             get_kpoint_info
   USE qs_environment_types,           ONLY: qs_environment_type, &
                                             get_qs_env
   USE qs_mo_types,                    ONLY: mo_set_type, &
                                             get_mo_set
   USE kpoint_k_r_trafo_simple,        ONLY: rs_to_kp
   USE mathlib,                        ONLY: gemm_square
   USE sbe_types,                      ONLY: sbe_env_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_methods'

   PUBLIC :: init_densmat, &
             get_damping_term, &
             add_occ_numbers_kp

CONTAINS
! TODO : For GPU parallelism, probably need to write a kernel instead of a subroutine and
!        handle the kpoint indices on same level as the matrix indices
! **************************************************************************************************
!> \brief Populates the density matrix
!> \param densmat Densitry matrix to be populated (specific spin and kpoint)
!> \param ispin Spin index of the given density matrix (needed to retrieve relevant occupation)
!> \param ham Kp hamiltonian
!> \param smat Kp overlap
!> \param qs_env QS environment, needed to get occupations
!> \param sbe_env SBE environment, provides workspace (densmat_work(1:2,:,:))
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE init_densmat(densmat, ispin, ham, smat, qs_env, sbe_env)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)         :: densmat
      INTEGER, INTENT(IN)                                   :: ispin
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: ham
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: smat
      TYPE(qs_environment_type), POINTER, INTENT(IN)        :: qs_env
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      TYPE(mo_set_type), DIMENSION(:), POINTER              :: mos
      REAL(kind=dp), DIMENSION(:), POINTER                  :: occupations
      INTEGER                                               :: n_ao, &
                                                               i, &
                                                               info

      CALL get_qs_env(qs_env, mos=mos)

      ! TODO : Get rid of sbe_env dependence by providing workspaces?
      ! Diagonalize the hamiltonian with overlap to get the MOs
      n_ao = sbe_env%n_ao
      ! TODO : Sanity checks
      ! Get the occupations from MOS
      CALL get_mo_set(mos(ispin), occupation_numbers=occupations)
      ! Copy input matrices
      sbe_env%densmat_work(1,:,:) = smat(:, :)
      sbe_env%densmat_work(2,:,:) = ham(:, :)
      CALL ZHEGV(1, 'V', 'U', n_ao, sbe_env%densmat_work(2,:,:), n_ao, &
                 sbe_env%densmat_work(1,:,:), n_ao, &
                 sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
      ! ham_copy is now overwritten with eigenvectors
      ! Density mat (in contravariant representation) is given as C * P_occ * C^H
      ! Save copy of eigenvectors
      sbe_env%densmat_work(1,:,:) = sbe_env%densmat_work(2,:,:)
      ! Scale C by P_occ (P_occ * C)
      DO i=1,n_ao
         ! ! DEBUG : Non-equilibrium density - transitions 1<->2
         ! IF (i == 1 .OR. i == 2) THEN
         !    sbe_env%densmat_work(2, :, i) = 0.5_dp * sbe_env%densmat_work(1,:,i) &
         !                                    + 0.5_dp * sbe_env%densmat_work(1,:,3-i)
         ! ELSE
         !    IF (i > SIZE(occupations)) THEN
         !       sbe_env%densmat_work(2,:,i) = 0.0_dp
         !    ELSE
         !       ! Set up according to occupations
         !       ! Scales occupations so that the density is still idempotent even in closed shell case
         !       sbe_env%densmat_work(2,:,i) = sbe_env%densmat_work(1,:,i) * &
         !          REAL(occupations(i), kind=dp) / REAL(sbe_env%n_spin, kind=dp)
         !    END IF
         ! END IF
         IF (i > SIZE(occupations)) THEN
            sbe_env%densmat_work(2,:,i) = 0.0_dp
         ELSE
            ! Set up according to occupations
            ! Scales occupations so that the density is still idempotent even in closed shell case
            ! TODO : Can n_spin be different than 1 or 2?
            sbe_env%densmat_work(2,:,i) = sbe_env%densmat_work(1,:,i) * &
               REAL(occupations(i), kind=dp) / REAL(3-sbe_env%n_spin, kind=dp)
         END IF
      END DO
      ! Multiply to get the densmat
      CALL gemm_square(sbe_env%densmat_work(2,:,:), 'N', &
                       sbe_env%densmat_work(1,:,:), 'C', densmat(:, :))

   END SUBROUTINE init_densmat

! **************************************************************************************************
!> \brief Determines the damping term using dephasing time, uses damping workspace
!> \param sbe_env SBE environment with current kp zero field hamiltonian and overlap
!> \param damping_term Complex matrix representing the damping term
!> \param densmat Current kp density matrix, which is needed to construct the damping term
!> \param ispin Spin index of the provided densmat (used to retreat other kpoint matrices)
!> \param ikp Kpoint index of the provided densmat (used to retreat other kpoint matrices)
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE get_damping_term(sbe_env, damping_term, densmat, ispin)
      TYPE(sbe_env_type), POINTER, INTENT(INOUT)            :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: damping_term
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)          :: densmat
      INTEGER, INTENT(IN)                                   :: ispin
      CHARACTER(len=*), PARAMETER                           :: routineN="get_damping_term"
      INTEGER                                               :: handle, &
                                                               info, &
                                                               i

      CALL timeset(routineN, handle)

      IF (sbe_env%dephasing_time == 0.0_dp) THEN
         damping_term(:,:) = 0.0_dp
         CALL timestop(handle)
         RETURN
      END IF

      ! Copy ham_kp and s_kp, in order to not lose them after diagonalization
      sbe_env%damping_work(1,:,:) = sbe_env%ham_zero_kp(ispin,:,:)
      sbe_env%damping_work(2,:,:) = sbe_env%overlap_kp(:,:)
      ! Allocate array for eigenvalues
      ! Diagonalize
      CALL ZHEGV(1, 'V', 'U', sbe_env%n_ao, sbe_env%damping_work(1,:,:), &
                 sbe_env%n_ao, sbe_env%damping_work(2,:,:), sbe_env%n_ao, &
                 sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
      ! With the diagonals, can get rho_MO
      CALL gemm_square(sbe_env%damping_work(1,:,:), 'C', &
                       sbe_env%overlap_kp(:,:), 'N', &
                       densmat, 'N', &
                       sbe_env%damping_work(2,:,:))
      CALL gemm_square(sbe_env%damping_work(2,:,:), 'N', &
                       sbe_env%overlap_kp(:,:), 'N', &
                       sbe_env%damping_work(1,:,:), 'N', &
                       damping_term)
      ! damping_term now contains rho^MO
      ! Set the diagonal to zero - no damping on the diagonal
      ! TODO : OMP Needed?
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(sbe_env,damping_term)
      DO i=1,sbe_env%n_ao
         damping_term(i,i) = 0.0_dp
      END DO
      !$OMP END PARALLEL DO

      ! To finish the damping term, divide by T_2 time, multiply by correct sign
      ! TODO : OMP?
      damping_term(:,:) = -damping_term(:,:) / sbe_env%dephasing_time

      CALL timestop(handle)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Adds to the KS occupation number, for the given kpoint
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE add_occ_numbers_kp(sbe_env, occ_numbers, densmat, ikp)
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      REAL(kind=dp), DIMENSION(:,:)                         :: occ_numbers ! dim 1 : spin, dim 2 : band
      COMPLEX(kind=dp), DIMENSION(:,:,:)                    :: densmat ! dim 1 : spin, dim 2 : rows, dim 3 : cols
      INTEGER                                               :: ikp
      INTEGER                                               :: n, info, ispin
      ! For current local kpoint, get the transformation from AO to get n_m(k)
      ! Get the molecular orbitals
      DO ispin=1, sbe_env%n_spin
         sbe_env%densmat_work(1,:,:) = sbe_env%ham_zero_kp(ispin,:,:)
         sbe_env%densmat_work(2,:,:) = sbe_env%overlap_kp(:,:)
         CALL ZHEGV(1, 'V', 'U', sbe_env%n_ao, sbe_env%densmat_work(1,:,:), &
                    sbe_env%n_ao, sbe_env%densmat_work(2,:,:), sbe_env%n_ao, &
                    sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
         ! With molecular orbitals in densmat_work(1,:,:), we can do the product
         CALL gemm_square(sbe_env%densmat_work(1,:,:), 'C', sbe_env%overlap_kp, 'N', &
                          densmat(ispin,:,:), 'N', sbe_env%densmat_work(2,:,:))
         CALL gemm_square(sbe_env%densmat_work(2,:,:), 'N', sbe_env%overlap_kp, 'N', &
                          sbe_env%densmat_work(1,:,:), 'N', sbe_env%densmat_work(3,:,:))
         ! densmat_work(3) now contains occupations on the diagonal
         DO n=1,SIZE(occ_numbers)
            ! TODO : Are occupation numbers always necceserily real?
            occ_numbers(ispin,n) = occ_numbers(ispin,n) + &
               REAL(sbe_env%densmat_work(3,n,n), kind=dp) * sbe_env%kpoints%wkp(ikp)
         END DO
      END DO
   END SUBROUTINE
END MODULE sbe_methods
