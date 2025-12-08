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

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_methods'

   PUBLIC :: init_densmat

CONTAINS
! TODO : For GPU parallelism, probably need to write a kernel instead of a subroutine and
!        handle the kpoint indices on same level as the matrix indices
! **************************************************************************************************
!> \brief Populates the density matrix
!> \param qs_env QS environment from the DFT calculation
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE init_densmat(densmat, ham_rs, smat_rs, ham, smat, kpoints, qs_env)
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), INTENT(OUT)     :: densmat
      REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(IN)         :: ham_rs
      REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(IN)         :: smat_rs
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), INTENT(INOUT)   :: ham
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), INTENT(INOUT)   :: smat
      TYPE(kpoint_type), POINTER, INTENT(IN)                :: kpoints
      TYPE(qs_environment_type), POINTER, INTENT(IN)        :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER              :: mos
      REAL(kind=dp), DIMENSION(:), POINTER                  :: occupations
      INTEGER                                               :: n_ao, &
                                                               ikp, &
                                                               ispin, &
                                                               info, &
                                                               lwork, &
                                                               mepos, &
                                                               i, &
                                                               j, &
                                                               j_local
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE              :: eigs
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE           :: work
      COMPLEX(kind=dp), DIMENSION(1)                        :: rec_work
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE              :: rwork
      COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: smat_copy, &
                                                               ham_copy
      ! Start by determining the kp matrices
      mepos = kpoints%para_env%mepos + 1
      DO j=kpoints%kp_dist(1,mepos), kpoints%kp_dist(2,mepos)
         j_local = j - kpoints%kp_dist(1, mepos) + 1
         CALL rs_to_kp(smat_rs(1,:,:,:), smat(1,:,:,j_local), kpoints%index_to_cell, kpoints%xkp(:, j))
         DO i=1,SIZE(ham_rs, 1)
            CALL rs_to_kp(ham_rs(i,:,:,:), ham(i,:,:,j_local), &
                          kpoints%index_to_cell, kpoints%xkp(:,j))
         END DO
      END DO

      CALL get_qs_env(qs_env, mos=mos)

      ! Diagonalize the hamiltonian with overlap to get the MOs
      n_ao = SIZE(densmat, 2)
      ALLOCATE (eigs(n_ao))
      ALLOCATE (rwork(3*n_ao-2))
      ALLOCATE (smat_copy(n_ao, n_ao))
      ALLOCATE (ham_copy(n_ao, n_ao))
      ! TODO : Sanity checks
      ! Query how much work is needed
      CALL ZHEGV(1, 'V', 'U', n_ao, ham(1, :, :, 1), n_ao, smat(1, :, :, 1), n_ao, &
         eigs, rec_work, -1, rwork, info)
      lwork = INT(rec_work(1))
      ALLOCATE (work(lwork))
      DO ispin=1,SIZE(ham, 1)
         ! Get the occupations from MOS
         CALL get_mo_set(mos(ispin), occupation_numbers=occupations)
         DO ikp=1,SIZE(ham, 4)
            ! Copy input matrices
            smat_copy(:,:) = smat(ispin, :, :, ikp)
            ham_copy(:,:) = ham(ispin, :, :, ikp)
            ! ZHEGV(itype, jobz, uplo, n, A, lda, B, ldb, &
            !       w, work, lwork, rwork, info)
            CALL ZHEGV(1, 'V', 'U', n_ao, ham_copy, n_ao, smat_copy, n_ao, &
                       eigs, work, lwork, rwork, info)
            ! ham_copy is now overwritten with eigenvectors
            ! Density mat (in contravariant representation) is given as C * P_occ * C^H
            ! Save copy of eigenvectors
            smat_copy(:,:) = ham_copy(:,:)
            ! Scale C by P_occ (P_occ * C)
            DO i=1,n_ao
               IF (i > SIZE(occupations)) THEN
                  ham_copy(:,i) = 0.0_dp
               ELSE
                  ! ! DEBUG : Non-equilibrium density - transitions 1<->2
                  ! IF (i == 1 .OR. i == 2) THEN
                  !    ham_copy(:, i) = 0.5_dp * ham_copy(:,i) + 0.5_dp * ham_copy(:,3-i)
                  ! ELSE
                  !    ham_copy(:,i) = ham_copy(:,i) * occupations(i)
                  ! END IF
                  ham_copy(:,i) = ham_copy(:,i) * occupations(i)
               END IF
            END DO
            ! Multiply to get the densmat
            CALL gemm_square(ham_copy, 'N', smat_copy, 'C', densmat(ispin, :, :, ikp))
         END DO
      END DO

      DEALLOCATE (eigs)
      DEALLOCATE (rwork)
      DEALLOCATE (work)

   END SUBROUTINE init_densmat
END MODULE sbe_methods
