!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains methods used in the context of density fitting
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      02.2013 moved from admm_methods
!> \author Manuel Guidon
! **************************************************************************************************
MODULE admm_utils
   USE admm_types,                      ONLY: admm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE input_constants,                 ONLY: do_admm_purify_cauchy,&
                                              do_admm_purify_cauchy_subspace,&
                                              do_admm_purify_mo_diag,&
                                              do_admm_purify_mo_no_diag,&
                                              do_admm_purify_none
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: admm_correct_for_eigenvalues, &
             admm_uncorrect_for_eigenvalues

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_utils'

!***

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param ks_matrix ...
! **************************************************************************************************
   SUBROUTINE admm_correct_for_eigenvalues(ispin, admm_env, ks_matrix)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: ks_matrix

      INTEGER                                            :: nao_aux_fit, nao_orb
      TYPE(dbcsr_type), POINTER                          :: work

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      IF (.NOT. admm_env%block_dm) THEN
         SELECT CASE (admm_env%purification_method)
         CASE (do_admm_purify_cauchy_subspace)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)

            ! ** calculate A^T*H_tilde*A
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%K(ispin), admm_env%A, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%H_corr(ispin))

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_diag)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            ! ** calculate A^T*H_tilde*A
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%K(ispin), admm_env%A, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%H_corr(ispin))

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_cauchy)
            ! do nothing
         END SELECT
      END IF

   END SUBROUTINE admm_correct_for_eigenvalues

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param ks_matrix ...
! **************************************************************************************************
   SUBROUTINE admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_matrix)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: ks_matrix

      INTEGER                                            :: nao_aux_fit, nao_orb
      TYPE(dbcsr_type), POINTER                          :: work

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      IF (.NOT. admm_env%block_dm) THEN
         SELECT CASE (admm_env%purification_method)
         CASE (do_admm_purify_cauchy_subspace)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_diag)
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_cauchy)
            ! do nothing
         END SELECT
      END IF
   END SUBROUTINE admm_uncorrect_for_eigenvalues

END MODULE admm_utils
