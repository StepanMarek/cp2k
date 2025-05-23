!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner_apply
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_release, dbcsr_type
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE input_constants,                 ONLY: ot_precond_full_all,&
                                              ot_precond_full_kinetic,&
                                              ot_precond_full_single,&
                                              ot_precond_full_single_inverse,&
                                              ot_precond_s_inverse,&
                                              ot_precond_solver_direct,&
                                              ot_precond_solver_inv_chol,&
                                              ot_precond_solver_update
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

   PUBLIC :: apply_preconditioner_fm, apply_preconditioner_dbcsr

CONTAINS

! **************************************************************************************************
!> \brief applies a previously created preconditioner to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Solver not implemented")
         END SELECT
      CASE DEFAULT
         CPABORT("Unknown preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_fm

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CPABORT("Apply_full_direct not supported with ot")
            !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Wrong solver")
         END SELECT
      CASE DEFAULT
         CPABORT("Wrong preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_dbcsr

! **************************************************************************************************
!> \brief apply to full matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_single'

      INTEGER                                            :: handle, k, n

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_out)
      CALL timestop(handle)

   END SUBROUTINE apply_full_single

! **************************************************************************************************
!> \brief apply to dbcsr matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_single'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(preconditioner_env%dbcsr_matrix)) &
         CPABORT("NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, matrix_in, &
                          0.0_dp, matrix_out)

      CALL timestop(handle)

   END SUBROUTINE apply_single

! **************************************************************************************************
!> \brief preconditioner contains the factorization, application done by
!>        solving the linear system
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_direct'

      INTEGER                                            :: handle, k, n
      TYPE(cp_fm_type)                                   :: work

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_fm_create(work, matrix_in%matrix_struct, name="apply_full_single", &
                        use_sp=matrix_in%use_sp)
      CALL cp_fm_cholesky_restore(matrix_in, k, preconditioner_env%fm, work,&
           &                      "SOLVE", transa="T")
      CALL cp_fm_cholesky_restore(work, k, preconditioner_env%fm, matrix_out,&
           &                      "SOLVE", transa="N")
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE apply_full_direct

! **************************************************************************************************
!> \brief full all to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_all'

      INTEGER                                            :: handle, i, j, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)

      CALL cp_fm_create(matrix_tmp, matrix_in%matrix_struct, name="apply_full_all")
      CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)

      !
      CALL parallel_gemm('T', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_tmp)

      ! do the right scaling
      DO j = 1, ncol_local
      DO i = 1, nrow_local
         dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                          preconditioner_env%full_evals(row_indices(i)) - preconditioner_env%occ_evals(col_indices(j)))
         local_data(i, j) = local_data(i, j)*dum
      END DO
      END DO

      ! mult back
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_tmp, 0.0_dp, matrix_out)

      CALL cp_fm_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE apply_full_all

! **************************************************************************************************
!> \brief full all to a dbcsr matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_all'

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            row, row_offset, row_size
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: DATA
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_tmp, matrix_in, name="apply_full_all")
      CALL dbcsr_multiply('T', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_in, 0.0_dp, matrix_tmp)
      ! do the right scaling
      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, DATA, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                             preconditioner_env%full_evals(row_offset + i - 1) &
                             - preconditioner_env%occ_evals(col_offset + j - 1))
            DATA(i, j) = DATA(i, j)*dum
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! mult back
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_tmp, 0.0_dp, matrix_out)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)

   END SUBROUTINE apply_all

END MODULE preconditioner_apply
