!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Main structure of the semiconductor Bloch equations (SBE) code
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
MODULE sbe
   USE input_constants,                 ONLY: rtp_method_sbe
   USE force_env_types,                 ONLY: force_env_type, &
                                              force_env_get
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE qs_environment_types,            ONLY: qs_environment_type, &
                                              get_qs_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE kpoint_types,                    ONLY: kpoint_type, &
                                              get_kpoint_info
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kinds,                           ONLY: dp
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE message_passing,                 ONLY: mp_para_env_type
   USE sbe_utils,                       ONLY: print_sbe_header_info, &
                                              distribute_rs_matrices
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
                                              dbcsr_type, &
                                              dbcsr_create, &
                                              dbcsr_release, &
                                              dbcsr_desymmetrize, &
                                              dbcsr_get_info
   USE cp_dbcsr_contrib,                ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cell_types,                      ONLY: cell_type
   USE sbe_methods,                     ONLY: integrate_rs_matrices

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe'

   PUBLIC :: run_propagation_sbe

CONTAINS
! **************************************************************************************************
!> \brief Runs the electron-only real time SBE propagation
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE run_propagation_sbe(force_env)
      TYPE(force_env_type), POINTER                     :: force_env
      TYPE(qs_environment_type), POINTER                :: qs_env
      TYPE(kpoint_type), POINTER                        :: kpoints
      CHARACTER(len=*), PARAMETER                       :: routineN="run_propagation_sbe"
      INTEGER                                           :: nkp, &
                                                           i, &
                                                           j, &
                                                           unit_nr, &
                                                           n_spin, &
                                                           n_images, &
                                                           mat_rows, &
                                                           mat_cols
      INTEGER, DIMENSION(2)                             :: dbcsr_shape
      REAL(kind=dp), DIMENSION(:, :), POINTER           :: xkp
      REAL(kind=dp), DIMENSION(3)                       :: single_k
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: matrix_s_kp, matrix_h_kp
      TYPE(cell_type), POINTER                          :: cell
      INTEGER, DIMENSION(:, :, :), POINTER              :: cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                 :: index_to_cell
      TYPE(dbcsr_type)                                  :: s_k, sc_k
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                   :: sab_nl
      INTEGER, DIMENSION(3)                             :: cells
      REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE, &
              TARGET                                    :: array
      REAL(kind=dp), DIMENSION(:,:,:,:), POINTER        :: array_p
      REAL(kind=dp), DIMENSION(:,:,:), POINTER          :: array_s_p
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, TARGET:: k_array_real
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, TARGET:: k_array_imag
      REAL(kind=dp), DIMENSION(:,:), POINTER            :: k_array_real_p
      REAL(kind=dp), DIMENSION(:,:), POINTER            :: k_array_imag_p
      TYPE(mp_para_env_type), POINTER                   :: para_env

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) CALL print_sbe_header_info(unit_nr)

      CALL force_env_get(force_env, qs_env=qs_env)

      ! Run the initial DFT calculation
      IF (unit_nr > 0) WRITE(unit_nr, "(A37)") " SBE| Running initial DFT calculation"
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)
      IF (unit_nr > 0) WRITE(unit_nr, "(A34)") " SBE| Initial DFT calculation done"


      ! Investigate h and S matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_h_kp, cell=cell, kpoints=kpoints)
      IF (unit_nr > 0) PRINT *, "SBE| ", "S matrix kp shape ", SHAPE(matrix_s_kp)
      IF (unit_nr > 0) PRINT *, "SBE| ", "h matrix kp shape ", SHAPE(matrix_h_kp)
      ! Get kpoint information
      ! NOTE : cell_to_index has curious index structure - includes negative indices
      !  - think about how to print it/iterate over it?
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell, cell_to_index=cell_to_index, xkp=xkp)
      ! === Transformation by the original method
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell, cell_to_index=cell_to_index, xkp=xkp, &
              sab_nl=sab_nl)
      PRINT *, "Checking matrix inverse", index_to_cell(:, 2), " vs ", index_to_cell(:, 3)
      CALL dbcsr_print(matrix_s_kp(1, 2)%matrix)
      CALL dbcsr_print(matrix_s_kp(1, 3)%matrix)
      CALL dbcsr_create(s_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      CALL dbcsr_create(sc_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(s_k, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(sc_k, sab_nl)
      single_k(:) = xkp(:,1)
      CALL rskp_transform(rmatrix=s_k, cmatrix=sc_k, rsmat=matrix_s_kp, ispin=1, xkp=single_k, &
              cell_to_index=cell_to_index, sab_nl=sab_nl)
      CALL dbcsr_print(s_k)
      CALL dbcsr_print(sc_k)
      CALL dbcsr_release(s_k)
      CALL dbcsr_release(sc_k)
      ! === Transformation by the duplicated method
      dbcsr_shape = SHAPE(matrix_s_kp)
      n_spin = dbcsr_shape(1)
      n_images = dbcsr_shape(2)
      CALL dbcsr_get_info(matrix_s_kp(1,1)%matrix, nfullrows_total=mat_rows, nfullcols_total=mat_cols)
      ALLOCATE (array(n_spin, n_images, mat_rows, mat_cols), source=0.0_dp)
      ! TODO : Update for the number of matrices stored per each rank
      IF (unit_nr > 0) WRITE(unit_nr, "(A52,E25.3,A3)") &
              " SBE| Approximate memory per rank for overlap matrix", REAL(n_spin, kind=dp) * &
              REAL(n_images, kind=dp) * REAL(mat_rows, kind=dp) * REAL(mat_cols, kind=dp) / &
              (1024.0_dp * 1024.0_dp * 1024.0_dp), "GiB"
      array_p => array
      CALL distribute_rs_matrices(matrix_s_kp, array_p, index_to_cell, cell_to_index)
      ! IF (unit_nr > 0) THEN
      !    DO i=1,mat_rows
      !       DO j=1,mat_cols
      !          PRINT *, i, j, array(1,2,i, j)
      !       END DO
      !    END DO
      ! END IF
      ! IF (unit_nr > 0) THEN
      !    DO i=1,mat_rows
      !       DO j=1,mat_cols
      !          PRINT *, i, j, array(1,3,i, j)
      !       END DO
      !    END DO
      ! END IF
      ALLOCATE (k_array_real(mat_rows, mat_cols), source=0.0_dp)
      ALLOCATE (k_array_imag(mat_rows, mat_cols), source=0.0_dp)
      k_array_real_p => k_array_real
      k_array_imag_p => k_array_imag
      array_s_p => array(1,:,:,:)
      CALL integrate_rs_matrices(k_array_real_p, k_array_imag_p, 1, array_s_p, kpoints, index_to_cell)
      IF (unit_nr > 0) THEN
         PRINT *, "=== Imag part"
         DO j=1,mat_cols
            DO i=1,mat_rows
               PRINT *, i, j, k_array_imag(i, j)
            END DO
         END DO
         PRINT *, "=== Real part"
         DO j=1,mat_cols
            DO i=1,mat_rows
               PRINT *, i, j, k_array_real(i, j)
            END DO
         END DO
      END IF
      DEALLOCATE(array)
      DEALLOCATE(k_array_real)
      DEALLOCATE(k_array_imag)
   END SUBROUTINE run_propagation_sbe
END MODULE sbe
