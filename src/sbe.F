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
                                              get_kpoint_info, &
                                              kpoint_create, &
                                              kpoint_release
   USE kpoint_methods,                  ONLY: rskp_transform, &
                                              kpoint_init_cell_index
   USE kinds,                           ONLY: dp
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE message_passing,                 ONLY: mp_para_env_type, &
                                              mp_para_env_type
   USE sbe_utils,                       ONLY: print_sbe_header_info, &
                                              iterate_nl
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
                                              dbcsr_type, &
                                              dbcsr_create, &
                                              dbcsr_release, &
                                              dbcsr_desymmetrize, &
                                              dbcsr_get_info
   USE cp_dbcsr_contrib,                ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cell_types,                      ONLY: cell_type
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices, &
                                              rs_to_kp
   USE cp_control_types,                  ONLY: dft_control_type

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
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints, &
                                                            kpoints_all
      CHARACTER(len=*), PARAMETER                        :: routineN="run_propagation_sbe"
      INTEGER                                            :: i, &
                                                            handle, &
                                                            j, &
                                                            unit_nr, &
                                                            n_spin, &
                                                            n_images, &
                                                            mat_rows, &
                                                            mat_cols
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: xkp
      REAL(kind=dp), DIMENSION(3)                        :: single_k
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_h_kp
      TYPE(dbcsr_type)                                   :: s_k, sc_k
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_all
      REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: overlap_repl
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: k_array_real
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: k_array_imag
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) CALL print_sbe_header_info(unit_nr)

      CALL force_env_get(force_env, qs_env=qs_env)

      ! Run the initial DFT calculation
      IF (unit_nr > 0) WRITE(unit_nr, "(A37)") " SBE| Running initial DFT calculation"
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)
      IF (unit_nr > 0) WRITE(unit_nr, "(A34)") " SBE| Initial DFT calculation done"

      ! Actual propagation
      CALL timeset(routineN, handle)

      ! Create full kpoints environment
      ! TODO : Also populate with user defined grid, without symmetry (due to breaking by nonlinearity)
      NULLIFY (kpoints_all)
      CALL kpoint_create(kpoints_all)
      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env, sab_all=sab_all)
      ! NOTE : This overwrites n_images in dft_control - not consistent with what is stored in qs_env
      CALL kpoint_init_cell_index(kpoints_all, sab_all, para_env, dft_control)

      ! DEBUG : Show the contents of neighbor list
      ! CALL iterate_nl(sab_all)

      ! Replicate h and S matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_h_kp, kpoints=kpoints)
      IF (unit_nr > 0) PRINT *, "SBE| ", "S matrix kp shape ", SHAPE(matrix_s_kp)
      IF (unit_nr > 0) PRINT *, "SBE| ", "h matrix kp shape ", SHAPE(matrix_h_kp)
      ! DEBUG : Comparing with original k-transform
      CALL dbcsr_create(s_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      CALL dbcsr_create(sc_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(s_k, kpoints%sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(sc_k, kpoints%sab_nl)
      CALL get_kpoint_info(kpoints, xkp=xkp)
      single_k(:) = xkp(:,1)
      CALL rskp_transform(rmatrix=s_k, cmatrix=sc_k, rsmat=matrix_s_kp, ispin=1, xkp=single_k, &
              cell_to_index=kpoints%cell_to_index, sab_nl=kpoints%sab_nl)
      CALL dbcsr_print(sc_k)
      CALL dbcsr_print(s_k)
      CALL dbcsr_release(sc_k)
      CALL dbcsr_release(s_k)

      n_spin = SIZE(matrix_s_kp, 1)
      n_images = SIZE(kpoints_all%index_to_cell, 2)
      IF (unit_nr > 0) PRINT *, "SBE| ", "kpoints_all images ", n_images
      CALL dbcsr_get_info(matrix_s_kp(1,1)%matrix, nfullrows_total=mat_rows, nfullcols_total=mat_cols)
      ALLOCATE (overlap_repl(n_spin, n_images, mat_rows, mat_cols), source=0.0_dp)
      ! TODO : Update for the number of matrices stored per each rank
      IF (unit_nr > 0) WRITE(unit_nr, "(A52,E25.3,A3)") &
              " SBE| Approximate memory per rank for overlap matrix", REAL(n_spin, kind=dp) * &
              REAL(n_images, kind=dp) * REAL(mat_rows, kind=dp) * REAL(mat_cols, kind=dp) / &
              (1024.0_dp * 1024.0_dp * 1024.0_dp), "GiB"
      CALL replicate_rs_matrices(matrix_s_kp, kpoints, overlap_repl, kpoints_all%cell_to_index)
      ! DEBUG : Explicit contents of RS matrices
      ! PRINT *, "RS matrix contents"
      ! PRINT *, "R-vector", index_to_cell(1,1), index_to_cell(2,1), index_to_cell(3,1)
      ! CALL dbcsr_print(matrix_s_kp(1,1)%matrix)
      ! IF (unit_nr > 0) THEN
      !    DO i=1,mat_rows
      !       DO j=1,mat_cols
      !          PRINT *, i, j, array(1,1,i, j)
      !       END DO
      !    END DO
      ! END IF
      ! PRINT *, "R-vector", index_to_cell(1,2), index_to_cell(2,2), index_to_cell(3,2)
      ! CALL dbcsr_print(matrix_s_kp(1,2)%matrix)
      ! IF (unit_nr > 0) THEN
      !    DO i=1,mat_rows
      !       DO j=1,mat_cols
      !          PRINT *, i, j, array(1,2,i, j)
      !       END DO
      !    END DO
      ! END IF
      ! PRINT *, "R-vector", index_to_cell(1,3), index_to_cell(2,3), index_to_cell(3,3)
      ! CALL dbcsr_print(matrix_s_kp(1,3)%matrix)
      ! IF (unit_nr > 0) THEN
      !    DO i=1,mat_rows
      !       DO j=1,mat_cols
      !          PRINT *, i, j, array(1,3,i, j)
      !       END DO
      !    END DO
      ! END IF
      ALLOCATE (k_array_real(mat_rows, mat_cols), source=0.0_dp)
      ALLOCATE (k_array_imag(mat_rows, mat_cols), source=0.0_dp)
      PRINT *, "kpoint", single_k
      CALL rs_to_kp(k_array_real, k_array_imag, single_k, overlap_repl(1,:,:,:), kpoints_all%index_to_cell)
      IF (unit_nr > 0) THEN
         PRINT *, "=== Imag part"
         DO i=1,mat_rows
            DO j=1,mat_cols
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
      DEALLOCATE (overlap_repl)
      DEALLOCATE (k_array_real)
      DEALLOCATE (k_array_imag)
      CALL kpoint_release(kpoints_all)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_sbe
END MODULE sbe
