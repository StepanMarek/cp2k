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
   USE input_constants,                ONLY: rtp_method_sbe
   USE force_env_types,                ONLY: force_env_type, &
                                             force_env_get
   USE force_env_methods,              ONLY: force_env_calc_energy_force
   USE qs_environment_types,           ONLY: qs_environment_type, &
                                             get_qs_env
   USE qs_neighbor_list_types,         ONLY: neighbor_list_set_p_type
   USE kpoint_types,                   ONLY: kpoint_type, &
                                             get_kpoint_info, &
                                             kpoint_create, &
                                             kpoint_release
   USE kpoint_methods,                 ONLY: rskp_transform, &
                                             kpoint_init_cell_index
   USE kinds,                          ONLY: dp
   USE cp_log_handling,                ONLY: cp_logger_get_default_io_unit
   USE message_passing,                ONLY: mp_para_env_type
   USE sbe_utils,                      ONLY: print_sbe_header_info
   USE sbe_types,                      ONLY: sbe_env_type, &
                                             init_sbe_env, &
                                             release_sbe_env
   USE cp_dbcsr_api,                   ONLY: dbcsr_p_type, &
                                             dbcsr_type, &
                                             dbcsr_create, &
                                             dbcsr_release, &
                                             dbcsr_desymmetrize, &
                                             dbcsr_get_info
   USE cp_dbcsr_contrib,               ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link,             ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cell_types,                     ONLY: cell_type
   USE kpoint_k_r_trafo_simple,        ONLY: replicate_rs_matrices, &
                                             rs_to_kp
   USE cp_control_types,               ONLY: dft_control_type

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
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(len=*), PARAMETER                        :: routineN="run_propagation_sbe"
      INTEGER                                            :: i, &
                                                            j, &
                                                            unit_nr
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: xkp
      REAL(kind=dp), DIMENSION(3)                        :: single_k
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_h_kp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      TYPE(dbcsr_type)                                   :: s_k, sc_k
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_nl
      COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE      :: k_array_complex
      INTEGER                                            :: handle
      TYPE(sbe_env_type), POINTER                        :: sbe_env

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) CALL print_sbe_header_info(unit_nr)

      CALL force_env_get(force_env, qs_env=qs_env)

      ! Run the initial DFT calculation
      IF (unit_nr > 0) WRITE(unit_nr, "(A37)") " SBE| Running initial DFT calculation"
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)
      IF (unit_nr > 0) WRITE(unit_nr, "(A34)") " SBE| Initial DFT calculation done"

      CALL timeset(routineN, handle)

      ! Init the environment
      CALL init_sbe_env(sbe_env, force_env)

      ! Investigate h and S matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, &
                      matrix_ks_kp=matrix_h_kp, kpoints=kpoints)
      IF (unit_nr > 0) PRINT *, "SBE| ", "S matrix kp shape ", SHAPE(sbe_env%overlap_rs)
      IF (unit_nr > 0) PRINT *, "SBE| ", "h matrix kp shape ", SHAPE(sbe_env%ham_zero_rs)
      ! Get kpoint information
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell, cell_to_index=cell_to_index, xkp=xkp, &
              sab_nl=sab_nl)

      ! === Transformation by the original method
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
      IF (unit_nr > 0) WRITE(unit_nr, "(A52,E25.3,A3)") &
              " SBE| Approximate memory per rank for overlap matrix", REAL(sbe_env%n_spin, kind=dp) * &
              REAL(sbe_env%n_images, kind=dp) * REAL(sbe_env%n_ao, kind=dp) * REAL(sbe_env%n_ao, kind=dp) / &
              (1024.0_dp * 1024.0_dp * 1024.0_dp), "GiB"
      CALL replicate_rs_matrices(matrix_s_kp, kpoints, sbe_env%overlap_rs, sbe_env%kpoints%cell_to_index)

      ALLOCATE (k_array_complex(sbe_env%n_ao, sbe_env%n_ao), source=CMPLX(0.0, 0.0, kind=dp))
      CALL rs_to_kp(sbe_env%overlap_rs(1,:,:,:), k_array_complex, sbe_env%kpoints%index_to_cell, single_k)
      IF (unit_nr > 0) THEN
         PRINT *, "=== Real part"
         DO j=1,sbe_env%n_ao
            DO i=1,sbe_env%n_ao
               PRINT *, i, j, REAL(k_array_complex(i, j))
            END DO
         END DO
         PRINT *, "=== Imag part"
         DO j=1,sbe_env%n_ao
            DO i=1,sbe_env%n_ao
               PRINT *, i, j, AIMAG(k_array_complex(i, j))
            END DO
         END DO
      END IF

      ! Deallocations
      DEALLOCATE (k_array_complex)
      CALL release_sbe_env(sbe_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_sbe
END MODULE sbe
