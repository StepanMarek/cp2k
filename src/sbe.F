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
   USE kpoint_types,                    ONLY: kpoint_type, &
                                              get_kpoint_info
   USE kinds,                           ONLY: dp
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE message_passing,                 ONLY: mp_para_env_type
   USE sbe_utils,                       ONLY: print_sbe_header_info
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cell_types,                      ONLY: cell_type

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
                                                           unit_nr
      REAL(kind=dp), DIMENSION(:, :), POINTER           :: xkp
      TYPE(mp_para_env_type), POINTER                   :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: matrix_s_kp, matrix_h_kp
      TYPE(cell_type), POINTER                          :: cell

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) CALL print_sbe_header_info(unit_nr)

      CALL force_env_get(force_env, qs_env=qs_env)
      ! CALL get_qs_env(qs_env, kpoints=kpoints, para_env=para_env)
      ! CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp)
      ! IF (unit_nr > 0) PRINT *, SHAPE(xkp)
      ! DO i=1,nkp
      !    PRINT *, i, xkp(1,i), xkp(2, i), xkp(3, i)
      ! END DO

      ! Run the initial DFT calculation
      IF (unit_nr > 0) WRITE(unit_nr, "(A37)") " SBE| Running initial DFT calculation"
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)
      IF (unit_nr > 0) WRITE(unit_nr, "(A34)") " SBE| Initial DFT calculation done"

      ! Investigate h and S matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_h_kp, cell=cell)
      IF (unit_nr > 0) PRINT *, "SBE| ", "S matrix kp shape ", SHAPE(matrix_s_kp)
      IF (unit_nr > 0) PRINT *, "SBE| ", "h matrix kp shape ", SHAPE(matrix_h_kp)
      IF (unit_nr > 0) PRINT *, "SBE| ", "cell perd ", cell%perd(1), cell%perd(2), cell%perd(3)
   END SUBROUTINE run_propagation_sbe
END MODULE sbe
