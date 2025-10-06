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
   SUBROUTINE run_propagation_sbe()
      CPABORT("Control in SBE module!")
   END SUBROUTINE run_propagation_sbe
END MODULE sbe
