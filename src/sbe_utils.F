!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Utility functions for SBE code
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
MODULE sbe_utils
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_utils'

   PUBLIC :: print_sbe_header_info

CONTAINS
! **************************************************************************************************
!> \brief Print basic information about the SBE to be run
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE print_sbe_header_info(unit_nr)
      INTEGER                                           :: unit_nr

      WRITE(unit_nr, "(A80)") " /--------------------------------------"// &
              "---------------------------------------\"
      WRITE(unit_nr, "(A2,A78)") " |", "|"
      WRITE(unit_nr, "(A2,A59,A19)") " |", "Real-time Semiconductor Bloch Equations (SBEs)", "|"
      WRITE(unit_nr, "(A2,A78)") " |", "|"
      WRITE(unit_nr, "(A80)") " \--------------------------------------"// &
              "---------------------------------------/"
   END SUBROUTINE print_sbe_header_info
END MODULE sbe_utils
