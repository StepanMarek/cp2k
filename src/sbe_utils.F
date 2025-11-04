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
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
                                              dbcsr_get_info, &
                                              dbcsr_type, &
                                              dbcsr_desymmetrize, &
                                              dbcsr_get_block_p, &
                                              dbcsr_release, &
                                              dbcsr_iterator_type, &
                                              dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, &
                                              dbcsr_iterator_blocks_left, &
                                              dbcsr_iterator_next_block
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE kpoint_types,                    ONLY: kpoint_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_utils'

   PUBLIC :: print_sbe_header_info, &
           distribute_rs_matrices

CONTAINS
! **************************************************************************************************
!> \brief Print basic information about the SBE to be run
!> \param unit_nr Output unit
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
