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

   USE qs_environment_types,           ONLY: qs_environment_type, &
                                             get_qs_env
   USE qs_neighbor_list_types,         ONLY: neighbor_list_set_p_type, &
                                             neighbor_list_iterator_p_type, &
                                             neighbor_list_iterator_create, &
                                             neighbor_list_iterator_release, &
                                             neighbor_list_iterate, &
                                             get_iterator_info

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_utils'

   PUBLIC :: print_sbe_header_info, &
      iterate_nl

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

   SUBROUTINE iterate_nl(sab)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iterator
      REAL(kind=dp), DIMENSION(3)                        :: r
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER                                            :: iatom, jatom

      CALL neighbor_list_iterator_create(iterator, sab)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell, r=r, iatom=iatom, jatom=jatom)
         PRINT *, "iatom", iatom, &
            "jatom", jatom, &
            "r", r, &
            "cell", cell
      END DO
      CALL neighbor_list_iterator_release(iterator)
   END SUBROUTINE iterate_nl
END MODULE sbe_utils
