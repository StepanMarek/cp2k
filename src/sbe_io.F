!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Input and output utilities for the SBE code
!> \author Stepan Marek
!> \date 01.2026
! **************************************************************************************************
MODULE sbe_io
   USE kinds, ONLY: dp
   USE input_section_types, ONLY: section_vals_type
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output
   USE physcon, ONLY: femtoseconds
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_io'

   CHARACTER(len=*), PARAMETER, PRIVATE :: dformat = 'E20.8E3'
   CHARACTER(len=*), PARAMETER, PRIVATE :: col_name_format = 'A20'

   PUBLIC :: output_field

CONTAINS
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param timestep current timestep
!> \param field field values
!> \param field_section input section corresponding to the field output
! **************************************************************************************************
   SUBROUTINE output_field(timestep, field, field_section, header_opt)
      REAL(kind=dp)                                      :: timestep
      REAL(kind=dp), DIMENSION(3)                        :: field
      TYPE(section_vals_type), POINTER                   :: field_section
      LOGICAL, OPTIONAL                                  :: header_opt
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                            :: header
      INTEGER                                            :: field_unit
      ! TODO : Handle outputs to std_out - same as in the rtbse case?

      logger => cp_get_default_logger()
      
      header = .FALSE.
      IF (PRESENT(header_opt)) header = header_opt

      field_unit = cp_print_key_unit_nr(logger, field_section, extension=".dat")

      IF (field_unit > 0) THEN
         IF (header) THEN
            WRITE (field_unit, '('//col_name_format//','//col_name_format//','// &
            col_name_format//','//col_name_format//')') &
            "#          Time [fs]", "E_x[at.u.]", "E_y[at.u.]", "E_z[at.u.]"
         END IF
         WRITE (field_unit, '('//dformat//','//dformat//','//dformat//','//dformat//')') &
            timestep*femtoseconds, field(1), field(2), field(3)
      END IF

      CALL cp_print_key_finished_output(field_unit, logger, field_section)
   END SUBROUTINE output_field
END MODULE sbe_io
