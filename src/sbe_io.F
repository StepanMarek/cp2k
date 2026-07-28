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
   USE rt_propagation_output, ONLY: print_rt_file
   USE sbe_methods, ONLY: nonherm_metric
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_io'

   CHARACTER(len=*), PARAMETER, PRIVATE :: dformat = 'E20.8E3'
   CHARACTER(len=*), PARAMETER, PRIVATE :: col_name_format = 'A20'

   PUBLIC :: output_field, &
            output_occupation_file, &
            print_mat, &
            output_nonherm_metric

CONTAINS
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param timestep current timestep
!> \param field field values
!> \param field_section input section corresponding to the field output
!> \param header_opt when true, prints the header to a file, defaults to false
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
! **************************************************************************************************
!> \brief Prints the occupation numbers - array of complex diagonal entries of the MO density
!> \param occ_numbers Occupation numbers at the given time
!> \param time Current simulation time
!> \param occupation_print_section Section specifying the filename, iterations etc.
!> \param header_opt when true, prints the header to a file, defaults to false
!> \param vacancy_opt Instead of occupations, print vacancies, i.e. 1-occupations. Extension change.
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE output_occupation_file(occ_numbers, time, occupation_print_section, &
         header_opt, vacancy_opt)
      COMPLEX(kind=dp), DIMENSION(:)                     :: occ_numbers
      REAL(kind=dp)                                      :: time
      TYPE(section_vals_type)                            :: occupation_print_section
      LOGICAL, OPTIONAL                                  :: header_opt
      LOGICAL, OPTIONAL                                  :: vacancy_opt

      LOGICAL                                            :: do_headers
      LOGICAL                                            :: do_vacancies
      CHARACTER(len=20), DIMENSION(:), ALLOCATABLE       :: headers
      TYPE(cp_logger_type), POINTER                      :: logger
      COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE      :: z_occ
      REAL(kind=dp), DIMENSION(1)                        :: a_time
      INTEGER                                            :: n_occ, i, occ_unit

      do_headers = .FALSE.
      IF (PRESENT(header_opt)) do_headers = header_opt
      do_vacancies = .FALSE.
      IF (PRESENT(vacancy_opt)) do_vacancies = vacancy_opt

      a_time(1) = time

      n_occ = SIZE(occ_numbers)
      ALLOCATE (z_occ(n_occ,1))

      logger => cp_get_default_logger()

      IF (do_vacancies) THEN
         z_occ(:,1) = CMPLX(1.0, 0.0, kind=dp)-occ_numbers(:)
         occ_unit = cp_print_key_unit_nr(logger, occupation_print_section, extension="_vacant.dat")
      ELSE
         z_occ(:,1) = occ_numbers(:)
         occ_unit = cp_print_key_unit_nr(logger, occupation_print_section, extension=".dat")
      END IF


      IF (do_headers) THEN
         ALLOCATE (headers(2*n_occ+1))

         WRITE (headers(1),'(A20)') "# Time [fs]         "
         DO i=2,2*n_occ+1
            IF (MODULO(i, 2) == 1) THEN
               WRITE (headers(i),'(A10,I8,A2)') "Im       n", (i-1)/2, "  "
            ELSE
               WRITE (headers(i),'(A10,I8,A2)') "Re       n", i/2, "  "
            END IF
         END DO

         CALL print_rt_file(occ_unit, headers, a_time, z_occ, xscale_opt=femtoseconds)

         DEALLOCATE (headers)
      ELSE
         CALL print_rt_file(occ_unit, xvals=a_time, yvals=z_occ, xscale_opt=femtoseconds)
      END IF

      DEALLOCATE (z_occ)

   END SUBROUTINE output_occupation_file
! **************************************************************************************************
!> \brief Outputs the maximum non-hermitian metric among several matrices
!> \param densmat 3D array, with two middle dimensions matrix rows and columns
!> \param unit_nr Number of the output unit (determined elsewhere)
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE output_nonherm_metric(densmat, unit_nr)
      COMPLEX(kind=dp), DIMENSION(:,:,:,:)               :: densmat
      INTEGER                                            :: unit_nr

      INTEGER                                            :: ispin, i, j, ikp_local, &
                                                            nspin, n_ao, nkp_local
      REAL(kind=dp)                                      :: buffer, max_buffer

      nspin = SIZE(densmat, 1)
      n_ao = SIZE(densmat, 2)
      nkp_local = SIZE(densmat, 4)
      buffer = 0.0_dp
      max_buffer = 0.0_dp

      ! Iterate over all dimensions
      DO ispin=1, nspin
         DO ikp_local=1, nkp_local
            buffer = nonherm_metric(densmat(ispin,:,:,ikp_local))
            IF (buffer > max_buffer) max_buffer = buffer
         END DO
      END DO
      ! Maximum found
      PRINT "(A28,E52.8E3)", " SBE| Max non-hermitian part", max_buffer
   END SUBROUTINE output_nonherm_metric
! **************************************************************************************************
!> \brief Pretty prints a complex matrix
!> \param mat Matrix to be pretty printed
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE print_mat(mat)
      COMPLEX(kind=dp), DIMENSION(:,:) :: mat

      INTEGER :: i,j,m,n

      m = SIZE(mat,1)
      n = SIZE(mat,2)

      DO i=1,m
         DO j=1,n-1
            WRITE (6,"(E15.3E3,E15.3E3)", advance="no") REAL(mat(i,j)),AIMAG(mat(i,j))
         END DO
            WRITE (6,"(E15.3E3,E15.3E3)") REAL(mat(i,j)),AIMAG(mat(i,j))
      END DO
   END SUBROUTINE print_mat
END MODULE sbe_io
