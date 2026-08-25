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
   USE kinds, ONLY: dp, &
                    default_path_length
   USE input_section_types, ONLY: section_vals_type
   USE kpoint_types, ONLY: kpoint_type
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output, &
                                 cp_print_key_generate_filename
   USE physcon, ONLY: femtoseconds
   USE rt_propagation_output, ONLY: print_rt_file, &
                                    rt_file_comp_real
   USE sbe_methods, ONLY: nonherm_metric, &
                          nonidempot_metric
   USE sbe_types, ONLY: sbe_env_type
   USE cp_files, ONLY: file_exists
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_io'

   CHARACTER(len=*), PARAMETER, PRIVATE :: dformat = 'E20.8E3'
   CHARACTER(len=*), PARAMETER, PRIVATE :: col_name_format = 'A20'

   PUBLIC :: output_field, &
            output_occupation_file, &
            output_occupations_kp_file, &
            print_mat, &
            output_nonidempot_metric, &
            gather_occupations, &
            output_ham_rs_restart, &
            read_ham_rs_restart, &
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
      IF (unit_nr > 0) WRITE (unit_nr, "(A28,E52.8E3)") " SBE| Max non-hermitian part", max_buffer
   END SUBROUTINE output_nonherm_metric
! **************************************************************************************************
!> \brief Outputs the maximum non-idempotent metric among several matrices
!> \param densmat 3D array, with two last dimensions matrix rows and columns
!> \param smat overlap matrix
!> \param work Density matrix workspace, at least 2 matrices
!> \param unit_nr Number of the output unit (determined elsewhere)
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE output_nonidempot_metric(densmat, smat, work, unit_nr)
      COMPLEX(kind=dp), DIMENSION(:,:,:)                 :: densmat
      COMPLEX(kind=dp), DIMENSION(:,:)                   :: smat
      COMPLEX(kind=dp), DIMENSION(:,:,:)                 :: work
      INTEGER                                            :: unit_nr

      INTEGER                                            :: ispin, i, j, ikp_local, &
                                                            nspin, n_ao, nkp_local
      REAL(kind=dp)                                      :: buffer, max_buffer

      nspin = SIZE(densmat, 1)
      n_ao = SIZE(densmat, 2)
      buffer = 0.0_dp
      max_buffer = 0.0_dp

      ! Iterate over all dimensions
      DO ispin=1, nspin
         buffer = nonidempot_metric(densmat(ispin,:,:), smat, work)
         IF (buffer > max_buffer) max_buffer = buffer
      END DO
      ! Maximum found
      WRITE (unit_nr, "(A29,E51.8E3)") " SBE| Max non-idempotent part", max_buffer
   END SUBROUTINE output_nonidempot_metric
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
! **************************************************************************************************
!> \brief Gathers the k-point resolved occupations from different MPI ranks
!> \param local_occupations Array of local occupations
!> \param occupations Memory allocated on IO rank only, where all the occupations are gathered
!> \param kpoints Kpoint type, which holds the distribution information
!> \param ispin Spin index, defaults to 1 (closed shell case)
!> \author Stepan Marek
!> \date 08.2026
! **************************************************************************************************
   SUBROUTINE gather_occupations(local_occupations, occupations, kpoints, ispin_opt)
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER        :: local_occupations ! spin,band,kpoint
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER        :: occupations
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, OPTIONAL                                  :: ispin_opt
      INTEGER                                            :: mepos, &
                                                            ispin, &
                                                            ikp_min, &
                                                            ikp_max
      ispin = 1
      IF (PRESENT(ispin_opt)) ispin = ispin_opt

      CALL kpoints%para_env%gather( &
         local_occupations(ispin,:,:), &
         occupations(ispin, :, :) &
      )
   END SUBROUTINE gather_occupations
! **************************************************************************************************
!> \brief Gathers the k-point resolved occupations from different MPI ranks and prints them to file
!> \param sbe_env SBE Environment (for parameters only)
!> \param local_occupations Explicitly passed local occupations, which are gathered on source rank
!> \param gathered_occupations Array of occupations at all kpoints, only valid on source rank
!> \author Stepan Marek
!> \date 08.2026
! **************************************************************************************************
   SUBROUTINE output_occupations_kp_file(sbe_env, local_occupations, gathered_occupations, headers_opt)
      TYPE(sbe_env_type), INTENT(IN)                     :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER, &
         INTENT(IN)                                      :: local_occupations
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER, &
         INTENT(OUT)                                     :: gathered_occupations
      LOGICAL, OPTIONAL, INTENT(IN)                      :: headers_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_occupations_kp_file"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, &
                                                            ispin, &
                                                            n_spin, &
                                                            ikp, &
                                                            nkp, &
                                                            iband, &
                                                            n_ao, &
                                                            occ_unit
      CHARACTER(len=20), DIMENSION(:), ALLOCATABLE       :: headers
      LOGICAL                                            :: do_headers
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: xvals

      CALL timeset(routineN, handle)

      n_spin = sbe_env%n_spin
      nkp = sbe_env%kpoints%nkp
      n_ao = sbe_env%n_ao
      do_headers = .FALSE.
      IF (PRESENT(headers_opt)) do_headers = headers_opt
      ! Do matrix gather for each spin channel
      DO ispin=1, n_spin
         CALL sbe_env%kpoints%para_env%gather( &
            local_occupations(ispin, :, :), &
            gathered_occupations(ispin, :, :) &
         )
      END DO

      ! Open the file
      logger => cp_get_default_logger()

      occ_unit = cp_print_key_unit_nr(logger, sbe_env%occ_kp_section, extension=".dat")

      IF (do_headers) THEN
         ! Prepare the headers
         ALLOCATE (headers(n_ao + 1))
         headers(1) = "#                ikp"
         DO iband=1, n_ao
            WRITE (headers(iband + 1), "(A14,I6)") &
               "     occ. band", iband
         END DO
      END IF

      ALLOCATE (xvals(nkp), source=REAL(0.0, kind=dp))
      DO ikp=1, nkp
         xvals(ikp) = REAL(ikp, kind=dp)
      END DO
      ! Print the reshaped array
      ! DEBUG : spin
      ispin = 1
      IF (do_headers) THEN
         IF (occ_unit > 0) CALL print_rt_file(occ_unit, xvals=xvals, &
            yvals=gathered_occupations(ispin, :, :), &
            headers=headers, comp_opt=rt_file_comp_real)
         DEALLOCATE (headers)
      ELSE
         IF (occ_unit > 0) CALL print_rt_file(occ_unit, xvals=xvals, &
            yvals=gathered_occupations(ispin, :, :), &
            comp_opt=rt_file_comp_real)
      END IF

      DEALLOCATE (xvals)

      CALL timestop(handle)
   END SUBROUTINE output_occupations_kp_file
! **************************************************************************************************
!> \brief Prints the RS KS Hamiltonian into unformatted matrix for restart
!> \author Stepan Marek
!> \date 08.2026
! **************************************************************************************************
   SUBROUTINE output_ham_rs_restart(restart_section, ham_rs)
      TYPE(section_vals_type)                            :: restart_section
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: ham_rs
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: n_spin, nrows, ncols, n_img, &
                                                            restart_unit, &
                                                            i_spin, irow, icol

      logger => cp_get_default_logger()
      restart_unit = cp_print_key_unit_nr(logger, restart_section, extension=".rs", &
         file_form="UNFORMATTED", file_position="REWIND")
      n_spin = SIZE(ham_rs, 1)
      nrows = SIZE(ham_rs, 2)
      ncols = SIZE(ham_rs, 3)
      n_img = SIZE(ham_rs, 4)
      ! The header contains the dimensions of the stored matrix - n_spin, nrows, ncols, n_img
      !  - present on all ranks, only io rank writes
      IF (restart_unit > 0) WRITE (restart_unit) n_spin, nrows, ncols, n_img
      ! rs matrices are also present on every rank, can output just on the io rank
      ! Note - records must be less than 2G, since record header is I4 -- TODO : check nimg?
      ! TODO : Is it worth adding a version integer?
      IF (restart_unit > 0) THEN
         DO i_spin=1,n_spin
            DO irow=1,nrows
               DO icol=1,ncols
                  WRITE (restart_unit) ham_rs(i_spin, irow, icol, :)
               END DO
            END DO
         END DO
      END IF
      CALL cp_print_key_finished_output(restart_unit, logger, restart_section)
   END SUBROUTINE output_ham_rs_restart
   SUBROUTINE read_ham_rs_restart(restart_section, ham_rs, success)
      TYPE(section_vals_type), POINTER                   :: restart_section
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: ham_rs
      LOGICAL                                            :: success
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER                                            :: n_spin, nrows, ncols, n_img, &
                                                            i_spin, irow, icol, &
                                                            target_n_spin, target_nrows, &
                                                            target_ncols, target_n_img, &
                                                            restart_unit

      logger => cp_get_default_logger()
      success = .TRUE.

      ! Check whether the file exists
      filename = cp_print_key_generate_filename(logger, restart_section, extension=".rs", &
         my_local=.FALSE.)
      IF (.NOT. file_exists(filename)) THEN
         success=.FALSE.
         RETURN
      END IF

      ! The file exists - read header and figure out whether we are compatible
      restart_unit = cp_print_key_unit_nr(logger, restart_section, extension=".rs", &
         file_form="UNFORMATTED", file_position="REWIND", file_action="READ")

      target_n_spin = SIZE(ham_rs, 1)
      target_nrows = SIZE(ham_rs, 2)
      target_ncols = SIZE(ham_rs, 3)
      target_n_img = SIZE(ham_rs, 4)
      ! Check compatibility on io rank only
      IF (restart_unit > 0) THEN
         READ (restart_unit) n_spin, nrows, ncols, n_img
      END IF
      ! Broadcast to other ranks
      CALL logger%para_env%bcast(n_spin, logger%para_env%source)
      CALL logger%para_env%bcast(nrows, logger%para_env%source)
      CALL logger%para_env%bcast(ncols, logger%para_env%source)
      CALL logger%para_env%bcast(n_img, logger%para_env%source)
      ! The comparison is done on all ranks
      IF (target_n_spin /= n_spin .OR. target_nrows /= nrows .OR. &
         target_ncols /= ncols .OR. target_n_img /= n_img) THEN
         success = .FALSE.
         RETURN
      END IF

      ! The data should be compatible, proceed with the read on io rank
      DO i_spin=1,n_spin
         DO irow=1,nrows
            DO icol=1,ncols
               IF (restart_unit > 0) READ (restart_unit) ham_rs(i_spin, irow, icol, :)
               ! Broadcast to other ranks
               CALL logger%para_env%bcast(ham_rs(i_spin, irow, icol, :), logger%para_env%source)
            END DO
         END DO
      END DO
   END SUBROUTINE read_ham_rs_restart
END MODULE sbe_io
