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
                              cp_logger_get_default_io_unit, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output, &
                                 cp_print_key_should_output, &
                                 cp_print_key_generate_filename
   USE physcon, ONLY: femtoseconds, &
                      evolt
   USE rt_propagation_output, ONLY: print_rt_file, &
                                    rt_file_comp_real
   USE sbe_methods, ONLY: nonherm_metric, &
                          nonidempot_metric
   USE sbe_utils, ONLY: polymorphic_gather
   USE sbe_types, ONLY: sbe_env_type
   USE cp_files, ONLY: file_exists
   USE kpoint_k_r_trafo_simple, ONLY: rs_to_kp
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
            gather_occupations, &
            output_ham_rs_restart, &
            read_ham_rs_restart, &
            output_run_info, &
            print_bandstructure

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
   SUBROUTINE output_occupation_file(occ_numbers, time, occupation_print_section, integrand_dist, &
         header_opt, vacancy_opt)
      COMPLEX(kind=dp), DIMENSION(:)                     :: occ_numbers
      REAL(kind=dp)                                      :: time
      TYPE(section_vals_type)                            :: occupation_print_section
      INTEGER, DIMENSION(:,:), INTENT(IN)                :: integrand_dist
      LOGICAL, OPTIONAL                                  :: header_opt
      LOGICAL, OPTIONAL                                  :: vacancy_opt

      LOGICAL                                            :: do_headers, &
                                                            do_vacancies, &
                                                            do_prefix
      CHARACTER(len=20), DIMENSION(:), ALLOCATABLE       :: headers
      CHARACTER(len=21)                                  :: prefix
      CHARACTER(len=5)                                   :: prefix_format
      TYPE(cp_logger_type), POINTER                      :: logger
      COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE      :: z_occ
      REAL(kind=dp), DIMENSION(1)                        :: a_time
      INTEGER                                            :: n_occ, &
                                                            i, &
                                                            occ_unit, &
                                                            integrand_nkp

      do_headers = .FALSE.
      IF (PRESENT(header_opt)) do_headers = header_opt
      do_vacancies = .FALSE.
      IF (PRESENT(vacancy_opt)) do_vacancies = vacancy_opt

      logger => cp_get_default_logger()

      ! TODO : Do determination of validity of output in a different way?
      integrand_nkp = 0
      DO i=1, SIZE(integrand_dist, 2)
         integrand_nkp = integrand_nkp + integrand_dist(2,i) - integrand_dist(1,i) + 1
      END DO

      IF (integrand_nkp <= 0) THEN
         IF (cp_print_key_should_output(logger%iter_info, occupation_print_section) /= 0) THEN
            ! Warn the user
            CPWARN("Printing of integrated observable (occupation) not possible -- no integrable k-grid.")
         END IF
         RETURN
      END IF

      a_time(1) = time

      n_occ = SIZE(occ_numbers)
      ALLOCATE (z_occ(n_occ,1))

      IF (do_vacancies) THEN
         z_occ(:,1) = CMPLX(1.0, 0.0, kind=dp)-occ_numbers(:)
         occ_unit = cp_print_key_unit_nr(logger, occupation_print_section, extension="_vacant.dat")
      ELSE
         z_occ(:,1) = occ_numbers(:)
         occ_unit = cp_print_key_unit_nr(logger, occupation_print_section, extension=".dat")
      END IF

      ! Determine, whether prefix is necessary
      prefix=""
      prefix_format="(A9)"
      do_prefix = .FALSE.
      IF (occ_unit == cp_logger_get_default_io_unit(logger)) THEN
         do_prefix = .TRUE.
         prefix=" SBE_OCC|"
         IF (do_vacancies) THEN
            prefix=" SBE_VAC|"
         END IF
      END IF

      IF (do_headers) THEN
         ALLOCATE (headers(2*n_occ+1))

         WRITE (headers(1),'(A20)') "# Time [fs]         "
         DO i=2,2*n_occ+1
            IF (MODULO(i, 2) == 1) THEN
               WRITE (headers(i),'(A12,I8)') "  Im       n", (i-1)/2
            ELSE
               WRITE (headers(i),'(A12,I8)') "  Re       n", i/2
            END IF
         END DO

         IF (do_prefix) THEN
            CALL print_rt_file(occ_unit, headers, a_time, z_occ, xscale_opt=femtoseconds, &
               prefix=prefix, prefix_format=prefix_format)
         ELSE
            CALL print_rt_file(occ_unit, headers, a_time, z_occ, xscale_opt=femtoseconds)
         END IF

         DEALLOCATE (headers)
      ELSE
         IF (do_prefix) THEN
            CALL print_rt_file(occ_unit, xvals=a_time, yvals=z_occ, xscale_opt=femtoseconds, &
               prefix=prefix, prefix_format=prefix_format)
         ELSE
            CALL print_rt_file(occ_unit, xvals=a_time, yvals=z_occ, xscale_opt=femtoseconds)
         END IF
      END IF

      DEALLOCATE (z_occ)

   END SUBROUTINE output_occupation_file
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
   SUBROUTINE output_occupations_kp_file(sbe_env, local_occupations, gathered_occupations, headers_opt, &
         vacancy_opt)
      TYPE(sbe_env_type), INTENT(IN)                     :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER, &
         INTENT(IN)                                      :: local_occupations
      COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER, &
         INTENT(OUT)                                     :: gathered_occupations
      LOGICAL, OPTIONAL, INTENT(IN)                      :: headers_opt, &
                                                            vacancy_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_occupations_kp_file"
      CHARACTER(len=5), PARAMETER                        :: prefix_format = "(A12)"
      CHARACTER(len=21)                                  :: prefix
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
      LOGICAL                                            :: do_headers, &
                                                            do_vacancies, &
                                                            do_prefix
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: xvals
      COMPLEX(kind=dp), DIMENSION(:,:,:), ALLOCATABLE    :: yvals
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: kset_length_dist

      CALL timeset(routineN, handle)

      n_spin = sbe_env%n_spin
      nkp = sbe_env%kset_nkp
      n_ao = sbe_env%n_ao
      do_headers = .FALSE.
      IF (PRESENT(headers_opt)) do_headers = headers_opt
      do_vacancies = .FALSE.
      IF (PRESENT(vacancy_opt)) do_vacancies = vacancy_opt

      ! Check whether we should output or not
      logger => cp_get_default_logger()
      IF (cp_print_key_should_output(logger%iter_info, sbe_env%occ_kp_section) == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! Exit and warn when no kset is provided
      IF (SIZE(gathered_occupations) == 0) THEN
         CPWARN("Cannot provide k-resolved observable -- no KPOINT_SET.")
         CALL timestop(handle)
         RETURN
      END IF

      ! Determine the length distribution
      ALLOCATE (kset_length_dist(SIZE(sbe_env%kset_dist, 2)), source=0)
      DO ikp=1, SIZE(kset_length_dist, 1)
         kset_length_dist(ikp) = sbe_env%kset_dist(2,ikp) - sbe_env%kset_dist(1,ikp) + 1
      END DO

      ! Carry out the gather
      DO ispin=1, n_spin
         DO iband=1, n_ao
            CALL polymorphic_gather(local_occupations(:, ispin, iband), &
               gathered_occupations(:, ispin, iband), &
               sbe_env%kpoints%para_env, &
               kset_length_dist)
         END DO
      END DO

      ! Open the file
      IF (do_vacancies) THEN
         occ_unit = cp_print_key_unit_nr(logger, sbe_env%occ_kp_section, extension="_vacant.dat")
         IF (occ_unit > 0) THEN
            gathered_occupations(:,:,:) = 1.0_dp - gathered_occupations(:,:,:)
         END IF
      ELSE
         occ_unit = cp_print_key_unit_nr(logger, sbe_env%occ_kp_section, extension=".dat")
      END IF

      IF (do_headers) THEN
         ! Prepare the headers
         ALLOCATE (headers(n_ao + 1))
         headers(1) = "#                ikp"
         DO iband=1, n_ao
            WRITE (headers(iband + 1), "(A14,I6)") &
               "     occ. band", iband
         END DO
      END IF

      do_prefix = .FALSE.
      prefix = ""
      IF (occ_unit == cp_logger_get_default_io_unit(logger)) THEN
         do_prefix = .TRUE.
         IF (do_vacancies) THEN
            prefix = " SBE_VAC_KP|"
         ELSE
            prefix = " SBE_OCC_KP|"
         END IF
      END IF

      ALLOCATE (xvals(nkp), source=REAL(0.0, kind=dp))
      ALLOCATE (yvals(n_ao, nkp, n_spin), source=CMPLX(0.0, 0.0, kind=dp))
      DO ikp=1, nkp
         xvals(ikp) = REAL(ikp, kind=dp)
         DO iband=1, n_ao
            DO ispin=1, n_spin
               yvals(iband, ikp, ispin) = gathered_occupations(ikp, ispin, iband)
            END DO
         END DO
      END DO
      ! Print the reshaped array
      ! DEBUG : spin
      ispin = 1
      IF (do_headers) THEN
         IF (do_prefix) THEN
            CALL print_rt_file(occ_unit, xvals=xvals, &
               yvals=yvals(:,:,ispin), &
               headers=headers, comp_opt=rt_file_comp_real, &
               prefix=prefix, prefix_format=prefix_format)
         ELSE
            CALL print_rt_file(occ_unit, xvals=xvals, &
               yvals=yvals(:,:,ispin), &
               headers=headers, comp_opt=rt_file_comp_real)
         END IF
         DEALLOCATE (headers)
      ELSE
         IF (do_prefix) THEN
            CALL print_rt_file(occ_unit, xvals=xvals, &
               yvals=yvals(:,:,ispin), &
               comp_opt=rt_file_comp_real, &
               prefix=prefix, prefix_format=prefix_format)
         ELSE
            CALL print_rt_file(occ_unit, xvals=xvals, &
               yvals=yvals(:,:,ispin), &
               comp_opt=rt_file_comp_real)
         END IF
      END IF

      DEALLOCATE (xvals)
      DEALLOCATE (yvals)
      DEALLOCATE (kset_length_dist)

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
   SUBROUTINE read_ham_rs_restart(restart_section, ham_rs, error)
      TYPE(section_vals_type), POINTER                   :: restart_section
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: ham_rs
      INTEGER                                            :: error
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER                                            :: n_spin, nrows, ncols, n_img, &
                                                            i_spin, irow, icol, &
                                                            target_n_spin, target_nrows, &
                                                            target_ncols, target_n_img, &
                                                            restart_unit

      logger => cp_get_default_logger()
      error = 0

      ! Check whether the file exists
      filename = cp_print_key_generate_filename(logger, restart_section, extension=".rs", &
         my_local=.FALSE.)
      IF (.NOT. file_exists(filename)) THEN
         error=1
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
         error = 2
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

! **************************************************************************************************
!> \brief Prints the initial bandstructure before propagation
!> \param ham_zero_kp The equilibrium Hamiltonian, evaluated on the kpoint grid for SBE
!> \author Stepan Marek
!> \date 09.2026
! **************************************************************************************************
   SUBROUTINE print_bandstructure(ham_zero_rs, overlap_rs, bandstructure_section, kpoints, &
         xkp, lwork, eigwork, eigrwork)
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: ham_zero_rs, &
                                                            overlap_rs
      TYPE(section_vals_type), POINTER                   :: bandstructure_section
      TYPE(kpoint_type)                                  :: kpoints
      REAL(kind=dp), DIMENSION(:,:)                      :: xkp
      INTEGER                                            :: lwork
      COMPLEX(kind=dp), DIMENSION(:)                     :: eigwork
      REAL(kind=dp), DIMENSION(:)                        :: eigrwork
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: unit_nr, &
                                                            i, j, k, &
                                                            info, &
                                                            local_nkp, &
                                                            i_local, &
                                                            nkp, &
                                                            n_ao, &
                                                            n_spin, &
                                                            mepos
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: local_eigvals
      REAL(kind=dp), DIMENSION(:,:,:), POINTER           :: gathered_eigvals
      COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE      :: h_kp, s_kp, s_work

      ! Create the logger to check whether the output is on
      logger => cp_get_default_logger()
      IF (cp_print_key_should_output(logger%iter_info, bandstructure_section) == 0) RETURN

      n_spin = SIZE(ham_zero_rs, 1)
      n_ao = SIZE(ham_zero_rs, 2)
      mepos = kpoints%para_env%mepos + 1
      local_nkp = kpoints%kp_dist(2, mepos) - kpoints%kp_dist(1,mepos) + 1
      nkp = kpoints%nkp

      ! Allocate local storage for eigenvalues
      ALLOCATE (local_eigvals(n_spin, n_ao, local_nkp))

      ! Temporary matrices for ZHEGV
      ALLOCATE (h_kp(n_ao, n_ao))
      ALLOCATE (s_kp(n_ao, n_ao))
      ALLOCATE (s_work(n_ao, n_ao))

      ! Iterate over local k-points
      DO i = kpoints%kp_dist(1, mepos), kpoints%kp_dist(2, mepos)
         i_local = i - kpoints%kp_dist(1, mepos) + 1

         ! Overlap matrix for this k-point
         CALL rs_to_kp(overlap_rs(1,:,:,:), s_kp, &
                       kpoints%index_to_cell, kpoints%xkp(:,i))

         DO j = 1, n_spin
            ! Hamiltonian matrix for this k-point and spin
            CALL rs_to_kp(ham_zero_rs(j,:,:,:), h_kp, &
                          kpoints%index_to_cell, kpoints%xkp(:,i))
            s_work(:,:) = s_kp(:,:)

            CALL ZHEGV(1, 'V', 'U', n_ao, h_kp, n_ao, &
                       s_work, n_ao, &
                       local_eigvals(j,:,i_local), eigwork, lwork, eigrwork, info)

            IF (info /= 0) THEN
               CPABORT("SBE| Bandstructure: ZHEGV failed for ikp=" // CHAR(i) // " ispin=" // CHAR(j))
            END IF
         END DO
      END DO

      ! Open the file to check if we are on IO rank
      unit_nr = cp_print_key_unit_nr(logger, bandstructure_section, extension=".dat", file_position="REWIND")

      ! Allocate gathered array only on IO rank
      NULLIFY (gathered_eigvals)
      IF (unit_nr > 0) THEN
         ALLOCATE (gathered_eigvals(n_spin, n_ao, nkp))
      ELSE
         ALLOCATE (gathered_eigvals(0,0,0))
      END IF

      ! Gather eigenvalues to the IO rank
      ! We gather spin by spin
      DO j = 1, n_spin
         CALL kpoints%para_env%gather( &
            local_eigvals(j,:,:), &
            gathered_eigvals(j,:,:) &
         )
      END DO

      ! Output on IO rank
      ! TODO : Prefix for STDOUT
      ! TODO : Check ordering
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(A10,A14,A14,A14,A5,A11,A20)") "# kp index", "k_x [b_x]", "k_y [b_y]", "k_z [b_z]", &
            "spin", "band index", "energy [eV]"

         DO i = 1, nkp
            DO j = 1, n_spin
               DO k = 1, n_ao
                  WRITE (unit_nr, "(I10,E14.6E2,E14.6E2,E14.6E2,I5,I11,E20.8E3)") &
                     i, xkp(:,i), j, k, gathered_eigvals(j, k, i)*evolt
               END DO
            END DO
         END DO
         CALL cp_print_key_finished_output(unit_nr, logger, bandstructure_section)
      END IF

      DEALLOCATE (local_eigvals)
      DEALLOCATE (h_kp)
      DEALLOCATE (s_kp)
      DEALLOCATE (s_work)
      IF (ASSOCIATED(gathered_eigvals)) DEALLOCATE (gathered_eigvals)
   END SUBROUTINE print_bandstructure
! **************************************************************************************************
!> \brief Prints various updated about the current loop
!> \param sbe_env SBE environment
!> \param max_idempot Maximum idempotency deviation within the kpoints controlled by this rank
!> \note The idempotency is given as argument since it requires S matrix, which is only kept on a
!>       per-kpoint basis. Hence, the maximum deviation is calculated during the propagation step,
!>       in order to not recalculate S matrix at each kpoint
!> \author Stepan Marek
!> \date 09.2026
! **************************************************************************************************
   SUBROUTINE output_run_info(sbe_env, max_idempot)
      TYPE(sbe_env_type), INTENT(IN)                     :: sbe_env
      REAL(kind=dp), INTENT(INOUT)                       :: max_idempot
      CHARACTER(len=*), PARAMETER                        :: routineN="output_run_info"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: ispin, &
                                                            ikp_local, &
                                                            unit_nr, &
                                                            handle
      REAL(kind=dp)                                      :: max_herm, buffer, total_n

      ! Get the logger
      logger => cp_get_default_logger()

      ! Only proceed if the print-key should output
      IF (cp_print_key_should_output(logger%iter_info, sbe_env%run_info_section) == 0) RETURN

      CALL timeset(routineN, handle)

      ! Get the output unit
      unit_nr = cp_print_key_unit_nr(logger, sbe_env%run_info_section, extension=".runinfo")

      ! Determine the maximum among all ranks
      ! NOTE : Idempotency deviation is from starting (not propagated) density
      CALL sbe_env%kpoints%para_env%max(max_idempot)

      ! Determine the Hermiticity deviation
      ! NOTE : Hermiticity deviation is from starting (not propagated) density
      ! Within this rank
      max_herm = 0.0_dp
      DO ispin=1, sbe_env%n_spin
         DO ikp_local=1, sbe_env%local_nkp
            buffer = nonherm_metric(sbe_env%densmat(ispin,:,:,ikp_local))
            IF (buffer > max_herm) max_herm = buffer
         END DO
      END DO

      ! Across all ranks
      CALL sbe_env%kpoints%para_env%max(max_herm)

      ! Get the total number of electrons
      total_n = 0.0_dp
      IF (sbe_env%integrand_nkp > 0) THEN
         DO ispin=1, sbe_env%n_spin

            ! Scale by spin degeneracy
            total_n = total_n + SUM(sbe_env%occ_numbers(ispin,:)) * REAL(3-sbe_env%n_spin, kind=dp)

         END DO
      END IF

      ! Proceed to printing
      IF (unit_nr > 0) WRITE(unit_nr, "(A15,I12,A6,"//dformat//",A7,"//dformat//")") &
         " SBE| Timestep ", sbe_env%sim_step, &
         " time ", sbe_env%sim_time * femtoseconds, &
         " dN ", total_n-REAL(sbe_env%n_elec_ref, kind=dp)
      IF (unit_nr > 0) WRITE(unit_nr, "(A20,"//dformat//",A20,"//dformat//")") &
         " SBE| idempot. dev. ", max_idempot, &
         " herm. dev. ", max_herm

      CALL timestop(handle)
   END SUBROUTINE output_run_info
END MODULE sbe_io
