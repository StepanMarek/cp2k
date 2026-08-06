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
   USE input_constants,                ONLY: rtp_method_sbe, &
                                             use_mom_ref_zero, &
                                             sbe_propagator_etrs, &
                                             sbe_propagator_rk4
   USE efield_utils,                   ONLY: make_field
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
   USE sbe_io,                         ONLY: output_field, &
                                             output_occupation_file, &
                                             print_mat, &
                                             output_nonidempot_metric, &
                                             output_nonherm_metric
   USE sbe_utils,                      ONLY: print_sbe_header_info
   USE sbe_types,                      ONLY: sbe_env_type, &
                                             init_sbe_env, &
                                             release_sbe_env, &
                                             get_moments_rs_dbcsr
   USE sbe_methods,                    ONLY: init_densmat, &
                                             nonherm_metric, &
                                             nonidempot_metric, &
                                             get_damping_term, &
                                             bch_exp_propagate, &
                                             add_occ_numbers_kp, &
                                             add_current_kp
   USE cp_dbcsr_api,                   ONLY: dbcsr_p_type, &
                                             dbcsr_type, &
                                             dbcsr_create, &
                                             dbcsr_release, &
                                             dbcsr_release_p, &
                                             dbcsr_desymmetrize, &
                                             dbcsr_get_info
   USE cp_dbcsr_contrib,               ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link,             ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cell_types,                     ONLY: cell_type
   USE kpoint_k_r_trafo_simple,        ONLY: replicate_rs_matrices, &
                                             rs_to_kp
   USE cp_control_types,               ONLY: dft_control_type
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_add_iter_level, &
                                 cp_rm_iter_level, &
                                 cp_iterate
   USE physcon,                        ONLY: femtoseconds
   USE mathlib,                        ONLY: invcmat, &
                                             gemm_square
   USE mathconstants,                  ONLY: twopi
   USE qs_moments,                     ONLY: qs_moment_kpoints

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
      INTEGER                                            :: unit_nr, &
                                                            i_step, &
                                                            i, &
                                                            i_local, &
                                                            j, k, &
                                                            info, &
                                                            handle
      ! REAL(kind=dp), DIMENSION(3)                        :: single_k
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_h_kp, moments_rs_img
      ! TYPE(dbcsr_type)                                   :: s_k, sc_k
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_nl
      TYPE(sbe_env_type), POINTER                        :: sbe_env
      REAL(kind=dp)                                      :: total_n, lambda, buffer, max_buffer
      COMPLEX(kind=dp), DIMENSION(3)                     :: total_current
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp), DIMENSION(:), POINTER               :: ref_custom

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) CALL print_sbe_header_info(unit_nr)

      CALL force_env_get(force_env, qs_env=qs_env)

      ! Run the initial DFT calculation
      IF (unit_nr > 0) WRITE(unit_nr, "(A37)") " SBE| Running initial DFT calculation"
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)
      IF (unit_nr > 0) WRITE(unit_nr, "(A34)") " SBE| Initial DFT calculation done"

      CALL timeset(routineN, handle)

      ! Init the environment
      NULLIFY (sbe_env)
      CALL init_sbe_env(sbe_env, force_env)
      ! DEBUG : Test the BCH propagation
      ! sbe_env%propagator_method = sbe_propagator_etrs
      ! CALL test_bch()

      ! ! DEBUG : Print the dipole moments and the bandstructure
      ! CALL qs_moment_kpoints(qs_env, 1, use_mom_ref_zero, ref_custom, sbe_env%n_ao, sbe_env%unit_nr)

      ! Investigate h and S matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, &
                      matrix_ks_kp=matrix_h_kp, kpoints=kpoints)
      ! Get kpoint information
      CALL get_kpoint_info(kpoints, sab_nl=sab_nl)

      ! ! === Transformation by the original method
      ! CALL dbcsr_create(s_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      ! CALL dbcsr_create(sc_k, name="S(k)", template=matrix_s_kp(1,1)%matrix)
      ! CALL cp_dbcsr_alloc_block_from_nbl(s_k, sab_nl)
      ! CALL cp_dbcsr_alloc_block_from_nbl(sc_k, sab_nl)
      ! single_k(:) = xkp(:,1)
      ! CALL rskp_transform(rmatrix=s_k, cmatrix=sc_k, rsmat=matrix_s_kp, ispin=1, xkp=single_k, &
      !         cell_to_index=cell_to_index, sab_nl=sab_nl)
      ! CALL dbcsr_print(s_k)
      ! CALL dbcsr_print(sc_k)
      ! CALL dbcsr_release(s_k)
      ! CALL dbcsr_release(sc_k)
      ! === Transformation by the duplicated method
      IF (unit_nr > 0) WRITE(unit_nr, "(A29,I51)") " SBE| Number of mirror images", sbe_env%n_images
      IF (unit_nr > 0) WRITE(unit_nr, "(A52,E25.3,A3)") &
              " SBE| Approximate memory per rank for overlap matrix", REAL(sbe_env%n_spin, kind=dp) * &
              REAL(sbe_env%n_images, kind=dp) * REAL(sbe_env%n_ao, kind=dp) * REAL(sbe_env%n_ao, kind=dp) / &
              (1024.0_dp * 1024.0_dp * 1024.0_dp), "GiB"
      CALL replicate_rs_matrices(matrix_s_kp, kpoints, sbe_env%overlap_rs, sbe_env%kpoints%cell_to_index, 1)
      CALL replicate_rs_matrices(matrix_h_kp, kpoints, sbe_env%ham_zero_rs, sbe_env%kpoints%cell_to_index)

      ! Generate the moments in RS
      ! TODO : Check dependence on reference point, also in RTBSE
      NULLIFY (moments_rs_img)
      ALLOCATE (moments_rs_img(3, sbe_env%n_images))
      CALL get_moments_rs_dbcsr(qs_env, moments_rs_img)
      ! ! DEBUG
      ! CALL dbcsr_print(moments_rs_img(1,2)%matrix)
      CALL replicate_rs_matrices(moments_rs_img, sbe_env%kpoints, &
                                 sbe_env%moments_rs, sbe_env%kpoints%cell_to_index)
      ! ! DEBUG
      ! DO i=1, sbe_env%n_ao
      !    DO j=1, sbe_env%n_ao
      !       PRINT *, i, j, sbe_env%moments_rs(1, i, j, 2)
      !    END DO
      ! END DO
      DO i=1,3
         DO j=1, sbe_env%n_images
            CALL dbcsr_release_p(moments_rs_img(i, j)%matrix)
         END DO
      END DO
      DEALLOCATE (moments_rs_img)

      ! Initialization of the density matrix
      ! Needs KS orbitals - need to get kpoint overlap and ham_zero
      DO i = sbe_env%kpoints%kp_dist(1, sbe_env%mepos), sbe_env%kpoints%kp_dist(2, sbe_env%mepos)
         ! Explicitly set the fields to zero
         ! TODO : Delta kick application by constant vector only?
         sbe_env%field(:) = 0.0_dp
         sbe_env%vector_potential(:) = 0.0_dp
         CALL update_kp_matrices(sbe_env, i)
         ! DEBUG : Print out dipole matrices
         i_local = i - sbe_env%kpoints%kp_dist(1, sbe_env%mepos) + 1
         CALL rs_to_kp(sbe_env%overlap_rs(1,:,:,:), sbe_env%overlap_kp(:,:), &
                       sbe_env%kpoints%index_to_cell, sbe_env%kpoints%xkp(:,i))
         DO j = 1, sbe_env%n_spin
            CALL rs_to_kp(sbe_env%ham_zero_rs(j,:,:,:), sbe_env%ham_zero_kp(j,:,:), &
                          sbe_env%kpoints%index_to_cell, sbe_env%kpoints%xkp(:,i))
            ! ! DEBUG : Elementary LUMO-HOMO excitation
            ! Need to get eigenvectors
            sbe_env%densmat_work(1,:,:) = sbe_env%overlap_kp(:, :)
            sbe_env%densmat_work(2,:,:) = sbe_env%ham_zero_kp(j,:, :)
            CALL ZHEGV(1, 'V', 'U', sbe_env%n_ao, sbe_env%densmat_work(2,:,:), sbe_env%n_ao, &
                       sbe_env%densmat_work(1,:,:), sbe_env%n_ao, &
                       sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
            ! ! DEBUG : Print out the bandstructure
            ! IF (sbe_env%unit_nr > 0) WRITE (sbe_env%unit_nr, "(A19)") " SBE| Bandstructure"
            ! PRINT "(A23,I6,A3,E16.8E2,A3,E16.8E2,A3,E16.8E2)", " SBE| Bandstructure ikp", i, &
            !    " kx", sbe_env%kpoints%xkp(1,i), " ky", sbe_env%kpoints%xkp(2,i), " kz", sbe_env%kpoints%xkp(3,i)
            ! DO k=1,sbe_env%n_ao
            !    PRINT "(A23,I6,A4,I6,A7,E20.8E3)", &
            !       " SBE| Bandstructure ikp", i, " nmo", k, " energy", REAL(sbe_env%eigvals(k), kind=dp)
            ! END DO
            ! Eigenvectors are in densmat_work(2,:,:)
            ! Save the eigenvectors
            sbe_env%current_mo_coeffs(:,:) = sbe_env%densmat_work(2,:,:)
            ! lambda = 0.2
            ! ! Occupations in MO
            ! ! CHECK : all relevant occupations set?
            ! sbe_env%densmat_work(1,1,1) = CMPLX(COS(lambda) * COS(lambda), 0.0, kind=dp)
            ! sbe_env%densmat_work(1,2,2) = CMPLX(SIN(lambda) * SIN(lambda), 0.0, kind=dp)
            ! sbe_env%densmat_work(1,1,2) = CMPLX(0.0, -SIN(lambda) * COS(lambda), kind=dp)
            ! sbe_env%densmat_work(1,2,1) = CMPLX(0.0, SIN(lambda) * COS(lambda), kind=dp)
            ! ! Multiply from right
            ! CALL gemm_square(sbe_env%densmat_work(2,:,:), "N", sbe_env%densmat_work(1,:,:), "N", sbe_env%densmat_work(3,:,:))
            ! ! Multiply from the left
            ! CALL gemm_square(sbe_env%densmat_work(3,:,:), "N", sbe_env%densmat_work(2,:,:), "C", sbe_env%densmat(j,:,:,i_local))
            ! PRINT *, "kpoint", i, "spin", j, "density AO", sbe_env%densmat(j, :, :, i_local)
            ! Initialize the density
            MARK_USED(lambda)
            CALL init_densmat(sbe_env%densmat(j,:,:,i_local), j, &
                              sbe_env%ham_zero_kp(j,:,:), sbe_env%overlap_kp, &
                              qs_env, sbe_env)
            ! CALL output_nonidempot_metric(sbe_env%densmat(:,:,:,i_local), sbe_env%overlap_kp, &
            !    sbe_env%densmat_work, sbe_env%unit_nr)
         END DO
      END DO

      ! Delta pulse application
      ! In velocity gauge, simply change the vector potential to non-equilibrium value
      IF (sbe_env%apply_delta_pulse) THEN
         sbe_env%vector_potential(:) = sbe_env%delta_pulse_scale * &
            REAL(sbe_env%delta_pulse_direction(:), kind=dp)
      END IF

      ! RT-schemes typically use "MD" as the iteration level
      logger => cp_get_default_logger()
      CALL cp_add_iter_level(logger%iter_info, "MD")

      ! Start the actual simulation
      DO i_step=sbe_env%sim_start, sbe_env%sim_nsteps
         ! Update timestep parameters
         sbe_env%sim_step = i_step
         sbe_env%sim_time = REAL(i_step, kind=dp) * sbe_env%sim_dt
         ! Append to the trace
         sbe_env%time_trace(i_step - sbe_env%sim_start + 1) = sbe_env%sim_time
         ! Update the field to current value
         ! Vector potential holds value from last iteration
         IF (sbe_env%dft_control%apply_efield_field) THEN
            CALL make_field(sbe_env%dft_control, sbe_env%field, sbe_env%sim_step, sbe_env%sim_time)
         END IF
         ! Output the field value to the file/to the __STDOUT__
         ! Default field is zero, so should be safe
         CALL output_field(sbe_env%sim_time, sbe_env%field, &
            sbe_env%field_section, i_step == sbe_env%sim_start_orig)

         ! Reset occupation numbers
         sbe_env%occ_numbers(:,:) = 0.0_dp
         total_current(:) = CMPLX(0.0, 0.0, kind=dp)

         ! Iteration over local kpoints
         ! TODO : Run OMP over local kpoints or not?
         buffer = 0.0_dp
         max_buffer = 0.0_dp
         DO i=sbe_env%kpoints%kp_dist(1,sbe_env%mepos), sbe_env%kpoints%kp_dist(2, sbe_env%mepos)
            i_local = i - sbe_env%kpoints%kp_dist(1, sbe_env%mepos) + 1
            ! Recalculate the k-matrices at new kpoints
            ! TODO : Indicate/change how the field is propagated
            ! TODO : Maybe implement field current and field next?
            ! TODO : For fields that are no longer changing, stop the update to speed up the simulation? Maybe via metric change?
            CALL update_kp_matrices(sbe_env, i)
            ! Actual propagation of the density rho(t) -> rho(t+dt)
            DO j=1,sbe_env%n_spin
               ! Log max idempotency violation in this k-point cycle
               ! DEBUG : Idempotency violation printing
               buffer = nonidempot_metric(sbe_env%densmat(j,:,:,i_local), sbe_env%overlap_kp, sbe_env%densmat_work)
               IF (buffer > max_buffer) max_buffer = buffer
               CALL propagate_density(sbe_env, sbe_env%densmat(j,:,:,i_local), &
                                      sbe_env%densmat_next(j,:,:,i_local), j, i, info)
            END DO
            ! Export observables
            ! Per-kpoint operations - add occ_numbers
            CALL add_occ_numbers_kp(sbe_env, sbe_env%occ_numbers, sbe_env%densmat_next(:,:,:,i_local), i)
            ! add current
            ! TODO : Spin
            CALL add_current_kp(sbe_env, sbe_env%densmat(1,:,:,i_local), sbe_env%kpoints%xkp(:,i), 1, total_current(:), i)
         END DO
         ! Info can report e.g. self-consistency failures
         IF (info > 0) THEN
            CALL cp_abort(__LOCATION__, "Failed to propagate the density matrix.")
         END IF

         ! Idempotency violation
         PRINT *, "SBE| Idempotency violation", max_buffer

         ! Increment iterations
         CALL cp_iterate(logger%iter_info, iter_nr=i_step, last=(i_step == sbe_env%sim_nsteps))
         ! Outputs
         ! Sum over the occupation numbers across all ranks
         CALL sbe_env%kpoints%para_env%sum(sbe_env%occ_numbers)
         CALL sbe_env%kpoints%para_env%sum(total_current)
         ! DEBUG : total occ number, TODO : Difference to equilibrium
         ! TODO : Idempotency metric
         total_n = SUM(sbe_env%occ_numbers(1,:)) * sbe_env%n_spin
         IF (sbe_env%unit_nr > 0) WRITE(sbe_env%unit_nr, "(A15,I6,A6,F7.4,A17,E20.8E3)") &
            " SBE| Timestep ", i_step, &
            " time ", sbe_env%sim_time * femtoseconds, " N ", total_n
         IF (sbe_env%unit_nr > 0) WRITE(sbe_env%unit_nr, "(A15,I6,A6,F7.4,A17,E20.8E3,E20.8E3,E20.8E3)") &
            " SBE| Timestep ", i_step, &
            " time ", sbe_env%sim_time * femtoseconds, " current ", &
            REAL(total_current(1)), REAL(total_current(2)), REAL(total_current(3))
         ! Output occupations/vacancies into files
         CALL output_occupation_file(sbe_env%occ_numbers(1,:), sbe_env%sim_time, sbe_env%occ_section, &
            header_opt = (i_step == sbe_env%sim_start_orig), vacancy_opt=.FALSE.)
         CALL output_occupation_file(sbe_env%occ_numbers(1,:), sbe_env%sim_time, sbe_env%occ_section, &
            header_opt = (i_step == sbe_env%sim_start_orig), vacancy_opt=.TRUE.)
         ! Output the Hermiticity violation
         CALL output_nonherm_metric(sbe_env%densmat_next, sbe_env%unit_nr)
         ! Set the current density matrix for the next step
         sbe_env%densmat(:,:,:,:) = sbe_env%densmat_next(:,:,:,:)
         ! Update the vector potential - A(t+dt) = A(t) - E(t) dt
         sbe_env%vector_potential(:) = sbe_env%vector_potential(:) - sbe_env%sim_dt * sbe_env%field(:)
      END DO

      ! Stop iteratons
      CALL cp_rm_iter_level(logger%iter_info, "MD")
      ! Deallocations
      CALL release_sbe_env(sbe_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_sbe
! **************************************************************************************************
!> \brief Propagates from one density to other, switch between different propagators here
!> \note Uses the previously calculated overlap, inverse overlap and hamiltonian in kspace
!> \note Field/vector potential calculation done in here, since there is explicit t dependence
!> \param sbe_env SBE environment
!> \param densmat_old Explicit old density
!> \param densmat_new Explicit new density
!> \param info Whether the propagation was successful, useful for self-consistent propagators
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE propagate_density(sbe_env, densmat_old, densmat_new, ispin, ikp, info)
      TYPE(sbe_env_type), POINTER, INTENT(INOUT)         :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)       :: densmat_old
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)      :: densmat_new
      INTEGER, INTENT(IN)                                :: ispin, ikp
      INTEGER, INTENT(OUT)                               :: info
      CHARACTER(len=*), PARAMETER                        :: routineN="propagate_density"
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! TODO : Switch/select for different propagators
      ! RK-4
      SELECT CASE (sbe_env%propagator_method)
      CASE (sbe_propagator_rk4)
         ! k1 term
         ! Vector potential is always set from before (previous timestep or delta kick or zero)
         CALL propagator_commutator(sbe_env%densmat_work(1,:,:), &
                                    sbe_env%inv_overlap_kp(:,:), &
                                    sbe_env%ham_kp(ispin,:,:), &
                                    densmat_old(:,:), &
                                    densmat_new(:,:))
         ! Add the non-commutator terms
         ! Add the damping term
         CALL get_damping_term(sbe_env, sbe_env%damping, densmat_old(:,:), ispin)
         sbe_env%densmat_work(1,:,:) = sbe_env%densmat_work(1,:,:) + sbe_env%damping(:,:)
         ! k2 and k3 terms are evaluated at t+Delta t/2 - recalculate h, S, S^(-1) and field
         ! k2 term
         ! Update the field, mid-point
         ! TODO : Sim-step interpolation?
         IF (sbe_env%dft_control%apply_efield_field) THEN
            CALL make_field(sbe_env%dft_control, sbe_env%field_next, sbe_env%sim_step, &
                            sbe_env%sim_time + 0.5_dp*sbe_env%sim_dt)
         END IF
         ! Add to vector potential A(t + dt/2) = A(t) - E(t) dt/2
         sbe_env%vector_potential_next(:) = sbe_env%vector_potential(:) - &
                                       0.5_dp*sbe_env%sim_dt * sbe_env%field(:)
         ! Since vector potential changes, update kp matrices
         CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential_next, sbe_env%field(:))
         ! CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential_next)
         CALL propagator_commutator(sbe_env%densmat_work(2,:,:), &
                                    sbe_env%inv_overlap_kp(:,:), &
                                    sbe_env%ham_kp(ispin,:,:), &
                                    densmat_old(:,:) + &
                                    sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(1,:,:), &
                                    densmat_new(:,:))
         ! Add the non-commutator terms - damping
         CALL get_damping_term(sbe_env, sbe_env%damping, &
                               densmat_old(:,:) + sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(1,:,:), &
                               ispin)
         sbe_env%densmat_work(2,:,:) = sbe_env%densmat_work(2,:,:) + sbe_env%damping(:,:)
         ! k3 term
         CALL propagator_commutator(sbe_env%densmat_work(3,:,:), &
                                    sbe_env%inv_overlap_kp(:,:), &
                                    sbe_env%ham_kp(ispin,:,:), &
                                    densmat_old(:,:) + &
                                    sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(2,:,:), &
                                    densmat_new(:,:))
         ! Add the non-commutator terms - damping and overlap derivative
         CALL get_damping_term(sbe_env, sbe_env%damping, &
                               densmat_old(:,:) + sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(2,:,:), &
                               ispin)
         sbe_env%densmat_work(3,:,:) = sbe_env%densmat_work(3,:,:) + sbe_env%damping(:,:)
         ! k4 term is evaluated at t+Delta t - recalculate h, S, S^(-1) and field
         ! k4 term
         ! Update the field, final time
         IF (sbe_env%dft_control%apply_efield_field) THEN
            CALL make_field(sbe_env%dft_control, sbe_env%field_next, sbe_env%sim_step + 1, &
                            sbe_env%sim_time + sbe_env%sim_dt)
         END IF
         ! Add to vector potential A(t+dt) = A(t) - E(t) dt
         sbe_env%vector_potential_next(:) = sbe_env%vector_potential(:) - &
                                       sbe_env%sim_dt*sbe_env%field(:) ! Forward
         ! Since vector potential changes, update kp matrices
         CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential_next, sbe_env%field(:))
         ! CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential_next)
         CALL propagator_commutator(sbe_env%densmat_work(4,:,:), &
                                    sbe_env%inv_overlap_kp(:,:), &
                                    sbe_env%ham_kp(ispin,:,:), &
                                    densmat_old(:,:) + &
                                    sbe_env%sim_dt * sbe_env%densmat_work(3,:,:), &
                                    densmat_new(:,:))
         ! Add the non-commutator terms - damping
         CALL get_damping_term(sbe_env, sbe_env%damping, &
                               densmat_old(:,:) + sbe_env%sim_dt * sbe_env%densmat_work(3,:,:), &
                               ispin)
         sbe_env%densmat_work(4,:,:) = sbe_env%densmat_work(4,:,:) + sbe_env%damping(:,:)
         ! Final determination of next step
         densmat_new(:,:) = densmat_old(:,:) + &
            sbe_env%sim_dt * (sbe_env%densmat_work(1,:,:) + 2.0_dp * sbe_env%densmat_work(2,:,:) + &
               2.0_dp * sbe_env%densmat_work(3,:,:) + sbe_env%densmat_work(4,:,:)) / 6.0_dp
         ! Reset the kp matrices if there are other spins to propagate
         IF (sbe_env%n_spin > 1 .AND. ispin == 1) THEN
            CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential, sbe_env%field)
         END IF
         info = 0
      CASE (sbe_propagator_etrs)
         ! rho(t + dt/2) = U^dagger[H(t),dt/2] rho(t) U[H(t),dt/2]
         ! Intermediate result stored in damping
         CALL bch_exp_propagate(sbe_env%damping, densmat_old, &
            sbe_env%inv_overlap_kp, sbe_env%ham_kp(ispin,:,:), 0.5_dp * sbe_env%sim_dt, &
            1.0E-15_dp, info, sbe_env%densmat_work, 50)
         ! rho(t + dt) = U^dagger[H(t+dt),dt/2] rho(t+dt/2) U[H(t+dt),dt/2]
         ! Update field and matrices
         IF (sbe_env%dft_control%apply_efield_field) THEN
            CALL make_field(sbe_env%dft_control, sbe_env%field_next, sbe_env%sim_step + 1, &
                            sbe_env%sim_time + sbe_env%sim_dt)
         END IF
         sbe_env%vector_potential_next(:) = sbe_env%vector_potential(:) - &
            sbe_env%sim_dt * sbe_env%field(:)
         CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential_next, &
            sbe_env%field_next(:))
         CALL bch_exp_propagate(densmat_new, sbe_env%damping, &
            sbe_env%inv_overlap_kp, sbe_env%ham_kp(ispin,:,:), 0.5_dp * sbe_env%sim_dt, &
            1.0E-15_dp, info, sbe_env%densmat_work, 50)
         IF (sbe_env%n_spin > 1 .AND. ispin == 1) THEN
            CALL update_kp_matrices(sbe_env, ikp, sbe_env%vector_potential, sbe_env%field)
         END IF
      CASE DEFAULT
         ! Matrices are updated from before
         ! Default to forwards Euler
         CALL propagator_commutator(sbe_env%densmat_work(1,:,:), &
                                    sbe_env%inv_overlap_kp(:,:), &
                                    sbe_env%ham_kp(ispin,:,:), &
                                    densmat_old(:,:), &
                                    densmat_new(:,:))
         ! Add the non-commutator terms - damping and overlap derivative
         CALL get_damping_term(sbe_env, sbe_env%damping, densmat_old(:,:), ispin)
         sbe_env%densmat_work(1,:,:) = sbe_env%densmat_work(1,:,:) + sbe_env%damping(:,:)
         ! Propagate
         densmat_new(:,:) = densmat_old(:,:) + &
            sbe_env%sim_dt * sbe_env%densmat_work(1,:,:)
         info = 0
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE propagate_density
! **************************************************************************************************
!> \brief Calculates the propagator commutator, used in different propagator methods
!> \param commutator Resulting commutator matrix
!> \param inv_overlap Inverse overlap matrix
!> \param ham Hamiltonian matrix
!> \param densmat density matrix
!> \param work workspace, same dimensions as the matrices
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE propagator_commutator(commutator, inv_overlap, ham, densmat, work)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)       :: densmat, &
                                                            inv_overlap, &
                                                            ham
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)      :: commutator, &
                                                            work

      ! S^(-1) h rho(t)
      CALL gemm_square(inv_overlap, 'N', ham, 'C', densmat, 'N', work)
      ! rho (t) h S^(-1)
      CALL gemm_square(densmat, 'N', ham, 'N', inv_overlap, 'N', commutator)

      commutator(:,:) = CMPLX(0.0, -1.0, kind=dp) * (work(:,:) - commutator(:,:))
   END SUBROUTINE propagator_commutator
! **************************************************************************************************
!> \brief Updates kpoint matrices based on the current value of field and vector potential
! **************************************************************************************************
   SUBROUTINE update_kp_matrices(sbe_env, ikp, current_vector_opt, field_opt)
      TYPE(sbe_env_type), POINTER, INTENT(INOUT)         :: sbe_env
      INTEGER, INTENT(IN)                                :: ikp
      REAL(kind=dp), DIMENSION(3), OPTIONAL              :: current_vector_opt, &
                                                            field_opt
      REAL(kind=dp), DIMENSION(3)                        :: xkp_current, &
                                                            field
      INTEGER                                            :: ikp_local, &
                                                            info, &
                                                            ispin, &
                                                            j

      ikp_local = ikp - sbe_env%kpoints%kp_dist(1, sbe_env%mepos) + 1
      ! Calculate the xkp based on equilibrium xkp and vector potential
      IF (PRESENT(current_vector_opt)) THEN
         xkp_current(:) = sbe_env%kpoints%xkp(:,ikp) - &
            MATMUL(TRANSPOSE(sbe_env%hmat), current_vector_opt(:)) / twopi
      ELSE
         xkp_current(:) = sbe_env%kpoints%xkp(:,ikp) - &
            MATMUL(TRANSPOSE(sbe_env%hmat), sbe_env%vector_potential(:)) / twopi
      END IF
      ! Determine the field
      IF (PRESENT(field_opt)) THEN
         field(:) = field_opt(:)
      ELSE
         field(:) = sbe_env%field(:)
      END IF
      ! Now, with the kpoint coordinate, calculate the matrices
      ! Start with overlap
      CALL rs_to_kp(sbe_env%overlap_rs(1,:,:,:), sbe_env%overlap_kp(:,:), &
                    sbe_env%kpoints%index_to_cell, xkp_current)
      ! Calculate the inverse overlap
      sbe_env%inv_overlap_kp(:,:) = sbe_env%overlap_kp(:,:)
      CALL invcmat(sbe_env%inv_overlap_kp, info)
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, "Overlap matrix inverse failed.")
      END IF
      ! Continue with spin dependent operators
      DO ispin=1, sbe_env%n_spin
         ! Equilibrium Hamiltonian
         CALL rs_to_kp(sbe_env%ham_zero_rs(ispin,:,:,:), sbe_env%ham_zero_kp(ispin,:,:), &
                    sbe_env%kpoints%index_to_cell, xkp_current)
      END DO
      ! Continue with direction dependent operators
      DO j = 1, 3
         ! Moments
         CALL rs_to_kp(sbe_env%moments_rs(j,:,:,:), sbe_env%moments_kp(j,:,:), &
                       sbe_env%kpoints%index_to_cell, xkp_current)
         ! Overlap gradient
         CALL rs_to_kp(sbe_env%overlap_rs(1,:,:,:), sbe_env%overlap_grad(j,:,:), &
                       sbe_env%kpoints%index_to_cell, xkp_current, &
                       j, sbe_env%hmat)
      END DO
      ! With moments in place, construct the effective Hamiltonian
      DO ispin = 1, sbe_env%n_spin
         sbe_env%ham_kp(ispin, :, :) = sbe_env%ham_zero_kp(ispin,:,:)
         DO j = 1, 3
            sbe_env%ham_kp(ispin, :, :) = sbe_env%ham_kp(ispin, :, :) - &
               field(j) * sbe_env%moments_kp(j,:,:)
         END DO
      END DO
   END SUBROUTINE update_kp_matrices
   SUBROUTINE test_bch()
      COMPLEX(kind=dp), DIMENSION(2,2)          :: densmat_new, densmat_old, hmat, smat
      COMPLEX(kind=dp), DIMENSION(4,2,2)        :: workspace
      INTEGER                                   :: converged
      REAL(kind=dp)                             :: eps, timestep

      densmat_old(:,:) = CMPLX(0.0, 0.0, kind=dp)
      densmat_old(1,1) = CMPLX(1.0, 0.0, kind=dp)
      densmat_old(1,2) = CMPLX(0.0, -0.5, kind=dp)
      densmat_old(2,1) = CMPLX(0.0, 0.5, kind=dp)

      timestep = 0.1_dp
      eps = 1.0E-10_dp

      smat(:,:) = CMPLX(0.0, 0.0, kind=dp)
      smat(1,1) = CMPLX(1.0, 0.0, kind=dp)
      smat(2,2) = CMPLX(1.0, 0.0, kind=dp)

      hmat(:,:) = CMPLX(0.0, 0.0, kind=dp)
      hmat(1,2) = CMPLX(1.0, 0.0, kind=dp)
      hmat(2,1) = CMPLX(1.0, 0.0, kind=dp)

      CALL bch_exp_propagate(densmat_new, densmat_old, smat, hmat, timestep, eps, converged, &
         workspace, 30)

      CALL print_mat(densmat_new)

      CPABORT("EXP test")
   END SUBROUTINE test_bch
END MODULE sbe
