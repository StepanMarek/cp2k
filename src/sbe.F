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
   USE input_constants,                ONLY: rtp_method_sbe
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
   USE sbe_utils,                      ONLY: print_sbe_header_info
   USE sbe_types,                      ONLY: sbe_env_type, &
                                             init_sbe_env, &
                                             release_sbe_env
   USE cp_dbcsr_api,                   ONLY: dbcsr_p_type, &
                                             dbcsr_type, &
                                             dbcsr_create, &
                                             dbcsr_release, &
                                             dbcsr_desymmetrize, &
                                             dbcsr_get_info
   USE cp_dbcsr_contrib,               ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link,             ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cell_types,                     ONLY: cell_type
   USE kpoint_k_r_trafo_simple,        ONLY: replicate_rs_matrices, &
                                             rs_to_kp
   USE cp_control_types,               ONLY: dft_control_type
   USE physcon,                        ONLY: femtoseconds
   USE mathlib,                        ONLY: invcmat, &
                                             gemm_square

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
                                                            j, &
                                                            info, &
                                                            handle
      ! REAL(kind=dp), DIMENSION(3)                        :: single_k
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_h_kp
      ! TYPE(dbcsr_type)                                   :: s_k, sc_k
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_nl
      TYPE(sbe_env_type), POINTER                        :: sbe_env

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
      IF (unit_nr > 0) WRITE(unit_nr, "(A52,E25.3,A3)") &
              " SBE| Approximate memory per rank for overlap matrix", REAL(sbe_env%n_spin, kind=dp) * &
              REAL(sbe_env%n_images, kind=dp) * REAL(sbe_env%n_ao, kind=dp) * REAL(sbe_env%n_ao, kind=dp) / &
              (1024.0_dp * 1024.0_dp * 1024.0_dp), "GiB"
      CALL replicate_rs_matrices(matrix_s_kp, kpoints, sbe_env%overlap_rs, sbe_env%kpoints%cell_to_index)
      CALL replicate_rs_matrices(matrix_h_kp, kpoints, sbe_env%ham_zero_rs, sbe_env%kpoints%cell_to_index)

      ! DEBUG
      sbe_env%densmat(1,1,1,1) = CMPLX(1.0,0.0,kind=dp)
      sbe_env%densmat(1,2,1,1) = CMPLX(0.0,1.0,kind=dp)
      sbe_env%densmat(1,1,2,1) = CMPLX(0.0,-1.0,kind=dp)

      ! Start the actual simulation
      DO i_step=sbe_env%sim_start, sbe_env%sim_nsteps
         ! Update timestep parameters
         sbe_env%sim_step = i_step
         sbe_env%sim_time = REAL(i_step, kind=dp) * sbe_env%sim_dt
         ! Append to the trace
         sbe_env%time_trace(i_step - sbe_env%sim_start + 1) = sbe_env%sim_time
         ! TODO : Update the kpoints to set the field

         ! Iteration over local kpoints
         ! TODO : Run OMP over local kpoints or not?
         DO i=sbe_env%kpoints%kp_dist(1,sbe_env%mepos), sbe_env%kpoints%kp_dist(2, sbe_env%mepos)
            i_local = i - sbe_env%kpoints%kp_dist(1, sbe_env%mepos) + 1
            ! TODO : Recalculate the k-matrices at new kpoints
            ! TODO : For fields that are no longer changing, stop the update to speed up the simulation? Maybe via metric change?
            ! Overlap
            CALL rs_to_kp(sbe_env%overlap_rs(1,:,:,:), sbe_env%overlap_kp(1,:,:,i_local), &
                          sbe_env%kpoints%index_to_cell, sbe_env%kpoints%xkp(:,i))
            ! Inverse overlap
            sbe_env%inv_overlap_kp(1,:,:,i_local) = sbe_env%overlap_kp(1,:,:,i_local)
            CALL invcmat(sbe_env%inv_overlap_kp(1,:,:,i_local), info)
            IF (info /= 0) THEN
               CALL cp_abort(__LOCATION__, "Overlap matrix inverse failed.")
            END IF
            ! Hamiltonian
            DO j=1,sbe_env%n_spin
               CALL rs_to_kp(sbe_env%ham_zero_rs(j,:,:,:), sbe_env%ham_zero_kp(j,:,:,i_local), &
                             sbe_env%kpoints%index_to_cell, sbe_env%kpoints%xkp(:,i))
            END DO

            ! TODO : Propagate the density
            DO j=1, sbe_env%n_spin
               CALL propagate_density(sbe_env, sbe_env%densmat, sbe_env%densmat_next, info)
            END DO
         END DO
         ! TODO : commit to timestep - check density criteria fulfilled in SC schemes
         IF (info > 0) THEN
            CALL cp_abort(__LOCATION__, "Failed to propagate the density matrix.")
         END IF
         sbe_env%densmat(:,:,:,:) = sbe_env%densmat_next(:,:,:,:)
         
         ! Outputs
         IF (sbe_env%unit_nr > 0) WRITE(sbe_env%unit_nr, "(A15,I6,A6,F7.4,A17,E12.3E3)") &
            " SBE| Timestep ", i_step, &
            " time ", sbe_env%sim_time * femtoseconds, " re rho(1,1,1,1) ", REAL(sbe_env%densmat_next(1,1,1,1))
      END DO

      ! Deallocations
      CALL release_sbe_env(sbe_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_sbe
! **************************************************************************************************
!> \brief Propagates from one density to other, switch between different propagators here
!> \param sbe_env SBE environment
!> \param densmat_old Explicit old density
!> \param densmat_new Explicit new density
!> \param info Whether the propagation was successful, useful for self-consistent propagators
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE propagate_density(sbe_env, densmat_old, densmat_new, info)
      TYPE(sbe_env_type), POINTER, INTENT(IN)            :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), INTENT(IN)   :: densmat_old
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), INTENT(OUT)  :: densmat_new
      INTEGER, INTENT(OUT)                               :: info
      CHARACTER(len=*), PARAMETER                        :: routineN="propagate_density"
      INTEGER                                            :: handle, &
                                                            ispin, &
                                                            ikp, &
                                                            propagator_method

      CALL timeset(routineN, handle)

      ! DEBUG
      propagator_method = 1
      ! TODO : Switch/select for different propagators
      DO ispin=1,sbe_env%n_spin
         DO ikp=1,sbe_env%local_nkp
            ! RK-4
            IF (propagator_method == 4) THEN
               CALL propagator_commutator(sbe_env%densmat_work(1,:,:,ikp), &
                                          sbe_env%inv_overlap_kp(1,:,:,ikp), &
                                          sbe_env%ham_zero_kp(ispin,:,:,ikp), &
                                          densmat_old(ispin,:,:,ikp), &
                                          densmat_new(ispin,:,:,ikp))
               CALL propagator_commutator(sbe_env%densmat_work(2,:,:,ikp), &
                                          sbe_env%inv_overlap_kp(1,:,:,ikp), &
                                          sbe_env%ham_zero_kp(ispin,:,:,ikp), &
                                          densmat_old(ispin,:,:,ikp) + &
                                          sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(1,:,:,ikp), &
                                          densmat_new(ispin,:,:,ikp))
               CALL propagator_commutator(sbe_env%densmat_work(3,:,:,ikp), &
                                          sbe_env%inv_overlap_kp(1,:,:,ikp), &
                                          sbe_env%ham_zero_kp(ispin,:,:,ikp), &
                                          densmat_old(ispin,:,:,ikp) + &
                                          sbe_env%sim_dt * 0.5_dp * sbe_env%densmat_work(2,:,:,ikp), &
                                          densmat_new(ispin,:,:,ikp))
               CALL propagator_commutator(sbe_env%densmat_work(4,:,:,ikp), &
                                          sbe_env%inv_overlap_kp(1,:,:,ikp), &
                                          sbe_env%ham_zero_kp(ispin,:,:,ikp), &
                                          densmat_old(ispin,:,:,ikp) + &
                                          sbe_env%sim_dt * sbe_env%densmat_work(3,:,:,ikp), &
                                          densmat_new(ispin,:,:,ikp))
               densmat_new(ispin,:,:,ikp) = densmat_old(ispin,:,:,ikp) + &
                  sbe_env%sim_dt * (sbe_env%densmat_work(1,:,:,ikp) + 2.0_dp * sbe_env%densmat_work(2,:,:,ikp) + &
                     2.0_dp * sbe_env%densmat_work(3,:,:,ikp) + sbe_env%densmat_work(4,:,:,ikp)) / 6.0_dp
            ELSE
               ! Default to forwards Euler
               CALL propagator_commutator(sbe_env%densmat_work(1,:,:,ikp), &
                                          sbe_env%inv_overlap_kp(1,:,:,ikp), &
                                          sbe_env%ham_zero_kp(ispin,:,:,ikp), &
                                          densmat_old(ispin,:,:,ikp), &
                                          densmat_new(ispin,:,:,ikp))
               densmat_new(ispin,:,:,ikp) = densmat_old(ispin,:,:,ikp) + &
                  sbe_env%sim_dt * sbe_env%densmat_work(1,:,:,ikp)
            END IF
         END DO
      END DO
      info = 0

      CALL timestop(handle)
   END SUBROUTINE propagate_density
   SUBROUTINE propagator_commutator(commutator, inv_overlap, ham, densmat, work)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)       :: densmat, &
                                                            inv_overlap, &
                                                            ham
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)      :: commutator, &
                                                            work

      ! S^(-1) h rho(t)
      CALL gemm_square(inv_overlap, 'N', ham, 'N', densmat, 'N', work)
      ! rho (t) h S^(-1)
      CALL gemm_square(densmat, 'N', ham, 'N', inv_overlap, 'N', commutator)

      commutator(:,:) = CMPLX(0.0, -1.0, kind=dp) * (work(:,:) - commutator(:,:))
   END SUBROUTINE propagator_commutator
END MODULE sbe
