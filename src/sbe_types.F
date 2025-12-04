!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Structures/types for Semiconductor Bloch Equations (SBE) code
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
MODULE sbe_types
   USE kpoint_types, ONLY: kpoint_type, &
                           kpoint_create, &
                           kpoint_release, &
                           read_kpoint_section
   USE kpoint_methods, ONLY: kpoint_init_cell_index
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type
   USE cp_control_types, ONLY: dft_control_type
   USE force_env_types, ONLY: force_env_type, &
                              force_env_get
   USE qs_environment_types, ONLY: qs_environment_type, &
                                   get_qs_env
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
   USE cp_dbcsr_api, ONLY: dbcsr_p_type, &
                           dbcsr_get_info
   USE cell_types, ONLY: cell_type
   USE input_section_types, ONLY: section_vals_type, &
                                  section_vals_get_subs_vals, &
                                  section_vals_val_get
   USE post_scf_bandstructure_utils, ONLY: compute_xkp
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit, &
                              cp_logger_type, &
                              cp_get_default_logger
   USE physcon, ONLY: femtoseconds
   USE input_constants, ONLY: sbe_propagator_rk4
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe'

! **************************************************************************************************
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   TYPE sbe_env_type
      TYPE(kpoint_type), POINTER                         :: kpoints => NULL()
      INTEGER                                            :: local_nkp = 1, &
                                                            mepos = 0
      ! dim 1 : spin/kind, dim 2 : rows, dim 3 : cols, dim 4 : kp
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE  :: overlap_kp, &
                                                            inv_overlap_kp, &
                                                            ham_zero_kp
      REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: overlap_rs, &
                                                            ham_zero_rs
      ! Propagation target - density matrix
      COMPLEX(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE  :: densmat, &
                                                            densmat_next, &
                                                            densmat_work
      INTEGER                                            :: n_ao = 0, &
                                                            n_images = 0, &
                                                            n_spin = 1, &
                                                            unit_nr = 6, &
                                                            sim_nsteps = -1, &
                                                            sim_start = 0, &
                                                            sim_start_orig = 0, &
                                                            sim_step = 0
      INTEGER                                            :: propagator_method = sbe_propagator_rk4
      REAL(kind=dp)                                      :: sim_dt = 1.0e-3_dp / femtoseconds, &
                                                            sim_time = 0.0_dp
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: time_trace
   END TYPE sbe_env_type

   PUBLIC :: sbe_env_type, &
             init_sbe_env, &
             release_sbe_env
CONTAINS
! **************************************************************************************************
!> \brief Initializes the sbe_env - allocates objects and retrieves parameters from DFT
!> \param sbe_env null pointer to the sbe_env type, which will be allocated and filled
!> \param force_env pointer to the force env type, which contains data from the DFT run
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE init_sbe_env(sbe_env, force_env)
      TYPE(sbe_env_type), POINTER                        :: sbe_env
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN="init_sbe_env"
      INTEGER                                            :: handle, &
                                                            i
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_all!, &
                                                            !sab_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s_kp, &
                                                            matrix_ks_kp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: kpoint_section, &
                                                            force_input
      REAL(kind=dp)                                      :: efficiency
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! Allocation done in this routine, release done in separate routine
      IF (ASSOCIATED(sbe_env)) THEN
         CALL cp_abort(__LOCATION__, "init_sbe_env expects uninitialized pointer to sbe_env.")
      END IF

      ALLOCATE (sbe_env)

      ! Start with kpoints allocation
      CALL force_env_get(force_env, qs_env=qs_env, input=force_input)
      CALL get_qs_env(qs_env, sab_all=sab_all, para_env=para_env, dft_control=dft_control, cell=cell)

      CALL kpoint_create(sbe_env%kpoints)
      CALL kpoint_init_cell_index(sbe_env%kpoints, sab_all, para_env, dft_control)
      sbe_env%kpoints%para_env => para_env
      CALL para_env%retain()
      sbe_env%n_images = SIZE(sbe_env%kpoints%index_to_cell, 2)

      ! TODO : Construct the (non-symmetric+full) kpoint grid
      ! Start by getting the info into kpoints
      ! TODO : Might be best to define a separate kpoint grid from the DFT calculation
      kpoint_section => section_vals_get_subs_vals(force_input, "DFT%KPOINTS")
      CALL read_kpoint_section(sbe_env%kpoints, kpoint_section, cell%hmat)
      SELECT CASE (sbe_env%kpoints%kp_scheme)
         CASE ("MONKHORST-PACK")
            ! Calculate the xkp
            ! TODO : Is any symmetry in kpoints available in non-linear SBE?
            sbe_env%kpoints%nkp = PRODUCT(sbe_env%kpoints%nkp_grid)
            ALLOCATE (sbe_env%kpoints%xkp(3,sbe_env%kpoints%nkp), &
                                          source=0.0_dp)
            ! The grid is accessible on all the ranks
            CALL compute_xkp(sbe_env%kpoints%xkp, 1, sbe_env%kpoints%nkp, sbe_env%kpoints%nkp_grid)
         CASE ("GENERAL")
            ! Already calculated in read_kpoint_section
         CASE DEFAULT
            ! Abort - unknown kpoint scheme
            CALL cp_abort(__LOCATION__, "Unknown kpoint scheme in SBE!")
      END SELECT
      ! TODO : Get the kpoint distribution (unless too hard, then just use modulos)
      CALL sbe_create_kpoint_dist(sbe_env%kpoints, efficiency)
      sbe_env%mepos = sbe_env%kpoints%para_env%mepos + 1
      sbe_env%local_nkp = sbe_env%kpoints%kp_dist(2, sbe_env%mepos) - &
                          sbe_env%kpoints%kp_dist(1, sbe_env%mepos) + 1
      logger => cp_get_default_logger()
      sbe_env%unit_nr = cp_logger_get_default_io_unit(logger)
      IF (sbe_env%unit_nr > 0) THEN
         WRITE (sbe_env%unit_nr, "(A36,F42.2,A2)") " SBE| Kpoint work balance efficiency", 100.0_dp * efficiency, " %"
         DO i=1, sbe_env%kpoints%para_env%num_pe
            WRITE (sbe_env%unit_nr, "(A11,I4,A13,I38,I4,A6,I4)") " SBE| Task ", i-1, &
               " kpoint range", sbe_env%kpoints%kp_dist(1,i), &
               sbe_env%kpoints%kp_dist(2,i), " total", sbe_env%local_nkp
         END DO
      END IF

      ! Continue with allocation of rs-matrices (these are defined on full, non-symmetric grid)
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_ks_kp)
      CALL dbcsr_get_info(matrix_s_kp(1,1)%matrix, nfullrows_total=sbe_env%n_ao)
      sbe_env%n_spin = SIZE(matrix_ks_kp, 1)
      ALLOCATE (sbe_env%overlap_rs(1, sbe_env%n_ao, sbe_env%n_ao, sbe_env%n_images), source=0.0_dp)
      ALLOCATE (sbe_env%ham_zero_rs(sbe_env%n_spin, sbe_env%n_ao, sbe_env%n_ao, sbe_env%n_images), source=0.0_dp)

      ! Allocate k-matrices
      ! TODO : Only allocate num_threads kpoints, to save memory?
      ALLOCATE (sbe_env%overlap_kp(1, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (sbe_env%inv_overlap_kp(1, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (sbe_env%ham_zero_kp(sbe_env%n_spin, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (sbe_env%densmat(sbe_env%n_spin, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (sbe_env%densmat_next(sbe_env%n_spin, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))
      ! TODO : Different numbers for different propagators - we start with euler/RK4, which is determined by RK4
      ALLOCATE (sbe_env%densmat_work(4, sbe_env%n_ao, &
                sbe_env%n_ao, sbe_env%local_nkp), source=CMPLX(0.0, 0.0, kind=dp))

      ! Get the other parameters for the simulation - time, number of steps etc
      CALL section_vals_val_get(force_env%root_section, "MOTION%MD%TIMESTEP", r_val=sbe_env%sim_dt)
      CALL section_vals_val_get(force_env%root_section, "MOTION%MD%STEPS", i_val=sbe_env%sim_nsteps)
      CALL section_vals_val_get(force_env%root_section, "MOTION%MD%STEP_START_VAL", i_val=sbe_env%sim_start)
      sbe_env%sim_start_orig = sbe_env%sim_start

      ! Propagator method
      CALL section_vals_val_get(force_input, "DFT%REAL_TIME_PROPAGATION%SBE%PROPAGATOR", &
                                i_val=sbe_env%propagator_method)

      ! Allocate time
      ALLOCATE (sbe_env%time_trace(sbe_env%sim_nsteps - sbe_env%sim_start + 1), source=0.0_dp)

      CALL timestop(handle)
   END SUBROUTINE init_sbe_env
! **************************************************************************************************
!> \brief Deallocates the memory and the structure
!> \param sbe_env pointer to the existing sbe_env structure
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE release_sbe_env(sbe_env)
      TYPE(sbe_env_type), POINTER                        :: sbe_env

      CPASSERT(ASSOCIATED(sbe_env))

      CALL kpoint_release(sbe_env%kpoints)
      DEALLOCATE (sbe_env%overlap_rs)
      DEALLOCATE (sbe_env%overlap_kp)
      DEALLOCATE (sbe_env%inv_overlap_kp)
      DEALLOCATE (sbe_env%ham_zero_rs)
      DEALLOCATE (sbe_env%ham_zero_kp)
      DEALLOCATE (sbe_env%time_trace)
      DEALLOCATE (sbe_env%densmat)
      DEALLOCATE (sbe_env%densmat_next)

      ! Finish by deallocating the pointer itself
      DEALLOCATE (sbe_env)
   END SUBROUTINE release_sbe_env
! **************************************************************************************************
!> \brief Calculates the kpoint distribution - tries put equal number of kpoints on each MPI rank
!> \param kpoints kpoint_type structure, containing nkp and destination for the distribution
!> \param efficiency_opt reports optional efficiency (actual work/(num_pe * max_work_per_pe))
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE sbe_create_kpoint_dist(kpoints, efficiency_opt)
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(kind=dp), OPTIONAL                            :: efficiency_opt
      INTEGER                                            :: num_pe, &
                                                            i, &
                                                            !mepos, &
                                                            per_task
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nkp_loc

      ! Get the number of tasks
      num_pe = kpoints%para_env%num_pe
      ! Allocate the distribution
      ALLOCATE (kpoints%kp_dist(2, num_pe), source=-1)
      per_task = kpoints%nkp / num_pe
      ALLOCATE (nkp_loc(num_pe), source=per_task)
      DO i=1, MODULO(kpoints%nkp, num_pe)
         nkp_loc(i) = per_task + 1
      END DO
      IF (kpoints%nkp < 1) THEN
         ! Edge case - invalid or zero number of kpoints
         ! In that case, set up so that all bounds lead to no iterations
         kpoints%kp_dist(1,:) = -1
         kpoints%kp_dist(2,:) = -2
      ELSE
         ! Standard case
         ! Start
         kpoints%kp_dist(1,1) = 1
         kpoints%kp_dist(2,1) = nkp_loc(1)
         ! Other tasks
         DO i=2, num_pe
            IF (nkp_loc(i) > 0) THEN
               kpoints%kp_dist(1,i) = kpoints%kp_dist(2,i-1) + 1
               kpoints%kp_dist(2,i) = kpoints%kp_dist(1,i) + nkp_loc(i) - 1
            ELSE
               ! No kpoints edge case - no iteration setup
               kpoints%kp_dist(1,i) = -1
               kpoints%kp_dist(2,i) = -2
            END IF
         END DO
      END IF
      ! The distribution is now setup
      ! If requested, report on efficiency of work balance distribution
      IF (PRESENT(efficiency_opt)) THEN
         efficiency_opt = REAL(kpoints%nkp, kind=dp) / ( &
            REAL(num_pe, kind=dp) * REAL(nkp_loc(1), kind=dp))
      END IF
      ! TODO : Put io task on the rank which has fewer kpoints, if such task exists
      DEALLOCATE (nkp_loc)
   END SUBROUTINE sbe_create_kpoint_dist
END MODULE sbe_types
