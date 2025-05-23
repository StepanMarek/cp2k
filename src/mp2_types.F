!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! **************************************************************************************************
MODULE mp2_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_finalize,&
                                              cp_eri_mme_param
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type,&
                                              pair_list_element_type
   USE input_constants,                 ONLY: &
        do_eri_mme, eri_default, gw_pade_approx, kp_weights_W_auto, mp2_method_direct, &
        mp2_method_gpw, mp2_method_none, mp2_ri_optimize_basis, ri_mp2_laplace, ri_mp2_method_gpw, &
        ri_rpa_g0w0_crossing_z_shot, ri_rpa_method_gpw, rpa_exchange_none, soc_none, &
        wfc_mm_style_gemm
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE local_gemm_api,                  ONLY: local_gemm_ctxt_type
   USE message_passing,                 ONLY: mp_request_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_types'

   PUBLIC :: mp2_type, &
             integ_mat_buffer_type, &
             integ_mat_buffer_type_2D, &
             mp2_method_none, &
             mp2_method_direct, &
             mp2_method_gpw, &
             mp2_ri_optimize_basis, &
             ri_mp2_method_gpw, &
             ri_rpa_method_gpw, &
             ri_mp2_laplace, &
             init_TShPSC_lmax

   PUBLIC :: mp2_env_create, &
             mp2_env_release, &
             mp2_biel_type, &
             pair_list_type_mp2, &
             one_dim_int_array, &
             two_dim_int_array, &
             one_dim_real_array, &
             two_dim_real_array, &
             three_dim_real_array

   INTEGER, SAVE :: init_TShPSC_lmax = -1

! TYPE definitions

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: array
   END TYPE one_dim_int_array

   TYPE two_dim_int_array
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: array
   END TYPE two_dim_int_array

   TYPE one_dim_real_array
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: array
   END TYPE one_dim_real_array

   TYPE two_dim_real_array
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: array
   END TYPE two_dim_real_array

   TYPE three_dim_real_array
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE     :: array
   END TYPE three_dim_real_array

   TYPE mp2_biel_type
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: index_table
   END TYPE mp2_biel_type

   TYPE mp2_laplace_type
      INTEGER                                            :: n_quadrature = -1, &
                                                            num_integ_groups = -1
   END TYPE mp2_laplace_type

   TYPE mp2_direct_type
      LOGICAL                                            :: big_send = .FALSE.
   END TYPE

   TYPE mp2_gpw_type
      REAL(KIND=dp)                                      :: eps_grid = 0.0_dp, &
                                                            eps_filter = 0.0_dp, &
                                                            eps_pgf_orb_S = 0.0_dp
      INTEGER                                            :: print_level = 0
      REAL(KIND=dp)                                      :: cutoff = 0.0_dp, &
                                                            relative_cutoff = 0.0_dp
      INTEGER                                            :: size_lattice_sum = 0
   END TYPE

   TYPE ri_mp2_type
      INTEGER                                            :: block_size = 0, &
                                                            number_integration_groups = 0
      LOGICAL                                            :: print_dgemm_info = .FALSE.
   END TYPE ri_mp2_type

   TYPE ri_rpa_type
      INTEGER                                            :: rpa_num_quad_points = -1, &
                                                            rpa_num_integ_groups = -1, &
                                                            mm_style = wfc_mm_style_gemm
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data => NULL()
      TYPE(section_vals_type), POINTER                   :: xc_section_primary => Null(), &
                                                            xc_section_aux => Null()
      LOGICAL                                            :: reuse_hfx = .FALSE., &
                                                            minimax_quad = .FALSE., &
                                                            do_ri_g0w0 = .FALSE., &
                                                            do_admm = .FALSE., &
                                                            do_rse = .FALSE., &
                                                            print_dgemm_info = .FALSE.

      ! GCC 8 has an issue with this being an ALLOCATABLE
      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: mo_coeff_o => NULL(), &
                                                            mo_coeff_v => NULL()
      INTEGER                                            :: exchange_correction = rpa_exchange_none, &
                                                            exchange_block_size = -1
      LOGICAL                                            :: use_hfx_implementation = .FALSE.
      REAL(KIND=dp)                                      :: ener_exchange = 0.0_dp, &
                                                            rse_corr_diag = 0.0_dp, &
                                                            rse_corr = 0.0_dp, &
                                                            scale_rpa = 0.0_dp

      !LOGICAL                                            :: do_sigma = .FALSE.
      INTEGER                                            :: sigma_param = 0.0_dp
      REAL(KIND=dp)                                      :: e_sigma_corr = 0.0_dp

   END TYPE ri_rpa_type

   TYPE ri_rpa_im_time_type
      INTEGER                                            :: cut_memory = 0
      LOGICAL                                            :: memory_info = .FALSE., &
                                                            make_chi_pos_definite = .FALSE., &
                                                            make_overlap_mat_ao_pos_definite = .FALSE., &
                                                            trunc_coulomb_ri_x = .FALSE., &
                                                            keep_quad = .FALSE., &
                                                            do_kpoints_from_Gamma = .FALSE., &
                                                            do_extrapolate_kpoints = .FALSE.
      REAL(KIND=dp)                                      :: eps_filter = 0.0_dp, &
                                                            eps_filter_factor = 0.0_dp, &
                                                            eps_compress = 0.0_dp, &
                                                            exp_tailored_weights = 0.0_dp, &
                                                            regularization_RI = 0.0_dp, &
                                                            eps_eigval_S = 0.0_dp, &
                                                            eps_eigval_S_Gamma = 0.0_dp, &
                                                            rel_cutoff_trunc_coulomb_ri_x = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tau_tj => NULL(), &
                                                            tau_wj => NULL(), &
                                                            tj => NULL(), &
                                                            wj => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: weights_cos_tf_t_to_w => NULL(), &
                                                            weights_cos_tf_w_to_t => NULL()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: Eigenval_Gamma, &
                                                            wkp_V
      INTEGER                                            :: group_size_P = 0, &
                                                            group_size_3c = 0, &
                                                            kpoint_weights_W_method = kp_weights_W_auto, &
                                                            k_mesh_g_factor = 0
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid => NULL()
      INTEGER, DIMENSION(3)                              :: kp_grid_extra = -1
      LOGICAL                                            :: do_im_time_kpoints = .FALSE.
      INTEGER                                            :: min_bsize = 0, &
                                                            min_bsize_mo = 0, &
                                                            nkp_orig = 0, &
                                                            nkp_extra = 0
      TYPE(kpoint_type), POINTER                         :: kpoints_G => NULL(), &
                                                            kpoints_Sigma => NULL(), &
                                                            kpoints_Sigma_no_xc => NULL()
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: starts_array_mc_RI, ends_array_mc_RI, &
                                                            starts_array_mc_block_RI, &
                                                            ends_array_mc_block_RI, &
                                                            starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block

   END TYPE ri_rpa_im_time_type

   TYPE ri_g0w0_type
      INTEGER                                            :: corr_mos_occ = 0, &
                                                            corr_mos_virt = 0, &
                                                            corr_mos_occ_beta = 0, &
                                                            corr_mos_virt_beta = 0, &
                                                            num_poles = 0, &
                                                            nparam_pade = 0, &
                                                            analytic_continuation = gw_pade_approx
      REAL(KIND=dp)                                      :: omega_max_fit = 0.0_dp
      INTEGER                                            :: crossing_search = ri_rpa_g0w0_crossing_z_shot
      REAL(KIND=dp)                                      :: fermi_level_offset = 0.0_dp
      INTEGER                                            :: iter_evGW = 0, &
                                                            iter_sc_GW0 = 0
      REAL(KIND=dp)                                      :: eps_iter = 0.0_dp
      LOGICAL                                            :: do_hedin_shift = .FALSE., &
                                                            do_ri_Sigma_x = .FALSE., &
                                                            do_periodic = .FALSE., &
                                                            print_self_energy = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_Sigma_x_minus_vxc_gw
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid => NULL(), &
                                                            kp_grid_Sigma => NULL()
      INTEGER                                            :: num_kp_grids = 0
      REAL(KIND=dp)                                      :: eps_kpoint = 0.0_dp
      LOGICAL                                            :: do_mo_coeff_gamma = .FALSE., &
                                                            do_average_deg_levels = .FALSE.
      REAL(KIND=dp)                                      :: eps_eigenval = 0.0_dp
      LOGICAL                                            :: do_extra_kpoints = .FALSE., &
                                                            do_aux_bas_gw = .FALSE.
      REAL(KIND=dp)                                      :: frac_aux_mos = 0.0_dp
      INTEGER                                            :: num_omega_points = 0
      LOGICAL                                            :: do_ic_model = .FALSE., &
                                                            print_ic_values = .FALSE.
      REAL(KIND=dp)                                      :: eps_dist = 0.0_dp
      TYPE(one_dim_real_array), DIMENSION(2)             :: ic_corr_list = one_dim_real_array(NULL())
      INTEGER                                            :: print_exx = 0
      LOGICAL                                            :: do_gamma_only_sigma = .FALSE.
      LOGICAL                                            :: update_xc_energy = .FALSE., &
                                                            do_kpoints_Sigma = .FALSE., &
                                                            print_local_bandgap = .FALSE.
      INTEGER                                            :: n_kp_in_kp_line = 0, &
                                                            n_special_kp = 0, &
                                                            nkp_self_energy = 0, &
                                                            nkp_self_energy_special_kp = 0, &
                                                            nkp_self_energy_monkh_pack = 0, &
                                                            soc_type = soc_none
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp_special_kp
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE      :: matrix_sigma_x_minus_vxc, &
                                                            matrix_ks
      REAL(KIND=dp)                                      :: broadening_print_loc_bandgap = 0.0_dp, &
                                                            energy_window_print_loc_bandgap = 0.0_dp, &
                                                            ldos_thresh_print_loc_bandgap = 0.0_dp, &
                                                            energy_spacing_print_loc_bandgap = 0.0_dp, &
                                                            regularization_minimax = 0.0_dp, &
                                                            soc_energy_window = 0.0_dp
      INTEGER, DIMENSION(:), POINTER                     :: stride_loc_bandgap => NULL()

      !GW_DOS
      REAL(KIND=dp)                                      :: dos_upper = 0.0_dp, &
                                                            dos_lower = 0.0_dp, &
                                                            dos_prec = 0.0_dp, &
                                                            dos_eta = 0.0_dp
      INTEGER                                            :: max_level_self_energy = 0, &
                                                            min_level_self_energy = 0, &
                                                            dos_min = 0, &
                                                            dos_max = 0
   END TYPE ri_g0w0_type

   TYPE ri_basis_opt
      REAL(KIND=dp)                                      :: DI_rel = 0.0_dp, &
                                                            DRI = 0.0_dp, &
                                                            eps_step = 0.0_dp
      INTEGER                                            :: max_num_iter = 0, &
                                                            basis_quality = 0
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: RI_nset_per_l
   END TYPE ri_basis_opt

   TYPE grad_util
      TYPE(two_dim_real_array), DIMENSION(2)             :: P_ij = two_dim_real_array(NULL()), &
                                                            P_ab = two_dim_real_array(NULL())
      TYPE(three_dim_real_array), DIMENSION(2)           :: Gamma_P_ia = three_dim_real_array(NULL())
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: operator_half, &
                                                            PQ_half, &
                                                            Gamma_PQ, &
                                                            Gamma_PQ_2
      TYPE(dbcsr_p_type), DIMENSION(:, :), ALLOCATABLE   :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE      :: mo_coeff_o, &
                                                            mo_coeff_v
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: P_mo, W_mo, L_jb
      REAL(KIND=dp)                                      :: cphf_eps_conv = 0.0_dp, &
                                                            scale_step_size = 0.0_dp
      INTEGER                                            :: cphf_max_num_iter = 0, &
                                                            z_solver_method = 0, &
                                                            cphf_restart = 0
      LOGICAL                                            :: enforce_decrease = .FALSE., &
                                                            recalc_residual = .FALSE., &
                                                            polak_ribiere = .FALSE.
      TYPE(qs_p_env_type), POINTER                       :: p_env => NULL()
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: mp2_force => NULL()
      REAL(KIND=dp), DIMENSION(3, 3)                     :: mp2_virial = 0.0_dp
      REAL(dp)                                           :: eps_canonical = 0.0_dp
      LOGICAL                                            :: free_hfx_buffer = .FALSE.
      INTEGER                                            :: dot_blksize = 0
      INTEGER                                            :: max_parallel_comm = 0
   END TYPE grad_util

   TYPE bse_type
      INTEGER                                            :: bse_spin_config = 0, &
                                                            bse_diag_method = 0, &
                                                            flag_tda = 0, &
                                                            num_exc_en = 0, &
                                                            num_print_exc = 0, &
                                                            num_print_exc_descr = 0, &
                                                            screening_method = 0, &
                                                            num_add_start_z_space = 0, &
                                                            fac_max_z_space = 0, &
                                                            num_new_t = 0, &
                                                            num_davidson_iter = 0, &
                                                            davidson_abort_cond = 0
      REAL(KIND=dp)                                      :: eps_res = 0.0_dp, &
                                                            eps_exc_en = 0.0_dp, &
                                                            eps_x = 0.0_dp, &
                                                            screening_factor = 0.0_dp, &
                                                            bse_cutoff_occ = 0.0_dp, &
                                                            bse_cutoff_empty = 0.0_dp, &
                                                            z_space_energy_cutoff = 0.0_dp
      LOGICAL                                            :: do_bse = .FALSE., &
                                                            bse_debug_print = .FALSE., &
                                                            print_directional_exc_descr = .FALSE., &
                                                            use_ks_energies = .FALSE.
      !BSE optical spectrum
      REAL(KIND=dp)                                      :: bse_spectrum_freq_step_size = 0.0_dp, &
                                                            bse_spectrum_freq_start = 0.0_dp, &
                                                            bse_spectrum_freq_end = 0.0_dp
      LOGICAL                                            :: bse_print_spectrum = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER               :: bse_eta_spectrum_list => NULL()

      !BSE NTO prints
      LOGICAL                                            :: do_nto_analysis = .FALSE., &
                                                            explicit_nto_list = .FALSE.
      REAL(KIND=dp)                                      :: eps_nto_eigval = 0.0_dp, &
                                                            eps_nto_osc_str = 0.0_dp
      INTEGER                                            :: num_print_exc_ntos = 0
      INTEGER, DIMENSION(:), POINTER                     :: bse_nto_state_list => NULL()
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: bse_nto_state_list_final
   END TYPE bse_type

   TYPE mp2_type
      INTEGER                                            :: method = mp2_method_none
      TYPE(mp2_laplace_type)                             :: ri_laplace = mp2_laplace_type()
      TYPE(mp2_direct_type)                              :: direct_canonical = mp2_direct_type()
      TYPE(libint_potential_type)                        :: potential_parameter = libint_potential_type()
      TYPE(mp2_gpw_type)                                 :: mp2_gpw = mp2_gpw_type()
      TYPE(ri_mp2_type)                                  :: ri_mp2 = ri_mp2_type()
      TYPE(ri_rpa_type)                                  :: ri_rpa = ri_rpa_type()
      ! There is a bug with some older compilers preventing requiring an explicit initialization of allocatable components
#if defined(FTN_NO_DEFAULT_INIT)
      TYPE(ri_rpa_im_time_type)                          :: ri_rpa_im_time = ri_rpa_im_time_type(Eigenval_Gamma=NULL(), &
                                                                                                 wkp_V=NULL(), &
                                                                                                 starts_array_mc_RI=NULL(), &
                                                                                                 ends_array_mc_RI=NULL(), &
                                                                                                 starts_array_mc_block_RI=NULL(), &
                                                                                                 ends_array_mc_block_RI=NULL(), &
                                                                                     starts_array_mc=NULL(), ends_array_mc=NULL(), &
                                                                                                 starts_array_mc_block=NULL(), &
                                                                                                 ends_array_mc_block=NULL())
      TYPE(ri_g0w0_type)                                 :: ri_g0w0 = ri_g0w0_type(vec_Sigma_x_minus_vxc_gw=NULL(), &
                                                                                   xkp_special_kp=NULL(), &
                                                                                   matrix_sigma_x_minus_vxc=NULL(), &
                                                                                   matrix_ks=NULL())
      TYPE(ri_basis_opt)                                 :: ri_opt_param = ri_basis_opt(RI_nset_per_l=NULL())
      TYPE(grad_util)                                    :: ri_grad = grad_util(operator_half=NULL(), &
                                                                                PQ_half=NULL(), &
                                                                                Gamma_PQ=NULL(), &
                                                                                Gamma_PQ_2=NULL(), &
                                                                                G_P_ia=NULL(), &
                                                                                mo_coeff_o=NULL(), &
                                                                                mo_coeff_v=NULL(), &
                                                                                P_mo=NULL(), W_mo=NULL(), L_jb=NULL())
      TYPE(bse_type)                                     :: bse = bse_type(bse_nto_state_list_final=NULL())
#else
      TYPE(ri_rpa_im_time_type)                          :: ri_rpa_im_time = ri_rpa_im_time_type()
      TYPE(ri_g0w0_type)                                 :: ri_g0w0 = ri_g0w0_type()
      TYPE(ri_basis_opt)                                 :: ri_opt_param = ri_basis_opt()
      TYPE(grad_util)                                    :: ri_grad = grad_util()
      TYPE(bse_type)                                     :: bse = bse_type()
#endif
      REAL(KIND=dp)                                      :: mp2_memory = 0.0_dp, &
                                                            scale_S = 0.0_dp, &
                                                            scale_T = 0.0_dp
      INTEGER                                            :: mp2_num_proc = 0
      INTEGER                                            :: block_size_row = 0
      INTEGER                                            :: block_size_col = 0
      LOGICAL                                            :: calc_PQ_cond_num = .FALSE.
      LOGICAL                                            :: hf_fail = .FALSE.
      LOGICAL                                            :: p_screen = .FALSE.
      LOGICAL                                            :: not_last_hfx = .FALSE.
      LOGICAL                                            :: do_im_time = .FALSE.
      INTEGER                                            :: eri_method = eri_default
      TYPE(cp_eri_mme_param), POINTER                    :: eri_mme_param => NULL()
      INTEGER, DIMENSION(:), POINTER                     :: eri_blksize => NULL()
      LOGICAL                                            :: do_svd = .FALSE.
      REAL(KIND=dp)                                      :: eps_svd = -1.0_dp
      REAL(KIND=dp)                                      :: eps_range = 0.0_dp
      TYPE(libint_potential_type)                        :: ri_metric = libint_potential_type()
      TYPE(local_gemm_ctxt_type)                         :: local_gemm_ctx = local_gemm_ctxt_type()
      REAL(dp)                                           :: e_gap = 0.0_dp, &
                                                            e_range = 0.0_dp
      LOGICAL                                            :: ri_aux_auto_generated = .FALSE.
   END TYPE mp2_type

   TYPE integ_mat_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: msg
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: sizes
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: indx
      INTEGER                                            :: proc = -1
      TYPE(mp_request_type)                              :: msg_req = mp_request_type()
   END TYPE integ_mat_buffer_type

   TYPE integ_mat_buffer_type_2D
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: msg
      INTEGER                                            :: proc = -1
      TYPE(mp_request_type)                              :: msg_req = mp_request_type()
   END TYPE integ_mat_buffer_type_2D

   TYPE pair_list_type_mp2
      TYPE(pair_list_element_type), DIMENSION(:), &
         ALLOCATABLE                                     :: elements
      INTEGER                                            :: n_element = 0
   END TYPE pair_list_type_mp2

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_release(mp2_env)
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! release the HFX section for the EXX calculation
      IF (.NOT. mp2_env%ri_rpa%reuse_hfx) THEN
         IF (ASSOCIATED(mp2_env%ri_rpa%x_data)) CALL hfx_release(mp2_env%ri_rpa%x_data)
      END IF
      IF (ASSOCIATED(mp2_env%ri_rpa%xc_section_aux)) CALL section_vals_release(mp2_env%ri_rpa%xc_section_aux)
      IF (ASSOCIATED(mp2_env%ri_rpa%xc_section_primary)) CALL section_vals_release(mp2_env%ri_rpa%xc_section_primary)

      IF (mp2_env%eri_method .EQ. do_eri_mme) CALL cp_eri_mme_finalize(mp2_env%eri_mme_param)
      IF (ASSOCIATED(mp2_env%eri_mme_param)) DEALLOCATE (mp2_env%eri_mme_param)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tau_tj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tau_tj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tau_wj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tau_wj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%wj)) DEALLOCATE (mp2_env%ri_rpa_im_time%wj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w)) DEALLOCATE (mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t)) DEALLOCATE (mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t)

      CALL mp2_env%local_gemm_ctx%destroy()

      CALL timestop(handle)

   END SUBROUTINE mp2_env_release

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_create(mp2_env)
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_env_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(mp2_env))

      ALLOCATE (mp2_env)

      NULLIFY (mp2_env%ri_rpa%x_data)

      CALL timestop(handle)

   END SUBROUTINE mp2_env_create

END MODULE mp2_types
