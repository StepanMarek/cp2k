!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Construction of the Exchange part of the Fock Matrix
!> \author Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
! **************************************************************************************************
MODULE se_fock_matrix_exchange
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              semi_empirical_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE input_constants,                 ONLY: do_se_IS_kdso,&
                                              do_se_IS_kdso_d
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE multipole_types,                 ONLY: do_multipole_none
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE se_fock_matrix_integrals,        ONLY: dfock2E,&
                                              fock1_2el,&
                                              fock2E
   USE semi_empirical_int_arrays,       ONLY: rij_threshold
   USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
   USE semi_empirical_types,            ONLY: get_se_param,&
                                              se_int_control_type,&
                                              se_taper_type,&
                                              semi_empirical_p_type,&
                                              semi_empirical_type,&
                                              setup_se_int_control_type
   USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                              initialize_se_taper
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_exchange'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: build_fock_matrix_exchange

CONTAINS

! **************************************************************************************************
!> \brief Construction of the Exchange part of the Fock matrix
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_p ...
!> \param calculate_forces ...
!> \param store_int_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE build_fock_matrix_exchange(qs_env, ks_matrix, matrix_p, calculate_forces, &
                                         store_int_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix, matrix_p
      LOGICAL, INTENT(in)                                :: calculate_forces
      TYPE(semi_empirical_si_type), POINTER              :: store_int_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_exchange'

      INTEGER                                            :: atom_a, atom_b, handle, iatom, icol, &
                                                            ikind, integral_screening, irow, &
                                                            jatom, jkind, natorb_a, nkind, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(2)                              :: size_p_block_a
      LOGICAL                                            :: anag, check, defined, found, switch, &
                                                            use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: se_defined
      REAL(KIND=dp)                                      :: delta, dr
      REAL(KIND=dp), DIMENSION(3)                        :: force_ab, rij
      REAL(KIND=dp), DIMENSION(45, 45)                   :: p_block_tot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ks_block_a, ks_block_b, p_block_a, &
                                                            p_block_b
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(se_int_control_type)                          :: se_int_control
      TYPE(se_taper_type), POINTER                       :: se_taper
      TYPE(semi_empirical_control_type), POINTER         :: se_control
      TYPE(semi_empirical_p_type), DIMENSION(:), POINTER :: se_kind_list
      TYPE(semi_empirical_type), POINTER                 :: se_kind_a, se_kind_b
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, cell, force, particle_set, se_control, se_taper)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper, &
                      para_env=para_env, virial=virial)

      CALL initialize_se_taper(se_taper, exchange=.TRUE.)
      se_control => dft_control%qs_control%se_control
      anag = se_control%analytical_gradients
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      nspins = dft_control%nspins

      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(SIZE(ks_matrix) > 0)

      ! Identify proper integral screening (according user requests)
      integral_screening = se_control%integral_screening
      IF ((integral_screening == do_se_IS_kdso_d) .AND. (.NOT. se_control%force_kdsod_EX)) THEN
         integral_screening = do_se_IS_kdso
      END IF
      CALL setup_se_int_control_type(se_int_control, shortrange=.FALSE., &
                                     do_ewald_r3=.FALSE., do_ewald_gks=.FALSE., integral_screening=integral_screening, &
                                     max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)

      CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, &
                      atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)

      nkind = SIZE(atomic_kind_set)
      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, force=force)
         delta = se_control%delta
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)
      END IF

      ALLOCATE (se_defined(nkind), se_kind_list(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), se_parameter=se_kind_a)
         se_kind_list(ikind)%se_param => se_kind_a
         CALL get_se_param(se_kind_a, defined=defined, natorb=natorb_a)
         se_defined(ikind) = (defined .AND. natorb_a >= 1)
      END DO

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)
         IF (.NOT. se_defined(ikind)) CYCLE
         IF (.NOT. se_defined(jkind)) CYCLE
         se_kind_a => se_kind_list(ikind)%se_param
         se_kind_b => se_kind_list(jkind)%se_param

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            switch = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            switch = .TRUE.
         END IF
         ! Retrieve blocks for KS and P
         CALL dbcsr_get_block_p(matrix=ks_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=ks_block_a, found=found)
         CPASSERT(ASSOCIATED(ks_block_a))
         CALL dbcsr_get_block_p(matrix=matrix_p(1)%matrix, &
                                row=irow, col=icol, BLOCK=p_block_a, found=found)
         CPASSERT(ASSOCIATED(p_block_a))
         size_p_block_a(1) = SIZE(p_block_a, 1)
         size_p_block_a(2) = SIZE(p_block_a, 2)
         p_block_tot(1:size_p_block_a(1), 1:size_p_block_a(2)) = 2.0_dp*p_block_a

         ! Handle more configurations
         IF (nspins == 2) THEN
            CALL dbcsr_get_block_p(matrix=ks_matrix(2)%matrix, &
                                   row=irow, col=icol, BLOCK=ks_block_b, found=found)
            CPASSERT(ASSOCIATED(ks_block_b))
            CALL dbcsr_get_block_p(matrix=matrix_p(2)%matrix, &
                                   row=irow, col=icol, BLOCK=p_block_b, found=found)
            CPASSERT(ASSOCIATED(p_block_b))
            check = (size_p_block_a(1) == SIZE(p_block_b, 1)) .AND. (size_p_block_a(2) == SIZE(p_block_b, 2))
            CPASSERT(check)
            p_block_tot(1:SIZE(p_block_a, 1), 1:SIZE(p_block_a, 2)) = p_block_a + p_block_b
         END IF

         dr = DOT_PRODUCT(rij, rij)
         IF (iatom == jatom .AND. dr < rij_threshold) THEN
            ! Once center - Two electron Terms
            IF (nspins == 1) THEN
               CALL fock1_2el(se_kind_a, p_block_tot, p_block_a, ks_block_a, factor=0.5_dp)
            ELSE IF (nspins == 2) THEN
               CALL fock1_2el(se_kind_a, p_block_tot, p_block_a, ks_block_a, factor=1.0_dp)
               CALL fock1_2el(se_kind_a, p_block_tot, p_block_b, ks_block_b, factor=1.0_dp)
            END IF
         ELSE
            ! Exchange Terms
            IF (nspins == 1) THEN
               CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_a, ks_block_a, &
                           factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                           store_int_env=store_int_env)
            ELSE IF (nspins == 2) THEN
               CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_a, ks_block_a, &
                           factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                           store_int_env=store_int_env)

               CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_b, ks_block_b, &
                           factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                           store_int_env=store_int_env)
            END IF
            IF (calculate_forces) THEN
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force_ab = 0.0_dp
               IF (nspins == 1) THEN
                  CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_a, &
                               factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, &
                               delta=delta)
               ELSE IF (nspins == 2) THEN
                  CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_a, &
                               factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, &
                               delta=delta)

                  CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_b, &
                               factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, &
                               delta=delta)
               END IF
               IF (switch) THEN
                  force_ab(1) = -force_ab(1)
                  force_ab(2) = -force_ab(2)
                  force_ab(3) = -force_ab(3)
               END IF
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1.0_dp, force_ab, rij)
               END IF

               force(ikind)%rho_elec(1, atom_a) = force(ikind)%rho_elec(1, atom_a) - force_ab(1)
               force(jkind)%rho_elec(1, atom_b) = force(jkind)%rho_elec(1, atom_b) + force_ab(1)

               force(ikind)%rho_elec(2, atom_a) = force(ikind)%rho_elec(2, atom_a) - force_ab(2)
               force(jkind)%rho_elec(2, atom_b) = force(jkind)%rho_elec(2, atom_b) + force_ab(2)

               force(ikind)%rho_elec(3, atom_a) = force(ikind)%rho_elec(3, atom_a) - force_ab(3)
               force(jkind)%rho_elec(3, atom_b) = force(jkind)%rho_elec(3, atom_b) + force_ab(3)
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (se_kind_list, se_defined)

      CALL finalize_se_taper(se_taper)

      CALL timestop(handle)

   END SUBROUTINE build_fock_matrix_exchange

END MODULE se_fock_matrix_exchange

