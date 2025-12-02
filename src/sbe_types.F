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
                           kpoint_release
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
      ! dim 1 : spin/kind, dim 2 : rows, dim 3 : cols
      COMPLEX(kind=dp), DIMENSION(:,:,:), ALLOCATABLE    :: overlap_k, &
                                                            inv_overlap_k, &
                                                            ham_zero_k
      REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: overlap_rs, &
                                                            ham_zero_rs
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      INTEGER                                            :: n_ao = 0, &
                                                            n_images = 0, &
                                                            n_spin = 1
   END TYPE sbe_env_type

   PUBLIC :: sbe_env_type, &
             init_sbe_env, &
             release_sbe_env
CONTAINS
   SUBROUTINE init_sbe_env(sbe_env, force_env)
      TYPE(sbe_env_type), POINTER                        :: sbe_env
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN="init_sbe_env"
      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
              POINTER                                    :: sab_all!, &
                                                            !sab_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s_kp, &
                                                            matrix_ks_kp

      CALL timeset(routineN, handle)

      ! Allocation done in this routine, release done in separate routine
      IF (ASSOCIATED(sbe_env)) THEN
         CALL cp_abort(__LOCATION__, "init_sbe_env expects uninitialized pointer to sbe_env.")
      END IF

      ALLOCATE (sbe_env)

      ! Start with kpoints allocation
      CALL force_env_get(force_env, qs_env=qs_env)
      CALL get_qs_env(qs_env, sab_all=sab_all, para_env=para_env, dft_control=dft_control)

      CALL kpoint_create(sbe_env%kpoints)
      CALL kpoint_init_cell_index(sbe_env%kpoints, sab_all, para_env, dft_control)
      sbe_env%n_images = SIZE(sbe_env%kpoints%index_to_cell, 2)

      ! TODO : Construct the (non-symmetric+full) kpoint grid

      ! TODO : Get the kpoint distribution (unless too hard, then just use modulos)

      ! Continue with allocation of rs-matrices (these are defined on full, non-symmetric grid)
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_ks_kp)
      CALL dbcsr_get_info(matrix_s_kp(1,1)%matrix, nfullrows_total=sbe_env%n_ao)
      sbe_env%n_spin = SIZE(matrix_ks_kp, 1)
      ALLOCATE (sbe_env%overlap_rs(1, sbe_env%n_ao, sbe_env%n_ao, sbe_env%n_images), source=0.0_dp)
      ALLOCATE (sbe_env%ham_zero_rs(sbe_env%n_spin, sbe_env%n_ao, sbe_env%n_ao, sbe_env%n_images), source=0.0_dp)

      CALL timestop(handle)
   END SUBROUTINE init_sbe_env
   SUBROUTINE release_sbe_env(sbe_env)
      TYPE(sbe_env_type), POINTER                        :: sbe_env

      CPASSERT(ASSOCIATED(sbe_env))

      CALL kpoint_release(sbe_env%kpoints)
      DEALLOCATE (sbe_env%overlap_rs)
      DEALLOCATE (sbe_env%ham_zero_rs)

      ! Finish by deallocating the pointer itself
      DEALLOCATE (sbe_env)
   END SUBROUTINE release_sbe_env
END MODULE sbe_types
