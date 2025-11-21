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
   USE kpoint_types, ONLY: kpoint_type
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type
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
   END TYPE sbe_env_type

   PUBLIC :: sbe_env_type
END MODULE sbe_types
