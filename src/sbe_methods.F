!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods for the SBE code (more model based, as opposed to sbe_utils, which are code based)
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
MODULE sbe_methods
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE kpoint_types,                    ONLY: kpoint_type, &
                                              get_kpoint_info

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_methods'

CONTAINS
! TODO : For GPU parallelism, probably need to write a kernel instead of a subroutine and
!        handle the kpoint indices on same level as the matrix indices
END MODULE sbe_methods
