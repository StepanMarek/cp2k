!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      Routine to initialize a real space grid from a given input section
!> \par History
!>      01.2014 moved routine from realspace_grid_types into separate file.
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_realspace_grid_init
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE realspace_grid_types,            ONLY: realspace_grid_input_type,&
                                              rsgrid_automatic,&
                                              rsgrid_replicated
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: init_input_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_realspace_grid_init'

CONTAINS

! **************************************************************************************************
!> \brief parses an input section to assign the proper values to the input type
!> \param input_settings ...
!> \param nsmax ...
!> \param rs_grid_section ...
!> \param ilevel ...
!> \param higher_grid_layout the layout of a higher level grid. layouts with
!>       negative or zero values are ignored
!> \par History
!>      01.2008 created [Joost VandeVondele]
!> \note
!>      if rs_grid_section is not present we setup for an replicated setup
! **************************************************************************************************
   SUBROUTINE init_input_type(input_settings, nsmax, rs_grid_section, ilevel, higher_grid_layout)
      TYPE(realspace_grid_input_type), INTENT(OUT)       :: input_settings
      INTEGER, INTENT(IN)                                :: nsmax
      TYPE(section_vals_type), OPTIONAL, POINTER         :: rs_grid_section
      INTEGER, INTENT(IN)                                :: ilevel
      INTEGER, DIMENSION(3), INTENT(IN)                  :: higher_grid_layout

      INTEGER                                            :: isection, max_distributed_level, nsection
      INTEGER, DIMENSION(:), POINTER                     :: tmp

      IF (PRESENT(rs_grid_section)) THEN
         input_settings%nsmax = nsmax
         ! we use the section corresponding to the level, or the largest available one
         ! i.e. the last section defines all following ones
         CALL section_vals_get(rs_grid_section, n_repetition=nsection)
         isection = MAX(1, MIN(ilevel, nsection))
         CALL section_vals_val_get(rs_grid_section, "DISTRIBUTION_TYPE", &
                                   i_rep_section=isection, &
                                   i_val=input_settings%distribution_type)
         CALL section_vals_val_get(rs_grid_section, "DISTRIBUTION_LAYOUT", &
                                   i_rep_section=isection, &
                                   i_vals=tmp)
         input_settings%distribution_layout = tmp
         CALL section_vals_val_get(rs_grid_section, "MEMORY_FACTOR", &
                                   i_rep_section=isection, &
                                   r_val=input_settings%memory_factor)
         CALL section_vals_val_get(rs_grid_section, "HALO_REDUCTION_FACTOR", &
                                   i_rep_section=isection, &
                                   r_val=input_settings%halo_reduction_factor)
         CALL section_vals_val_get(rs_grid_section, "LOCK_DISTRIBUTION", &
                                   i_rep_section=isection, &
                                   l_val=input_settings%lock_distribution)
         CALL section_vals_val_get(rs_grid_section, "MAX_DISTRIBUTED_LEVEL", &
                                   i_rep_section=isection, &
                                   i_val=max_distributed_level)

         ! multigrids that are to coarse are not distributed in the automatic scheme
         IF (input_settings%distribution_type == rsgrid_automatic) THEN
            IF (ilevel > max_distributed_level) THEN
               input_settings%distribution_type = rsgrid_replicated
            END IF
         END IF
      END IF
      IF (input_settings%lock_distribution) THEN
         IF (ALL(higher_grid_layout > 0)) input_settings%distribution_layout = higher_grid_layout
      END IF
   END SUBROUTINE init_input_type

END MODULE cp_realspace_grid_init
