!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Utility functions for SBE code
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
MODULE sbe_utils
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
                                              dbcsr_get_info, &
                                              dbcsr_type, &
                                              dbcsr_desymmetrize, &
                                              dbcsr_get_block_p, &
                                              dbcsr_release, &
                                              dbcsr_iterator_type, &
                                              dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, &
                                              dbcsr_iterator_blocks_left, &
                                              dbcsr_iterator_next_block
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE kpoint_types,                    ONLY: kpoint_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_utils'

   PUBLIC :: print_sbe_header_info, &
           distribute_rs_matrices

CONTAINS
! **************************************************************************************************
!> \brief Print basic information about the SBE to be run
!> \param unit_nr Output unit
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE print_sbe_header_info(unit_nr)
      INTEGER                                           :: unit_nr

      WRITE(unit_nr, "(A80)") " /--------------------------------------"// &
              "---------------------------------------\"
      WRITE(unit_nr, "(A2,A78)") " |", "|"
      WRITE(unit_nr, "(A2,A59,A19)") " |", "Real-time Semiconductor Bloch Equations (SBEs)", "|"
      WRITE(unit_nr, "(A2,A78)") " |", "|"
      WRITE(unit_nr, "(A80)") " \--------------------------------------"// &
              "---------------------------------------/"
   END SUBROUTINE print_sbe_header_info
! TODO : Think whether to do arrays or whether to transfer to special distribution of fm?
! **************************************************************************************************
!> \brief Convert dbcsr matrices representing operators in real-space image cells to arrays
!> \param rs_dbcsr Array of dbcsr matrices
!> \param array Multidimensional array - matrices are duplicated on each MPI rank
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE distribute_rs_matrices(rs_dbcsr, array, index_to_cell, cell_to_index)
      ! dimension 1 : spin, dimension 2 : image cell index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: rs_dbcsr
      TYPE(dbcsr_type)                                  :: nosym_dbcsr
      ! dimension 1 : spin, dimension 2 : image cell index, dimension 3 : rows, dimension 4 : cols
      REAL(kind=dp), DIMENSION(:,:,:,:), POINTER        :: array
      ! TODO : Work with kpoints structure
      INTEGER, DIMENSION(:, :), POINTER                 :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER              :: cell_to_index
      INTEGER                                           :: block_rows, &
                                                           block_cols, &
                                                           i_b, j_b, &
                                                           i, j, &
                                                           k, l, &
                                                           i_g, j_g, &
                                                           local_rows, &
                                                           local_cols, &
                                                           n_images, &
                                                           n_spin, &
                                                           my_prow, &
                                                           my_pcol, &
                                                           row_offset, &
                                                           col_offset
      INTEGER, DIMENSION(:), POINTER                    :: block_row_offset, &
                                                           block_col_offset, &
                                                           proc_row_dist, &
                                                           proc_col_dist
      INTEGER, DIMENSION(2)                             :: dbcsr_shape
      REAL(kind=dp), DIMENSION(:, :), POINTER           :: block
      LOGICAL                                           :: found, owned
      TYPE(mp_comm_type)                                :: group
      TYPE(dbcsr_iterator_type)                         :: iterator

      ! TODO : Sanity checks - at least 1 image exists, created array has the correct dimensions
      ! Get the global number of rows and columns
      dbcsr_shape = SHAPE(rs_dbcsr)
      n_spin = dbcsr_shape(1)
      n_images = dbcsr_shape(2)
      CALL dbcsr_get_info(rs_dbcsr(1,1)%matrix, &
              my_prow=my_prow, my_pcol=my_pcol, group=group)
      ! Iterate over the blocks - get the block coordinates and save to a relevant section locally
      DO i=1,n_images
         CALL dbcsr_iterator_start(iterator, rs_dbcsr(1,i)%matrix)
         DO WHILE(dbcsr_iterator_blocks_left(iterator))
            ! Get the block
            CALL dbcsr_iterator_next_block(iterator, row=i_b, column=j_b, block=block, &
                    row_size=local_rows, col_size=local_cols, &
                    row_offset=row_offset, col_offset=col_offset)
            ! Copy the block to local storage
            DO k=1, local_rows
               DO l=1, local_cols
                  i_g = k + row_offset - 1
                  j_g = l + col_offset - 1
                  ! TODO : OMP, Spin
                  array(1, i, i_g, j_g) = block(k, l)
                  IF (i_b /= j_b) THEN
                     ! TODO : Check for mirror image existence
                     array(1, &
                             cell_to_index(-index_to_cell(1,i), -index_to_cell(2,i), -index_to_cell(3,i)), &
                             j_g, i_g) = block(k,l)
                  END IF
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iterator)
      END DO
      ! Sum only after the entire R-space has been covered on individual ranks
      ! TODO : Spin
      CALL group%sum(array(1,:,:,:))
   END SUBROUTINE distribute_rs_matrices
END MODULE sbe_utils
