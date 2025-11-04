!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Implements transformations from k-space to R-space for Fortran array matrices
!> \note This code is less performant/more memory consuming than the methods in kpoint_methods.F
!>       Use only when transformations are not the computational bottleneck.
!> \par History
!>       2025.11 created [Stepan Marek]
!> \author Stepan Marek
! **************************************************************************************************
MODULE kpoint_k_R_trafo_simple
   USE kinds,                             ONLY: dp
   USE mathconstants,                     ONLY: twopi
   USE cp_dbcsr_api,                      ONLY: dbcsr_p_type, &
                                                dbcsr_type, &
                                                dbcsr_iterator_type, &
                                                dbcsr_get_info, &
                                                dbcsr_iterator_start, &
                                                dbcsr_iterator_blocks_left, &
                                                dbcsr_iterator_next_block, &
                                                dbcsr_iterator_stop
   USE message_passing,                   ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_k_R_trafo_simple'

   PUBLIC :: replicate_rs_matrices, &
      rs_to_kp

CONTAINS
! **************************************************************************************************
!> \brief Convert dbcsr matrices representing operators in real-space image cells to arrays
!> \param rs_dbcsr Array of dbcsr matrices
!> \param array Multidimensional array - matrices are duplicated on each MPI rank
!> \param index_to_cell Gets the image cell coordinates from the rs_dbcsr index
!> \param cell_to_index Gets the rs_dbcsr index from image cell coordinates
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE replicate_rs_matrices(rs_dbcsr, array, index_to_cell, cell_to_index)
      ! dimension 1 : spin, dimension 2 : image cell index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rs_dbcsr
      TYPE(dbcsr_type)                                   :: nosym_dbcsr
      ! dimension 1 : spin, dimension 2 : image cell index, dimension 3 : rows, dimension 4 : cols
      REAL(kind=dp), DIMENSION(:,:,:,:), POINTER         :: array
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      CHARACTER(len=*), PARAMETER                        :: routineN = "replicate_rs_matrices"
      INTEGER                                            :: n_spin, &
                                                            n_images, &
                                                            handle, &
                                                            j, i, k, l, &
                                                            i_b, j_b, &
                                                            local_rows, local_cols, &
                                                            row_offset, col_offset, &
                                                            i_g, j_g
      TYPE(mp_comm_type)                                 :: group
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL timeset(routineN, handle)
      ! TODO : Sanity checks - at least 1 image exists, created array has the correct dimensions
      ! Get the global number of rows and columns
      n_spin = SIZE(rs_dbcsr, 1)
      n_images = SIZE(rs_dbcsr, 2)
      CALL dbcsr_get_info(rs_dbcsr(1,1)%matrix, group=group)
      ! Iterate over the blocks - get the block coordinates and save to a relevant section locally
      DO j=1,n_spin
         DO i=1,n_images
            CALL dbcsr_iterator_start(iterator, rs_dbcsr(j,i)%matrix)
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
                     array(j, i, i_g, j_g) = block(k, l)
                     IF (i_b /= j_b) THEN
                        ! TODO : Check for mirror image existence
                        array(j, &
                                cell_to_index(-index_to_cell(1,i), -index_to_cell(2,i), -index_to_cell(3,i)), &
                                j_g, i_g) = block(k,l)
                     END IF
                  END DO
               END DO
            END DO
            CALL dbcsr_iterator_stop(iterator)
         END DO
      END DO
      ! Sum only after the entire R-space has been covered on individual ranks
      ! TODO : Spin
      CALL group%sum(array(:,:,:,:))
      CALL timestop(handle)
   END SUBROUTINE replicate_rs_matrices
! **************************************************************************************************
!> \brief Integrate RS matrices (stored as Fortran array) into a kpoint matrix at given kp
!> \param ks_real Array of resulting transformed 
!> \param rs_real Multidimensional array of real parts of the matrix
!> \param array Multidimensional array - matrices are duplicated on each MPI rank
!> \param index_to_cell Gets the image cell coordinates from the rs_dbcsr index
!> \param cell_to_index Gets the rs_dbcsr index from image cell coordinates
!> \param deriv_direction Derivative direction - vector to take scalar product with the gradient
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE rs_to_kp(ks_real, ks_imag, kpoint, rs_real, index_to_cell, deriv_direction)
      ! Index 1 : rows, Index 2 : columns of the matrix (AO or MO indices)
      REAL(kind=dp), DIMENSION(:,:), POINTER, &
         INTENT(OUT)                                     :: ks_real, &
                                                            ks_imag
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: kpoint
      ! Index 1 : image cell index, index 2 : rows, index 3 : columns of the matrix
      REAL(kind=dp), DIMENSION(:,:,:), POINTER, &
         INTENT(IN)                                      :: rs_real
      ! Index 1 : direction, index 2 : image cell index
      INTEGER, DIMENSION(:,:), POINTER                   :: index_to_cell
      REAL(kind=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: deriv_direction
      CHARACTER(len=*), PARAMETER                        :: routineN="rs_to_kp"
      INTEGER                                            :: handle, &
                                                            row_size, &
                                                            col_size, &
                                                            n_images, &
                                                            i, j, k
      REAL(kind=dp)                                      :: cosfactor, &
                                                            sinfactor
      LOGICAL, DIMENSION(3)                              :: derivs_required

      CALL timeset(routineN, handle)
      ! Get matrix constants
      row_size = SIZE(ks_real, 1)
      col_size = SIZE(ks_real, 2)
      n_images = SIZE(rs_real, 1)
      ! Get the required derivatives for the deriv direction
      IF (PRESENT(deriv_direction)) THEN
         DO k=1,3
            derivs_required(k) = (deriv_direction(k) /= 0.0_dp)
         END DO
      ELSE
         derivs_required(:) = .FALSE.
      END IF
      ! TODO : Implement derivatives
      IF (ANY(derivs_required)) THEN
         CPABORT("rs_to_kp kpoint derivatives not yet implemented")
      END IF
      ! Now, iterate over realspace and build the sum
      ! TODO : OMP
      ks_real(:,:) = 0.0_dp
      ks_imag(:,:) = 0.0_dp
      DO i=1,n_images
         ! Construct the phase and factors
         cosfactor = COS(twopi*(kpoint(1) * index_to_cell(1,i) + &
                 kpoint(2) * index_to_cell(2,i) + kpoint(3) * index_to_cell(3,i)))
         sinfactor = SIN(twopi*(kpoint(1) * index_to_cell(1,i) + &
                 kpoint(2) * index_to_cell(2,i) + kpoint(3) * index_to_cell(3,i)))
         DO j=1,row_size
            DO k=1,col_size
               ks_real(j,k) = ks_real(j,k) + cosfactor * rs_real(i, j, k)
               ks_imag(j,k) = ks_imag(j,k) + sinfactor * rs_real(i, j, k)
            END DO
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE rs_to_kp
END MODULE kpoint_k_R_trafo_simple
