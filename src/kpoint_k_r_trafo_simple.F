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
MODULE kpoint_k_r_trafo_simple
   USE kinds,                             ONLY: dp
   USE mathconstants,                     ONLY: twopi
   USE cp_dbcsr_api,                      ONLY: dbcsr_p_type, &
                                                dbcsr_type, &
                                                dbcsr_get_info, &
                                                dbcsr_get_block_p
   USE message_passing,                   ONLY: mp_comm_type
   USE qs_neighbor_list_types,            ONLY: neighbor_list_set_p_type, &
                                                neighbor_list_iterator_p_type, &
                                                get_iterator_info, &
                                                neighbor_list_iterate, &
                                                neighbor_list_iterator_create, &
                                                neighbor_list_iterator_release
   USE kpoint_types,                      ONLY: kpoint_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_k_r_trafo_simple'

   PUBLIC :: replicate_rs_matrices, &
      rs_to_kp

CONTAINS
! **************************************************************************************************
!> \brief Convert dbcsr matrices representing operators in real-space image cells to arrays
!> \param rs_dbcsr_src Array of dbcsr matrices
!> \param kpoint_src The kpoint environment of the source matrix (providing neighbor list and cell_to_index)
!> \param rs_array_dst Multidimensional array - matrices are duplicated on each MPI rank
!> \param cell_to_index_dst Cell to index array for the destination array
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE replicate_rs_matrices(rs_dbcsr_src, kpoint_src, rs_array_dst, cell_to_index_dst)
      ! dimension 1 : spin, dimension 2 : image cell index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rs_dbcsr_src
      TYPE(kpoint_type), POINTER                         :: kpoint_src
      ! dimension 1 : spin/derivs, dimension 2 : image cell index, dimension 3 : rows, dimension 4 : cols
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: rs_array_dst
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_dst
      CHARACTER(len=*), PARAMETER                        :: routineN = "replicate_rs_matrices"
      INTEGER                                            :: n_spin, ispin, &
                                                            handle, &
                                                            i, j, &
                                                            row_offset, col_offset, &
                                                            iatom, jatom, &
                                                            dst_index, src_index, &
                                                            col_size, row_size
      INTEGER, DIMENSION(:), POINTER                     :: row_offsets, col_offsets
      TYPE(mp_comm_type)                                 :: group
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: found
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp_src

      CALL timeset(routineN, handle)

      ! Start by constructing the cell_to_index_src
      sab_kp_src => kpoint_src%sab_nl
      ! TODO : Sanity checks - array to be filled has the correct dimensions
      ! NOTE : The first index in matrix_s_kp is not spin index, but number of derivatives.
      !        But, for matrix_ks_kp, this is indeed the spin index.
      n_spin = SIZE(rs_dbcsr_src, 1)
      IF (SIZE(rs_dbcsr_src, 2) < 1) THEN
         CALL cp_abort(__LOCATION__, "No source image cells provided!")
      END IF
      CALL dbcsr_get_info(rs_dbcsr_src(1, 1)%matrix, group=group)
      DO ispin=1,n_spin
         CALL neighbor_list_iterator_create(iterator, sab_kp_src)
         DO WHILE (neighbor_list_iterate(iterator) == 0)
            CALL get_iterator_info(iterator, cell=cell, iatom=iatom, jatom=jatom)
            src_index = kpoint_src%cell_to_index(cell(1), cell(2), cell(3))
            IF (src_index == 0) THEN
               CALL cp_abort(__LOCATION__, "Image not found in the source array.")
            END IF
            CALL dbcsr_get_info(rs_dbcsr_src(ispin, src_index)%matrix, &
                                row_blk_offset=row_offsets, col_blk_offset=col_offsets)
            ! Get the source block
            IF (iatom > jatom) THEN
               ! More complicated swaps
               CALL dbcsr_get_block_p(rs_dbcsr_src(ispin,src_index)%matrix, &
                                      row=jatom, col=iatom, &
                                      block=block, found=found, &
                                      row_size=row_size, col_size=col_size)
               row_offset = row_offsets(jatom)
               col_offset = col_offsets(iatom)
               ! If not found, some other rank will handle this
               IF (.NOT. found) CYCLE
               dst_index = cell_to_index_dst(-cell(1), -cell(2), -cell(3))
               IF (dst_index == 0) THEN
                  CALL cp_abort(__LOCATION__, "Mirror image not found in destination array.")
               END IF
               DO i = 1,row_size
                  DO j = 1,col_size
                     rs_array_dst(ispin, dst_index, i+row_offset-1, j+col_offset-1) = block(i, j)
                  END DO
               END DO
               dst_index = cell_to_index_dst(cell(1), cell(2), cell(3))
               IF (dst_index == 0) THEN
                  CALL cp_abort(__LOCATION__, "Image not found in destination array.")
               END IF
               ! Save the transposed part
               DO i = 1,row_size
                  DO j = 1,col_size
                     rs_array_dst(ispin, dst_index, j+col_offset-1, i+row_offset-1) = block(i, j)
                  END DO
               END DO
            ELSE
               ! Standard case, just save the block directly and if off-diagonal, save the transpose to negative
               CALL dbcsr_get_block_p(rs_dbcsr_src(ispin, src_index)%matrix, &
                                      row=iatom, col=jatom, block=block, found=found,&
                                      row_size=row_size, col_size=col_size)
               row_offset = row_offsets(iatom)
               col_offset = col_offsets(jatom)
               IF (.NOT. found) CYCLE
               dst_index = cell_to_index_dst(cell(1), cell(2), cell(3))
               IF (dst_index == 0) THEN
                  CALL cp_abort(__LOCATION__, "Image not found in destination array.")
               END IF
               DO i=1, row_size
                  DO j=1, col_size
                     rs_array_dst(ispin, dst_index, i+row_offset-1, j+col_offset-1) = block(i, j)
                  END DO
               END DO
               IF ( iatom /= jatom ) THEN
                  ! Also save the transpose
                  dst_index = cell_to_index_dst(-cell(1), -cell(2), -cell(3))
                  IF (dst_index == 0) THEN
                     CALL cp_abort(__LOCATION__, "Mirror image not found in destination array.")
                  END IF
                  DO i=1, row_size
                     DO j=1, col_size
                        rs_array_dst(ispin, dst_index, j+col_offset-1, i+row_offset-1) = block(i, j)
                     END DO
                  END DO
               END IF
            END IF
         END DO
         CALL neighbor_list_iterator_release(iterator)
      END DO
      CALL group%sum(rs_array_dst(:,:,:,:))
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
   SUBROUTINE rs_to_kp(ks_real, ks_imag, xkp, rs_real, index_to_cell, deriv_direction, hmat)
      ! Index 1 : rows, Index 2 : columns of the matrix (AO or MO indices)
      REAL(kind=dp), DIMENSION(:,:), INTENT(OUT)         :: ks_real, &
                                                            ks_imag
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: xkp
      ! Index 1 : image cell index, index 2 : rows, index 3 : columns of the matrix
      REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN)        :: rs_real
      ! Index 1 : direction, index 2 : image cell index
      INTEGER, DIMENSION(:,:), POINTER                   :: index_to_cell
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(kind=dp), DIMENSION(3)                        :: cell_vectors
      REAL(kind=dp)                                      :: deriv_factor
      REAL(kind=dp), DIMENSION(3,3), INTENT(IN), OPTIONAL:: hmat
      CHARACTER(len=*), PARAMETER                        :: routineN="rs_to_kp"
      INTEGER                                            :: handle, &
                                                            row_size, &
                                                            col_size, &
                                                            n_images, &
                                                            i, j, k
      REAL(kind=dp)                                      :: cosfactor, &
                                                            sinfactor
      LOGICAL                                            :: derivs_required

      CALL timeset(routineN, handle)
      ! Get matrix constants
      row_size = SIZE(ks_real, 1)
      col_size = SIZE(ks_real, 2)
      n_images = SIZE(rs_real, 1)
      ! Get the required derivatives for the deriv direction
      derivs_required = .FALSE.
      IF (PRESENT(deriv_direction)) THEN
         derivs_required = (deriv_direction /= 0)
         IF (.NOT. PRESENT(hmat)) THEN
            CALL cp_abort(__LOCATION__, "derivative requested but h matrix not provided")
         END IF
      END IF
      ! Now, iterate over realspace and build the sum
      ! TODO : OMP
      ks_real(:,:) = 0.0_dp
      ks_imag(:,:) = 0.0_dp
      DO i=1,n_images
         ! Construct the phase and factors
         deriv_factor = 1.0_dp
         IF (derivs_required) THEN
            cell_vectors = MATMUL(hmat, REAL(index_to_cell(1:3, i), dp))
            deriv_factor = cell_vectors(deriv_direction)
         END IF
         cosfactor = COS(twopi * SUM(xkp * index_to_cell(:,i)))
         sinfactor = SIN(twopi * SUM(xkp * index_to_cell(:,i)))
         DO j=1,row_size
            DO k=1,col_size
               ks_real(j,k) = ks_real(j,k) + deriv_factor * cosfactor * rs_real(i, j, k)
               ks_imag(j,k) = ks_imag(j,k) + deriv_factor * sinfactor * rs_real(i, j, k)
            END DO
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE rs_to_kp
END MODULE kpoint_k_r_trafo_simple
