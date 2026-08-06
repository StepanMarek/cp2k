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
   USE kinds,                          ONLY: dp
   USE mathconstants,                  ONLY: twopi
   USE kpoint_types,                   ONLY: kpoint_type, &
                                             get_kpoint_info
   USE qs_environment_types,           ONLY: qs_environment_type, &
                                             get_qs_env
   USE qs_mo_types,                    ONLY: mo_set_type, &
                                             get_mo_set
   USE kpoint_k_r_trafo_simple,        ONLY: rs_to_kp
   USE mathlib,                        ONLY: gemm_square
   USE sbe_types,                      ONLY: sbe_env_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_methods'

   PUBLIC :: init_densmat, &
             nonherm_metric, &
             nonidempot_metric, &
             get_damping_term, &
             get_convective_term, &
             bch_exp_propagate, &
             add_occ_numbers_kp, &
             add_current_kp

CONTAINS
! TODO : For GPU parallelism, probably need to write a kernel instead of a subroutine and
!        handle the kpoint indices on same level as the matrix indices
! **************************************************************************************************
!> \brief Populates the density matrix
!> \param densmat Densitry matrix to be populated (specific spin and kpoint)
!> \param ispin Spin index of the given density matrix (needed to retrieve relevant occupation)
!> \param ham Kp hamiltonian
!> \param smat Kp overlap
!> \param qs_env QS environment, needed to get occupations
!> \param sbe_env SBE environment, provides workspace (densmat_work(1:2,:,:))
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE init_densmat(densmat, ispin, ham, smat, qs_env, sbe_env)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)         :: densmat
      INTEGER, INTENT(IN)                                   :: ispin
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: ham
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: smat
      TYPE(qs_environment_type), POINTER, INTENT(IN)        :: qs_env
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      TYPE(mo_set_type), DIMENSION(:), POINTER              :: mos
      REAL(kind=dp), DIMENSION(:), POINTER                  :: occupations
      INTEGER                                               :: n_ao, &
                                                               i, &
                                                               info

      CALL get_qs_env(qs_env, mos=mos)

      ! TODO : Get rid of sbe_env dependence by providing workspaces?
      ! Diagonalize the hamiltonian with overlap to get the MOs
      n_ao = sbe_env%n_ao
      ! TODO : Sanity checks
      ! Get the occupations from MOS
      CALL get_mo_set(mos(ispin), occupation_numbers=occupations)
      ! Copy input matrices
      sbe_env%densmat_work(1,:,:) = smat(:, :)
      sbe_env%densmat_work(2,:,:) = ham(:, :)
      CALL ZHEGV(1, 'V', 'U', n_ao, sbe_env%densmat_work(2,:,:), n_ao, &
                 sbe_env%densmat_work(1,:,:), n_ao, &
                 sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
      ! ham_copy is now overwritten with eigenvectors
      ! Density mat (in contravariant representation) is given as C * P_occ * C^H
      ! Save copy of eigenvectors
      sbe_env%densmat_work(1,:,:) = sbe_env%densmat_work(2,:,:)
      ! Scale C by P_occ (P_occ * C)
      DO i=1,n_ao
         ! ! DEBUG : Non-equilibrium density - transitions 1<->2
         ! IF (i == 1 .OR. i == 2) THEN
         !    sbe_env%densmat_work(2, :, i) = 0.5_dp * sbe_env%densmat_work(1,:,i) &
         !                                    + 0.5_dp * sbe_env%densmat_work(1,:,3-i)
         ! ELSE
         !    IF (i > SIZE(occupations)) THEN
         !       sbe_env%densmat_work(2,:,i) = 0.0_dp
         !    ELSE
         !       ! Set up according to occupations
         !       ! Scales occupations so that the density is still idempotent even in closed shell case
         !       sbe_env%densmat_work(2,:,i) = sbe_env%densmat_work(1,:,i) * &
         !          REAL(occupations(i), kind=dp) / REAL(sbe_env%n_spin, kind=dp)
         !    END IF
         ! END IF
         IF (i > SIZE(occupations)) THEN
            sbe_env%densmat_work(2,:,i) = 0.0_dp
         ELSE
            ! Set up according to occupations
            ! Scales occupations so that the density is still idempotent even in closed shell case
            ! TODO : Can n_spin be different than 1 or 2?
            sbe_env%densmat_work(2,:,i) = sbe_env%densmat_work(1,:,i) * &
               REAL(occupations(i), kind=dp) / REAL(3-sbe_env%n_spin, kind=dp)
         END IF
      END DO
      ! Multiply to get the densmat
      CALL gemm_square(sbe_env%densmat_work(2,:,:), 'N', &
                       sbe_env%densmat_work(1,:,:), 'C', densmat(:, :))

   END SUBROUTINE init_densmat
! **************************************************************************************************
!> \brief Pretty prints a complex matrix
!> \param mat Matrix to be pretty printed
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE print_mat(mat)
      COMPLEX(kind=dp), DIMENSION(:,:) :: mat

      INTEGER :: i,j,m,n

      m = SIZE(mat,1)
      n = SIZE(mat,2)

      DO i=1,m
         DO j=1,n-1
            WRITE (6,"(E15.3E3,E15.3E3)", advance="no") REAL(mat(i,j)),AIMAG(mat(i,j))
         END DO
            WRITE (6,"(E15.3E3,E15.3E3)") REAL(mat(i,j)),AIMAG(mat(i,j))
      END DO
   END SUBROUTINE print_mat

! **************************************************************************************************
!> \brief Determines the damping term using dephasing time, uses damping workspace
!> \param sbe_env SBE environment with current kp zero field hamiltonian and overlap
!> \param damping_term Complex matrix representing the damping term
!> \param densmat Current kp density matrix, which is needed to construct the damping term
!> \param ispin Spin index of the provided densmat (used to retreat other kpoint matrices)
!> \param ikp Kpoint index of the provided densmat (used to retreat other kpoint matrices)
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE get_damping_term(sbe_env, damping_term, densmat, ispin)
      TYPE(sbe_env_type), POINTER, INTENT(INOUT)            :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: damping_term
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)          :: densmat
      INTEGER, INTENT(IN)                                   :: ispin
      CHARACTER(len=*), PARAMETER                           :: routineN="get_damping_term"
      INTEGER                                               :: handle, &
                                                               info, &
                                                               i

      CALL timeset(routineN, handle)

      IF (sbe_env%dephasing_time == 0.0_dp) THEN
         damping_term(:,:) = 0.0_dp
         CALL timestop(handle)
         RETURN
      END IF

      ! Copy ham_kp and s_kp, in order to not lose them after diagonalization
      sbe_env%damping_work(1,:,:) = sbe_env%ham_zero_kp(ispin,:,:)
      sbe_env%damping_work(2,:,:) = sbe_env%overlap_kp(:,:)
      ! Allocate array for eigenvalues
      ! Diagonalize
      CALL ZHEGV(1, 'V', 'U', sbe_env%n_ao, sbe_env%damping_work(1,:,:), &
                 sbe_env%n_ao, sbe_env%damping_work(2,:,:), sbe_env%n_ao, &
                 sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
      ! With the diagonals, can get rho_MO
      CALL gemm_square(sbe_env%damping_work(1,:,:), 'C', &
                       sbe_env%overlap_kp(:,:), 'N', &
                       densmat, 'N', &
                       sbe_env%damping_work(2,:,:))
      CALL gemm_square(sbe_env%damping_work(2,:,:), 'N', &
                       sbe_env%overlap_kp(:,:), 'N', &
                       sbe_env%damping_work(1,:,:), 'N', &
                       damping_term)
      ! damping_term now contains rho^MO
      ! Set the diagonal to zero - no damping on the diagonal
      ! TODO : OMP Needed?
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(sbe_env,damping_term)
      DO i=1,sbe_env%n_ao
         damping_term(i,i) = 0.0_dp
      END DO
      !$OMP END PARALLEL DO

      ! To finish the damping term, divide by T_2 time, multiply by correct sign
      ! TODO : OMP?
      damping_term(:,:) = -damping_term(:,:) / sbe_env%dephasing_time

      CALL timestop(handle)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Adds to the KS occupation number, for the given kpoint
!> \author Stepan Marek
!> \date 12.2025
! **************************************************************************************************
   SUBROUTINE add_occ_numbers_kp(sbe_env, occ_numbers, densmat, ikp)
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:)                         :: occ_numbers ! dim 1 : spin, dim 2 : band
      COMPLEX(kind=dp), DIMENSION(:,:,:)                    :: densmat ! dim 1 : spin, dim 2 : rows, dim 3 : cols
      INTEGER                                               :: ikp
      INTEGER                                               :: n, info, ispin
      ! ! DEBUG : orthogonality of MOs
      ! INTEGER                                               :: i, j
      ! For current local kpoint, get the transformation from AO to get n_m(k)
      ! Get the molecular orbitals
      DO ispin=1, sbe_env%n_spin
         sbe_env%densmat_work(1,:,:) = sbe_env%ham_zero_kp(ispin,:,:)
         sbe_env%densmat_work(2,:,:) = sbe_env%overlap_kp(:,:)
         CALL ZHEGV(1, 'V', 'U', sbe_env%n_ao, sbe_env%densmat_work(1,:,:), &
                    sbe_env%n_ao, sbe_env%densmat_work(2,:,:), sbe_env%n_ao, &
                    sbe_env%eigvals, sbe_env%eigwork, sbe_env%lwork, sbe_env%eigrwork, info)
         ! With molecular orbitals in densmat_work(1,:,:), we can do the product
         CALL gemm_square(sbe_env%densmat_work(1,:,:), 'C', sbe_env%overlap_kp, 'N', &
                          densmat(ispin,:,:), 'N', sbe_env%densmat_work(2,:,:))
         CALL gemm_square(sbe_env%densmat_work(2,:,:), 'N', sbe_env%overlap_kp, 'N', &
                          sbe_env%densmat_work(1,:,:), 'N', sbe_env%densmat_work(3,:,:))
         ! densmat_work(3) now contains occupations on the diagonal
         DO n=1,SIZE(occ_numbers)
            ! TODO : Are occupation numbers always necceserily real?
            occ_numbers(ispin,n) = occ_numbers(ispin,n) + &
               REAL(sbe_env%densmat_work(3,n,n), kind=dp) * sbe_env%kpoints%wkp(ikp)
         END DO
         ! ! DEBUG : Check precision in orthogonality of the MOs
         ! CALL gemm_square(sbe_env%densmat_work(1,:,:), 'C', sbe_env%overlap_kp, 'N', &
         !    sbe_env%densmat_work(1,:,:), 'N', sbe_env%densmat_work(3,:,:))
         ! DO i=1,sbe_env%n_ao
         !    DO j=1, sbe_env%n_ao
         !       PRINT *, ikp, i, j, sbe_env%densmat_work(3,i,j)
         !    END DO
         ! END DO
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculates the convective term (including the sign). Spinless.
!> \param sbe_env SBE environment with current kp overlap gradient and kp inverse overlap
!> \param convective_term Complex value of the convective term
!> \param densmat Current kp density matrix, which is needed to construct the convective term
!> \param conv_workspace Complex matrix for intermediate term
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE get_convective_term(sbe_env, convective_term, densmat, conv_workspace, field_opt)
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)         :: convective_term, &
                                                               conv_workspace
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)          :: densmat
      INTEGER                                               :: k
      REAL(kind=dp), DIMENSION(3), OPTIONAL                 :: field_opt
      REAL(kind=dp), DIMENSION(3)                           :: field

      IF (PRESENT(field_opt)) THEN
         field(:) = field_opt(:)
      ELSE
         field(:) = sbe_env%field(:)
      END IF

      ! Start by multiplying field with grad of overlap
      conv_workspace(:,:) = CMPLX(-field(1), 0.0,kind=dp) * sbe_env%overlap_grad(1,:,:)
      DO k=2,3
         conv_workspace(:,:) = conv_workspace(:,:) + &
            CMPLX(-field(k), 0.0,kind=dp) * sbe_env%overlap_grad(k,:,:)
      END DO
      ! Now, multiply by S^(-1) and rho to get
      ! - S^(-1) (E dot grad(S)) rho
      CALL gemm_square(sbe_env%inv_overlap_kp, 'N', conv_workspace, 'N', densmat, 'N', &
         convective_term)
   END SUBROUTINE get_convective_term
! **************************************************************************************************
!> \brief Calculates the current contributions from a given k-point as k-derivative of H
!> \note Interprocess k-sum has to be done outside
!> \param sbe_env SBE environment with current RS Hamiltonian
!> \param densmat Current kp density matrix
!> \param kpoint Current kpoint
!> \param ispin Current spin index
!> \param current Current vector contribution from the current kpoint
!> \author Stepan Marek
!> \date 02.2025
! **************************************************************************************************
   SUBROUTINE add_current_kp(sbe_env, densmat, xkp, ispin, current, ikp)
      TYPE(sbe_env_type), POINTER, INTENT(IN)               :: sbe_env
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)       :: densmat
      REAL(kind=dp), DIMENSION(3), INTENT(IN)               :: xkp
      INTEGER, INTENT(IN)                                   :: ispin
      COMPLEX(kind=dp), DIMENSION(3), INTENT(OUT)           :: current
      CHARACTER(len=*), PARAMETER                           :: routineN="add_current_kp"
      INTEGER                                               :: i, &
                                                               handle, &
                                                               k, j, ikp

      CALL timeset(routineN, handle)

      ! Get the k-derivative of the Hamiltonian, store it in the densmat workspace
      ! should be done after propagation?
      ! TODO : Loop for cartesian directions
      ! TODO : This might be a very good place to test the reallocation price
      DO k=1, 3
         CALL rs_to_kp(sbe_env%ham_zero_rs(ispin,:,:,:), sbe_env%densmat_work(1,:,:), &
            sbe_env%kpoints%index_to_cell, xkp, k, sbe_env%hmat)

         ! TODO : OMP Parallelize the trace? Could be unstable...
         DO i=1, SIZE(densmat, 2)
            DO j=1, SIZE(densmat, 1)
               current(k) = current(k) + sbe_env%densmat_work(1,i,j) * densmat(j,i) * &
                  sbe_env%kpoints%wkp(ikp)
            END DO
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE add_current_kp
! **************************************************************************************************
!> \brief Determines the maximum absolute value of element of 2 times non-Hermitian part of a matrix
!> \param mat Matrix whose non-Hermitian part is being investigated
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   FUNCTION nonherm_metric(mat) RESULT(max_buffer)
      COMPLEX(kind=dp), DIMENSION(:,:)                   :: mat

      INTEGER                                            :: i, j, n_ao
      REAL(kind=dp)                                      :: buffer, max_buffer

      n_ao = SIZE(mat, 1)

      max_buffer = 0.0_dp
      buffer = 0.0_dp

      DO i=1, n_ao
         DO j=1,n_ao
            buffer = ABS(mat(i,j) - CONJG(mat(j,i)))
            IF (buffer > max_buffer) THEN
               max_buffer = buffer
            END IF
         END DO
      END DO

   END FUNCTION nonherm_metric
! **************************************************************************************************
!> \brief Determines the non-idempotency metric - max abs value of element of mat * S * mat - mat
!> \param mat Matrix whose non-idempotent part is being investigated
!> \param smat Overlap matrix for the given basis
!> \param work Working area, at least 2 matrices
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   FUNCTION nonidempot_metric(mat, smat, work) RESULT(max_buffer)
      COMPLEX(kind=dp), DIMENSION(:,:)                   :: mat, smat
      COMPLEX(kind=dp), DIMENSION(:,:,:)                 :: work

      INTEGER                                            :: i, j, n_ao
      REAL(kind=dp)                                      :: buffer, max_buffer

      n_ao = SIZE(mat, 1)

      max_buffer = 0.0_dp
      buffer = 0.0_dp

      ! Multiply S * rho * S * rho * S
      ! S * rho
      CALL gemm_square(smat, "N", mat, "N", work(1,:,:))
      ! (S * rho) * S
      CALL gemm_square(work(1,:,:), "N", smat, "N", work(2,:,:))
      ! (S * rho * S ) * rho
      CALL gemm_square(work(2,:,:), "N", mat, "N", work(1,:,:))
      ! (S * rho * S * rho ) * S
      CALL gemm_square(work(1,:,:), "N", smat, "N", work(2,:,:))

      ! Multiply S * rho * S
      CALL gemm_square(smat, "N", mat, "N", smat, "N", work(1,:,:))

      DO i=1, n_ao
         DO j=1,n_ao
            buffer = ABS(work(2,i,j) - work(1,i,j))
            IF (buffer > max_buffer) THEN
               max_buffer = buffer
            END IF
         END DO
      END DO

   END FUNCTION nonidempot_metric
   FUNCTION abs_metric(mat) RESULT(max_buffer)
      COMPLEX(kind=dp), DIMENSION(:,:)                      :: mat
      REAL(kind=dp)                                         :: max_buffer

      INTEGER                                               :: i, j, n_ao
      REAL(kind=dp)                                         :: buffer

      n_ao = SIZE(mat, 1)

      max_buffer = 0.0_dp
      buffer = 0.0_dp

      DO i=1, n_ao
         DO j=1, n_ao
            buffer = ABS(mat(i,j))
            IF (buffer > max_buffer) max_buffer = buffer
         END DO
      END DO
   END FUNCTION abs_metric
! **************************************************************************************************
!> \brief Produces the effect of the exponential propagation of densmat using hamiltonian and smat
!> \note Does not assume Hamiltonian Hermiticity, works with Fortran native arrays
!> \author Stepan Marek
!> \date 07.2026
! **************************************************************************************************
   SUBROUTINE bch_exp_propagate(densmat_new, densmat_old, inv_overlap, hamiltonian, timestep, eps, &
                                converged, workspace, max_iter_opt)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)         :: densmat_new
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(IN)          :: densmat_old, &
                                                               inv_overlap, &
                                                               hamiltonian
      REAL(kind=dp), INTENT(IN)                             :: timestep, eps
      INTEGER, INTENT(OUT)                                  :: converged
      !COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: workspace_opt ! Save some time on reallocations
      COMPLEX(kind=dp), DIMENSION(:,:,:)                    :: workspace ! Save some time on reallocations
      INTEGER, INTENT(IN), OPTIONAL                         :: max_iter_opt
      CHARACTER(len=*), PARAMETER                           :: routineN="bch_exp_propagate"
      INTEGER                                               :: handle, &
                                                               n_ao, &
                                                               max_iter, &
                                                               i
      ! COMPLEX(kind=dp), DIMENSION(:,:,:), POINTER           :: workspace
      COMPLEX(kind=dp)                                      :: prefactor

      CALL timeset(routineN, handle)

      n_ao = SIZE(densmat_old, 1)

      converged = 1

      max_iter = 30
      IF (PRESENT(max_iter_opt)) max_iter = max_iter_opt

      ! IF (PRESENT(workspace_opt)) THEN
      !    workspace => workspace_opt
      ! ELSE
      !    NULLIFY (workspace)
      !    ALLOCATE (workspace(4, n_ao, n_ao))
      ! END IF

      ! Store the cumulant in workspace(1,:,:),
      ! previous commutator in workspace(2,:,:),
      ! intermediary results in workspace(3,:,:) and workspace(4,:,:)
 
      ! Initialize cummulant
      ! C_0 (S^-1 H , H^dagger S^-1; rho)
      workspace(1,:,:) = densmat_old(:,:)
      ! Initialize commutator
      workspace(2,:,:) = densmat_old(:,:)
      prefactor = CMPLX(1.0, 0.0, kind=dp)

      DO i=1, max_iter
         ! Determine C_i(S^-1 H^dagger, H S^-1; rho)
         ! as C_1(S^-1 H^dagger, H S^-1); C_(i-1)(S^-1 H^dagger , H S^-1; rho)) =
         ! = (S^-1 H^dagger C_(i-1) -  C_(i-1) H S^-1)
         CALL gemm_square(inv_overlap, "N", hamiltonian, "C", workspace(2,:,:), "N", workspace(3,:,:))
         CALL gemm_square(workspace(2,:,:), "N", hamiltonian, "N", inv_overlap, "N", workspace(4,:,:))

         ! Update the commutator
         workspace(2,:,:) = workspace(3,:,:) - workspace(4,:,:)
         ! Update the cummulant
         prefactor = prefactor * CMPLX(0.0, -timestep/REAL(i, kind=dp), kind=dp)
         workspace(1,:,:) = workspace(1,:,:) + prefactor * workspace(2,:,:)

         ! If the added term is below the convergence threshold, exit
         IF (abs_metric(prefactor * workspace(2,:,:)) < eps) THEN
            converged = 0
            EXIT
         END IF
      END DO

      ! Always copy the best estimate to new matrix, even if not converged
      densmat_new(:,:) = workspace(1,:,:)

      ! IF (.NOT. PRESENT(workspace_opt)) DEALLOCATE (workspace)

      CALL timestop(handle)
   END SUBROUTINE bch_exp_propagate
END MODULE sbe_methods
