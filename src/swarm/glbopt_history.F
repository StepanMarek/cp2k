!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief History of minima, calculates, stores and compares fingerprints of minima.
!>        Used by Minima Hopping and Minima Crawling.
!> \author Ole Schuett
! **************************************************************************************************
MODULE glbopt_history
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE history_fingerprint_type
      PRIVATE
      REAL(KIND=dp)                            :: Epot = 0.0
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: goedecker
   END TYPE history_fingerprint_type

   TYPE history_entry_type
      TYPE(history_fingerprint_type), POINTER :: p => Null()
      INTEGER                                 :: id = -1
   END TYPE history_entry_type

   TYPE history_type
      PRIVATE
      TYPE(history_entry_type), DIMENSION(:), POINTER :: entries => Null()
      INTEGER                              :: length = 0
      INTEGER                              :: iw = -1
      REAL(KIND=dp)                        :: E_precision = 0.0
      REAL(KIND=dp)                        :: FP_precision = 0.0
   END TYPE history_type

   PUBLIC :: history_type, history_fingerprint_type
   PUBLIC :: history_init, history_finalize
   PUBLIC :: history_add, history_lookup
   PUBLIC :: history_fingerprint
   PUBLIC :: history_fingerprint_match

   LOGICAL, PARAMETER                     :: debug = .FALSE.
   INTEGER, PARAMETER                     :: history_grow_unit = 1000
CONTAINS

! **************************************************************************************************
!> \brief Initializes a history.
!> \param history ...
!> \param history_section ...
!> \param iw ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE history_init(history, history_section, iw)
      TYPE(history_type), INTENT(INOUT)                  :: history
      TYPE(section_vals_type), POINTER                   :: history_section
      INTEGER                                            :: iw

      ALLOCATE (history%entries(history_grow_unit))
      history%iw = iw
      CALL section_vals_val_get(history_section, "ENERGY_PRECISION", &
                                r_val=history%E_precision)
      CALL section_vals_val_get(history_section, "FINGERPRINT_PRECISION", &
                                r_val=history%FP_precision)

      IF (iw > 0) THEN
         WRITE (iw, '(A,T66,E15.3)') &
            " GLBOPT| History energy precision", history%E_precision
         WRITE (iw, '(A,T66,E15.3)') &
            " GLBOPT| History fingerprint precision", history%FP_precision
      END IF
   END SUBROUTINE history_init

! **************************************************************************************************
!> \brief Calculates a fingerprint for a given configuration.
!> \param Epot ...
!> \param pos ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION history_fingerprint(Epot, pos) RESULT(fp)
      REAL(KIND=dp), INTENT(IN)                          :: Epot
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: pos
      TYPE(history_fingerprint_type)                     :: fp

      INTEGER                                            :: handle
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp

      CALL timeset("glbopt_history_fingerprint", handle)

      NULLIFY (tmp)
      fp%Epot = Epot
      CALL goedecker_fingerprint(pos, tmp)

      !copy pointer to allocatable
      ALLOCATE (fp%goedecker(SIZE(tmp)))
      fp%goedecker(:) = tmp
      DEALLOCATE (tmp)

      CALL timestop(handle)
   END FUNCTION history_fingerprint

! **************************************************************************************************
!> \brief Helper routine for history_fingerprint.
!>        Calculates a fingerprint based on inter-atomic distances.
!> \param pos ...
!> \param res ...
!> \author Stefan Goedecker
! **************************************************************************************************
   SUBROUTINE goedecker_fingerprint(pos, res)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: pos
      REAL(KIND=dp), DIMENSION(:), POINTER               :: res

      INTEGER                                            :: i, info, j, N
      REAL(KIND=dp)                                      :: d2, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: matrix, work
      REAL(KIND=dp), DIMENSION(3)                        :: d

      IF (ASSOCIATED(res)) CPABORT("goedecker_fingerprint: res already allocated")
      N = SIZE(pos)/3 ! number of atoms

      ALLOCATE (matrix(N, N), work(N, N))
      DO i = 1, N
         matrix(i, i) = 1.0
         DO j = i + 1, N
            d = pos(3*i - 2:3*i) - pos(3*j - 2:3*j)
            d2 = SUM(d**2)
            t = EXP(-0.5*d2)
            matrix(i, j) = t
            matrix(j, i) = t
         END DO
      END DO
      ALLOCATE (res(N))
      ! matrix values are garbage on exit because of jobz='N'
      CALL dsyev('N', 'U', N, matrix, N, res, work, N**2, info)
      IF (info /= 0) CPABORT("goedecker_fingerprint: DSYEV failed")
   END SUBROUTINE goedecker_fingerprint

! **************************************************************************************************
!> \brief Checks if two given fingerprints match.
!> \param history ...
!> \param fp1 ...
!> \param fp2 ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION history_fingerprint_match(history, fp1, fp2) RESULT(res)
      TYPE(history_type), INTENT(IN)                     :: history
      TYPE(history_fingerprint_type), INTENT(IN)         :: fp1, fp2
      LOGICAL                                            :: res

      res = (ABS(fp1%Epot - fp2%Epot) < history%E_precision) .AND. &
            (fingerprint_distance(fp1, fp2) < history%fp_precision)

   END FUNCTION history_fingerprint_match

! **************************************************************************************************
!> \brief Helper routine for history_fingerprint_match
!>        Calculates the distance between two given fingerprints.
!> \param fp1 ...
!> \param fp2 ...
!> \return ...
!> \author Stefan Goedecker
! **************************************************************************************************
   PURE FUNCTION fingerprint_distance(fp1, fp2) RESULT(res)
      TYPE(history_fingerprint_type), INTENT(IN)         :: fp1, fp2
      REAL(KIND=dp)                                      :: res

      res = SQRT(SUM((fp1%goedecker - fp2%goedecker)**2)/SIZE(fp1%goedecker))
   END FUNCTION fingerprint_distance

! **************************************************************************************************
!> \brief Addes a new fingerprints to the history.
!>        Optionally, an abitrary id can be stored alongside the fingerprint.
!> \param history ...
!> \param fingerprint ...
!> \param id ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE history_add(history, fingerprint, id)
      TYPE(history_type), INTENT(INOUT)                  :: history
      TYPE(history_fingerprint_type), INTENT(IN)         :: fingerprint
      INTEGER, INTENT(IN), OPTIONAL                      :: id

      INTEGER                                            :: handle, i, k, n
      TYPE(history_entry_type), DIMENSION(:), POINTER    :: tmp

      CALL timeset("glbopt_history_add", handle)

      n = SIZE(history%entries)
      IF (n == history%length) THEN
         ! grow history%entries array
         tmp => history%entries
         ALLOCATE (history%entries(n + history_grow_unit))
         history%entries(1:n) = tmp(:)
         DEALLOCATE (tmp)
         n = n + history_grow_unit
      END IF

      k = interpolation_search(history, fingerprint%Epot)

      !history%entries(k+1:) = history%entries(k:n-1)
      !Workaround for an XLF bug - pointer array copy does
      !not work correctly
      DO i = n, k + 1, -1
         history%entries(i) = history%entries(i - 1)
      END DO

      ALLOCATE (history%entries(k)%p)
      history%entries(k)%p = fingerprint
      IF (PRESENT(id)) &
         history%entries(k)%id = id
      history%length = history%length + 1

      IF (debug) THEN
         ! check history for correct order
         DO k = 1, history%length
            !WRITE(*,*) "history: ", k, "Epot",history%entries(k)%p%Epot
            IF (k > 1) THEN
               IF (history%entries(k - 1)%p%Epot > history%entries(k)%p%Epot) &
                  CPABORT("history_add: history in wrong order")
            END IF
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE history_add

! **************************************************************************************************
!> \brief Checks if a given fingerprints is contained in the history.
!> \param history ...
!> \param fingerprint ...
!> \param found ...
!> \param id ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE history_lookup(history, fingerprint, found, id)
      TYPE(history_type), INTENT(IN)                     :: history
      TYPE(history_fingerprint_type), INTENT(IN)         :: fingerprint
      LOGICAL, INTENT(OUT)                               :: found
      INTEGER, INTENT(OUT), OPTIONAL                     :: id

      INTEGER                                            :: found_i, handle, i, k, k_max, k_min
      REAL(KIND=dp)                                      :: best_match, dist, Epot

      CALL timeset("glbopt_history_lookup", handle)

      found = .FALSE.
      IF (PRESENT(id)) id = -1
      best_match = HUGE(1.0_dp)

      IF (history%length > 0) THEN
         Epot = fingerprint%Epot
         k = interpolation_search(history, fingerprint%Epot)

         DO k_min = k - 1, 1, -1
            IF (history%entries(k_min)%p%Epot < Epot - history%E_precision) EXIT
         END DO

         DO k_max = k, history%length
            IF (history%entries(k_max)%p%Epot > Epot + history%E_precision) EXIT
         END DO

         k_min = MAX(k_min + 1, 1)
         k_max = MIN(k_max - 1, history%length)

         IF (debug) found_i = -1

         DO i = k_min, k_max
            dist = fingerprint_distance(fingerprint, history%entries(i)%p)
            !WRITE(*,*) "entry ", i, " dist: ",dist
            IF (dist < history%fp_precision .AND. dist < best_match) THEN
               best_match = dist
               found = .TRUE.
               IF (PRESENT(id)) id = history%entries(i)%id
               IF (debug) found_i = i
            END IF
         END DO

         IF (debug) CALL verify_history_lookup(history, fingerprint, found_i)
      END IF

      CALL timestop(handle)

   END SUBROUTINE history_lookup

! **************************************************************************************************
!> \brief Helper routine for history_lookup
!> \param history ...
!> \param Efind ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION interpolation_search(history, Efind) RESULT(res)
      TYPE(history_type), INTENT(IN)                     :: history
      REAL(KIND=dp), INTENT(IN)                          :: Efind
      INTEGER                                            :: res

      INTEGER                                            :: high, low, mid
      REAL(KIND=dp)                                      :: slope

      low = 1
      high = history%length

      DO WHILE (low < high)
         !linear interpolation
         slope = REAL(high - low, KIND=dp)/(history%entries(high)%p%Epot - history%entries(low)%p%Epot)
         mid = low + INT(slope*(Efind - history%entries(low)%p%Epot))
         mid = MIN(MAX(mid, low), high)

         IF (history%entries(mid)%p%Epot < Efind) THEN
            low = mid + 1
         ELSE
            high = mid - 1
         END IF
      END DO

      IF (0 < low .AND. low <= history%length) THEN
         IF (Efind > history%entries(low)%p%Epot) low = low + 1
      END IF

      res = low
   END FUNCTION interpolation_search

! **************************************************************************************************
!> \brief Debugging routine, performs a slow (but robust) linear search.
!> \param history ...
!> \param fingerprint ...
!> \param found_i_ref ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE verify_history_lookup(history, fingerprint, found_i_ref)
      TYPE(history_type), INTENT(IN)                     :: history
      TYPE(history_fingerprint_type), INTENT(IN)         :: fingerprint
      INTEGER, INTENT(IN)                                :: found_i_ref

      INTEGER                                            :: found_i, i
      REAL(KIND=dp)                                      :: best_fp_match, Epot_dist, fp_dist

      found_i = -1
      best_fp_match = HUGE(1.0_dp)

      DO i = 1, history%length
         Epot_dist = ABS(fingerprint%Epot - history%entries(i)%p%Epot)
         IF (Epot_dist > history%E_precision) CYCLE
         fp_dist = fingerprint_distance(fingerprint, history%entries(i)%p)
         !WRITE(*,*) "entry ", i, " dist: ",dist
         IF (fp_dist < history%fp_precision .AND. fp_dist < best_fp_match) THEN
            best_fp_match = fp_dist
            found_i = i
         END IF
      END DO

      IF (found_i /= found_i_ref) THEN
         WRITE (*, *) found_i, found_i_ref
         CPABORT("verify_history_lookup failed")
      END IF

   END SUBROUTINE verify_history_lookup

! **************************************************************************************************
!> \brief Finalizes a history.
!> \param history ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE history_finalize(history)
      TYPE(history_type)                                 :: history

      INTEGER                                            :: i

      DO i = 1, history%length
         IF (ASSOCIATED(history%entries(i)%p)) &
            DEALLOCATE (history%entries(i)%p)
      END DO

      DEALLOCATE (history%entries)

   END SUBROUTINE history_finalize

END MODULE glbopt_history
