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
   USE message_passing,                 ONLY: mp_comm_type, &
                                              mp_request_type, &
                                              mp_para_env_type
   USE kpoint_types,                    ONLY: kpoint_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sbe_utils'

   PUBLIC :: print_sbe_header_info, &
             polymorphic_gather

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
! **************************************************************************************************
!> \brief Print basic information about the SBE to be run
!> \param local_vec Locally stored elements of the global vector
!> \param global_vec Destination vector of the stored elements
!> \param para_env Relevant parallel environment
!> \param dist Distribution of vector lengths among the processes
!> \author Stepan Marek
!> \date 10.2026
! **************************************************************************************************
   SUBROUTINE polymorphic_gather(local_vec, global_vec, para_env, dist)
      COMPLEX(kind=dp), DIMENSION(:), INTENT(IN)         :: local_vec
      COMPLEX(kind=dp), DIMENSION(:), INTENT(OUT)        :: global_vec
      TYPE(mp_para_env_type), POINTER                    :: para_env
      ! Determine whether lengths only is better than classic
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist
      CHARACTER(len=*), PARAMETER                        :: routineN="polymorphic_gather"
      TYPE(mp_request_type), DIMENSION(:), ALLOCATABLE   :: requests
      INTEGER                                            :: offset, i, handle

      CALL timeset(routineN, handle)

      offset = 0

      ALLOCATE (requests(para_env%num_pe))

      ! TODO : Might be able to redo via non-separated irecv/isend and waitall
      IF (para_env%mepos == para_env%source) THEN
         ! Receive everything on source
         DO i=1,para_env%num_pe
            IF (i == para_env%mepos+1) THEN
               ! Just copy
               global_vec(offset+1:offset+dist(i)) = local_vec(1:dist(i))
            ELSE
               ! MP ireceive
               CALL para_env%irecv(global_vec(offset+1:offset+dist(i)), & ! message
                  i-1, &! source
                  requests(i)) ! request
            END IF
            offset = offset + dist(i)
         END DO
      ELSE
         ! Send from this rank to source
         CALL para_env%isend(local_vec, & ! message
            para_env%source, & ! destination
            requests(para_env%mepos+1)) ! request
      END IF

      ! Sync
      IF (para_env%mepos == para_env%source) THEN
         DO i=1, para_env%num_pe

            ! We are not waiting for anything on source
            IF (i-1 == para_env%source) CYCLE

            CALL requests(i)%wait()

         END DO
      ELSE
         CALL requests(para_env%mepos+1)%wait()
      ENDIF

      DEALLOCATE (requests)

      CALL timestop(handle)

   END SUBROUTINE polymorphic_gather
END MODULE sbe_utils
