!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculate optx
!> \note
!>      will need proper testing / review
!> \author Joost VandeVondele [03.2004]
! **************************************************************************************************
MODULE xc_optx
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE xc_derivative_desc,              ONLY: deriv_norm_drho,&
                                              deriv_norm_drhoa,&
                                              deriv_norm_drhob,&
                                              deriv_rho,&
                                              deriv_rhoa,&
                                              deriv_rhob
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_get_derivative
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_optx'

   PUBLIC :: optx_lda_info, optx_lda_eval, optx_lsd_info, optx_lsd_eval
CONTAINS

! **************************************************************************************************
!> \brief info about the optx functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv implemented derivative of the xc functional
!> \author Joost
! **************************************************************************************************
   SUBROUTINE optx_lda_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "OPTX, Handy NC and Cohen AJ,  JCP 116, p. 5411 (2002) (LDA)"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "OPTX exchange (LDA)"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         needs%norm_drho = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1
   END SUBROUTINE optx_lda_info

! **************************************************************************************************
!> \brief info about the optx functional (LSD)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv implemented derivative of the xc functional
!> \author Joost
! **************************************************************************************************
   SUBROUTINE optx_lsd_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "OPTX, Handy NC and Cohen AJ,  JCP 116, p. 5411 (2002), (LSD) "
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "OPTX exchange (LSD)"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho_spin = .TRUE.
         needs%norm_drho_spin = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1
   END SUBROUTINE optx_lsd_info

! **************************************************************************************************
!> \brief evaluates the optx functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param optx_params input parameter (scaling)
!> \par History
!>      01.2007 added scaling [Manuel Guidon]
!> \author Joost
! **************************************************************************************************
   SUBROUTINE optx_lda_eval(rho_set, deriv_set, grad_deriv, optx_params)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: grad_deriv
      TYPE(section_vals_type), POINTER                   :: optx_params

      INTEGER                                            :: npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: a1, a2, epsilon_drho, epsilon_rho, gam, &
                                                            sx
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: e_0, e_ndrho, e_rho, norm_drho, rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      NULLIFY (e_0, e_ndrho, e_rho, norm_drho, rho)

      CALL section_vals_val_get(optx_params, "scale_x", r_val=sx)
      CALL section_vals_val_get(optx_params, "a1", r_val=a1)
      CALL section_vals_val_get(optx_params, "a2", r_val=a2)
      CALL section_vals_val_get(optx_params, "gamma", r_val=gam)

      CALL xc_rho_set_get(rho_set, rho=rho, &
                          norm_drho=norm_drho, local_bounds=bo, rho_cutoff=epsilon_rho, &
                          drho_cutoff=epsilon_drho)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_0)
      deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_rho)
      deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
      IF (grad_deriv > 1 .OR. grad_deriv < -1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

      CALL optx_lda_calc(rho=rho, norm_drho=norm_drho, &
                         e_0=e_0, e_rho=e_rho, e_ndrho=e_ndrho, &
                         npoints=npoints, epsilon_rho=epsilon_rho, &
                         epsilon_drho=epsilon_drho, sx=sx, &
                         a1=a1, a2=a2, gam=gam)
   END SUBROUTINE optx_lda_eval

! **************************************************************************************************
!> \brief evaluates the optx functional for lsd
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param optx_params input parameter (scaling)
!> \par History
!>      01.2007 added scaling [Manuel Guidon]
!> \author Joost
! **************************************************************************************************
   SUBROUTINE optx_lsd_eval(rho_set, deriv_set, grad_deriv, optx_params)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: grad_deriv
      TYPE(section_vals_type), POINTER                   :: optx_params

      INTEGER                                            :: ispin, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: a1, a2, epsilon_drho, epsilon_rho, gam, &
                                                            sx
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: e_0
      TYPE(cp_3d_r_cp_type), DIMENSION(2)                :: e_ndrho, e_rho, ndrho, rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      NULLIFY (e_0)
      DO ispin = 1, 2
         NULLIFY (e_rho(ispin)%array)
         NULLIFY (e_ndrho(ispin)%array)
         NULLIFY (rho(ispin)%array)
         NULLIFY (ndrho(ispin)%array)
      END DO

      CALL section_vals_val_get(optx_params, "scale_x", r_val=sx)
      CALL section_vals_val_get(optx_params, "a1", r_val=a1)
      CALL section_vals_val_get(optx_params, "a2", r_val=a2)
      CALL section_vals_val_get(optx_params, "gamma", r_val=gam)

      CALL xc_rho_set_get(rho_set, rhoa=rho(1)%array, rhob=rho(2)%array, &
                          norm_drhoa=ndrho(1)%array, &
                          norm_drhob=ndrho(2)%array, rho_cutoff=epsilon_rho, &
                          drho_cutoff=epsilon_drho, local_bounds=bo)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_0)
      deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_rho(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_rho(2)%array)

      deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_ndrho(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_ndrho(2)%array)

      IF (grad_deriv > 1 .OR. grad_deriv < -1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF
      DO ispin = 1, 2
         CALL optx_lsd_calc(rho=rho(ispin)%array, norm_drho=ndrho(ispin)%array, &
                            e_0=e_0, e_rho=e_rho(ispin)%array, e_ndrho=e_ndrho(ispin)%array, &
                            npoints=npoints, epsilon_rho=epsilon_rho, &
                            epsilon_drho=epsilon_drho, sx=sx, &
                            a1=a1, a2=a2, gam=gam)
      END DO
   END SUBROUTINE optx_lsd_eval

! **************************************************************************************************
!> \brief optx exchange functional
!> \param rho the full density
!> \param norm_drho the norm of the gradient of the full density
!> \param e_0 the value of the functional in that point
!> \param e_rho the derivative of the functional wrt. rho
!> \param e_ndrho the derivative of the functional wrt. norm_drho
!> \param epsilon_rho the cutoff on rho
!> \param epsilon_drho ...
!> \param npoints ...
!> \param sx scaling-parameter for exchange
!> \param a1 a1 coefficient of the OPTX functional
!> \param a2 a2 coefficient of the OPTX functional
!> \param gam gamma coefficient of the OPTX functional
!> \par History
!>      01.2007 added scaling [Manuel Guidon]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE optx_lda_calc(rho, norm_drho, e_0, e_rho, e_ndrho, &
                            epsilon_rho, epsilon_drho, npoints, sx, &
                            a1, a2, gam)
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, norm_drho
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0, e_rho, e_ndrho
      REAL(kind=dp), INTENT(in)                          :: epsilon_rho, epsilon_drho
      INTEGER, INTENT(in)                                :: npoints
      REAL(kind=dp), INTENT(in)                          :: sx, a1, a2, gam

      REAL(KIND=dp), PARAMETER                           :: cx = 0.930525736349100_dp, &
                                                            o43 = 4.0_dp/3.0_dp

      INTEGER                                            :: ii
      REAL(KIND=dp)                                      :: denom, ex, gamxsxs, myndrho, myrho, &
                                                            rho43, tmp, xs

!$OMP     PARALLEL DO DEFAULT (NONE) &
!$OMP                 SHARED(rho, norm_drho, e_0, e_rho, e_ndrho) &
!$OMP                 SHARED(epsilon_rho, epsilon_drho, sx, npoints) &
!$OMP                 SHARED(a1, a2, gam) &
!$OMP                 PRIVATE(ii, myrho, myndrho, rho43, xs, gamxsxs) &
!$OMP                 PRIVATE(denom, ex, tmp)

      DO ii = 1, npoints
         ! we get the full density and need spin parts -> 0.5
         myrho = 0.5_dp*rho(ii)
         myndrho = 0.5_dp*MAX(norm_drho(ii), epsilon_drho)
         IF (myrho > 0.5_dp*epsilon_rho) THEN
            rho43 = myrho**o43
            xs = (myndrho/rho43)
            gamxsxs = gam*xs*xs
            denom = 1.0_dp/(1.0_dp + gamxsxs)
            ex = rho43*(a1*cx + a2*(gamxsxs*denom)**2)
            ! 2.0 for both spins
            e_0(ii) = e_0(ii) - (2.0_dp*ex)*sx
            tmp = rho43*2.0_dp*a2*gamxsxs*denom**2*(1.0_dp - gamxsxs*denom)
            ! derive e_0 wrt to rho (full) and ndrho (also full)
            e_rho(ii) = e_rho(ii) - ((o43*ex + tmp*gamxsxs*(-2.0_dp*o43))/myrho)*sx
            e_ndrho(ii) = e_ndrho(ii) - ((tmp*gam*2.0_dp*myndrho/rho43**2))*sx
         END IF
      END DO

!$OMP     END PARALLEL DO

   END SUBROUTINE optx_lda_calc

! **************************************************************************************************
!> \brief optx exchange functional
!> \param rho the *spin* density
!> \param norm_drho the norm of the gradient of the *spin* density
!> \param e_0 the value of the functional in that point
!> \param e_rho the derivative of the functional wrt. rho
!> \param e_ndrho the derivative of the functional wrt. norm_drho
!> \param epsilon_rho the cutoff on rho
!> \param epsilon_drho ...
!> \param npoints ...
!> \param sx scaling parameter for exchange
!> \param a1 a1 coefficient of the OPTX functional
!> \param a2 a2 coefficient of the OPTX functional
!> \param gam gamma coefficient of the OPTX functional
!> \par History
!>      01.2007 added scaling [Manuel Guidon]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE optx_lsd_calc(rho, norm_drho, e_0, e_rho, e_ndrho, &
                            epsilon_rho, epsilon_drho, npoints, sx, &
                            a1, a2, gam)
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, norm_drho
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0, e_rho, e_ndrho
      REAL(kind=dp), INTENT(in)                          :: epsilon_rho, epsilon_drho
      INTEGER, INTENT(in)                                :: npoints
      REAL(kind=dp), INTENT(in)                          :: sx, a1, a2, gam

      REAL(KIND=dp), PARAMETER                           :: cx = 0.930525736349100_dp, &
                                                            o43 = 4.0_dp/3.0_dp

      INTEGER                                            :: ii
      REAL(KIND=dp)                                      :: denom, ex, gamxsxs, myndrho, myrho, &
                                                            rho43, tmp, xs

!$OMP     PARALLEL DO DEFAULT(NONE) &
!$OMP                 SHARED(rho, norm_drho, e_0, e_rho, e_ndrho) &
!$OMP                 SHARED(epsilon_rho, epsilon_drho, npoints, sx) &
!$OMP                 SHARED(a1, a2, gam) &
!$OMP                 PRIVATE(ii, denom, ex, gamxsxs, myndrho, myrho) &
!$OMP                 PRIVATE(rho43, tmp, xs)

      DO ii = 1, npoints
         ! we do have the spin density already
         myrho = rho(ii)
         myndrho = MAX(norm_drho(ii), epsilon_drho)
         IF (myrho > epsilon_rho) THEN
            rho43 = myrho**o43
            xs = (myndrho/rho43)
            gamxsxs = gam*xs*xs
            denom = 1.0_dp/(1.0_dp + gamxsxs)
            ex = rho43*(a1*cx + a2*(gamxsxs*denom)**2)
            ! for a single spin
            e_0(ii) = e_0(ii) - ex*sx
            tmp = rho43*2.0_dp*a2*gamxsxs*denom**2*(1.0_dp - gamxsxs*denom)
            ! derive e_0 wrt to rho and ndrho
            e_rho(ii) = e_rho(ii) - ((o43*ex + tmp*gamxsxs*(-2.0_dp*o43))/myrho)*sx
            e_ndrho(ii) = e_ndrho(ii) - ((tmp*gam*2.0_dp*myndrho/rho43**2))*sx
         END IF
      END DO

!$OMP     END PARALLEL DO

   END SUBROUTINE optx_lsd_calc

END MODULE xc_optx
