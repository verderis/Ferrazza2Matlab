! *************************************************************************************
! MOD_PERMITTIVITY contains the algorithms that compute the dielectric constants
!
! VEGEDIEL
! SOILDIEL
! DOBSOIL
! *************************************************************************************
MODULE MOD_PERMITTIVITY

USE MOD_GLOBAL_PARAMETERS

PUBLIC	VEGEDIEL   ,&
				SOILDIEL   ,&
				DOBSOIL

CONTAINS

! SUBROUTINE VEGEDIEL
! CALLED BY: RAYGANS, PHO, CYL, CORNER
!
! Computes vegetation permittivity uSINg:
!   -   empirical fitting of experimental data in [10] for volumetric moisture <= 50%;
!   -   The model of [11] for volumetric moisture > 50%.
!
! Input  : gravimetric moisture (amoi)
! output : real (er) and imaginary (ej) parts of permittivity.
!
SUBROUTINE VEGEDIEL(dsdw,amoi,er,ej,select)

! Dummy variable
REAL, INTENT(IN)		  :: amoi,dsdw
INTEGER, INTENT(IN) 	:: select
REAL, INTENT(OUT)		 	:: er,ej

! Local variables
REAL	  :: eps0,fhz,t,sw,eoo,eo,atsw,eos,btsw,tau2PI,tau2PIs,dd,ffi,ss25,ss,erw,ejw,amd,zz,vw,er0
 
 eps0 = 8.854E-12 ! vacuum absolute permittivity
 fhz = F*1.E+09   ! frequency in Hz
 t = 20.          ! tyPIcal vegetation temperature (in C)
 sw = 10.         ! tyPIcal salinity (parts per thousand)
 !
 ! Liquid water permittivity by [12]
 !
 eoo = 4.9																                                          ! e16
 eo = 87.134 - .1949*t - .01276*t*t + .0002491*t*t*t							                  ! e23
 atsw = 1 + t*sw*1.613E-05 - sw*3.656E-03 + sw*sw*3.21E-05 - sw*sw*sw*4.232E-07     ! e24
 eos = eo*atsw															                                        ! e22
 btsw = 1 + t*sw*2.282E-05 - sw*7.638E-04 - sw*sw*7.76E-06 + sw*sw*sw*1.105E-08	    ! e26
 tau2PI = 1.1109E-10 - t*3.824E-12 + t*t*6.938E-14 - t*t*t*5.096E-16			          ! e17
 tau2PIs = tau2PI*btsw													                                    ! e25
 dd = 25 - t
 ffi = dd*(2.033E-02 + dd*1.266E-04 + dd*dd*2.464E-06 - sw*(1.849E-05 - dd*2.551E-07 + dd*dd*2.551E-08)) ! e28b
 ss25 = sw*(0.18252 - sw*1.4619E-03 + sw*sw*2.093E-05 - sw*sw*sw*1.282E-07)		      ! e28a
 ss = ss25*EXP( -ffi)														                                    ! e27
 erw = eoo + (eos - eoo)/(1 + (fhz*tau2PIs)**2)									                ! liquid water pemittivity (real)
 ejw = fhz*tau2PIs*(eos - eoo)/(1 + (fhz*tau2PIs)**2) + ss/(2*PI*fhz*eps0)		  ! liquid water permittivity (imaginary)
 !
 ! The model of [11] is used for higher moistures.
 !
 IF(select == 1) THEN
 ! Matzler
   amd = 1. - amoi
   zz = 0.522*(1 - 1.32*amd)
   er = zz*erw + 0.51 + 3.84*amd
   ej = zz*ejw
 !
 ! For lower moistures an emPIrical formula is used, obtained by
 ! fitting data published in [10] and impoSINg the liquid water
 ! permittivity to be obtained when the moisture is 100%.
 !
 ELSE
 ! Ulaby
   vw = amoi*dsdw/(1 - amoi*(1 - dsdw))  ! volumetric moisture
   er0 = 1.75
   er = er0 + (erw - er0)*(1 - EXP(-0.55*vw/(1. - vw)))
   ej = ejw*(1. - EXP(-0.4*vw/(1. - vw))) 
 END IF

END SUBROUTINE VEGEDIEL
! SUBROUTINE SOILDIEL
! CALLED BY: IEM, CORNER
!
! This SUBROUTINE computes soil permittivity uSINg
! the semiemPIrical model described in [12], pp. 2101 - 2104.
! Inputs: volumetric moisture, frequency.
! Outputs: real (egr) and imaginary (egj) parts of soil permittivity.
SUBROUTINE SOILDIEL(vmoi,egr,egj)

! Dummy variables
REAL, INTENT(IN)  :: vmoi
REAL, INTENT(OUT) :: egr,egj

! Local variables
COMPLEX           :: g2,g3,g4,egc
REAL			        :: lambda,alpha,beta,ROB,ross,epsss,g1,f0,fhz,g2b

 lambda = 30./F        ! wavelength, in cm
 alpha = 0.65          ! pag. 2103 of [12]
 beta = 1.1            ! pag. 2103 of [12] (tyPIcal value)
 ROB = 1.1             ! tyPIcal bulk density(g/cm^3)
 ross = 2.65           ! soil porosity (g/cm^3)
 epsss = 4.7           ! solid soil permittivity
 g1 = 1 + (ROB/ross)*(epsss**alpha - 1)
 f0 = F/18.64
 fhz = F*1.E+09

 IF(F <= 4) THEN
   g2 = 4.9 + 74.1/(CMPLX(1.,f0))
   g2b = .107/(2.*PI*fhz*8.854E-12)*(ross - ROB)/(ross*vmoi)
   g2 = (g2 - CMPLX(0.,g2b))**alpha
 ELSE
   g2 = (4.9 + 74.1/(CMPLX(1.,f0)))**alpha
 END IF

 g3 = g1 + (vmoi)**beta*(g2 - 1)
 g4 = (CLOG(g3))/alpha
 egc = CEXP(g4)
 egr = REAL(egc)
 egj = IMAG(egc)

END SUBROUTINE SOILDIEL
!
! SUBROUTINE SOILDIEL
! CALLED BY: IEM, CORNER
!
! This SUBROUTINE computes soil permittivity using
! the semiempirical model described in Dobson , 85, en IEEE, vol.23, no.1
! Inputs 
!					TSOL	: ground temperature in Kelvin
!					xmv	  : volumetric soil moisture
!					SAND	: SAND percentage 
!					CLAY	: caly percentage
!					ROB	  : dry matter density
!					epsol : dielectric constant
!
! Outputs 
!				  eg    : complex soil permittivity.
!
! Modello 
SUBROUTINE DOBSOIL(xmv,eg)
! Dummy argument declaration
	REAL, INTENT(IN)	  :: xmv
	COMPLEX,INTENT(OUT)	:: eg
!
! Local variables declaration
!
	COMPLEX   ::  epsol,epfw,cx
	REAL	    ::  permit,alp,ros,epwi,salsol,ts,epx,epy,epsi,epss,xa,xb,f0
	REAL	    ::  epw0,toPI,seff,bet1,bet2,eps,x,y,fr

	fr = F*1.e9
	! Air permitivitty
	permit = 8.854E-12
	! Alfa (in refractive model, ulaby appendix E.62)
	alp = 0.65	
	! Density of the solid soil material
	ros = 2.66 
	! e in the high freq limit (indep of salinity, Stogryn, 71, ulaby appendix E.112)
	epwi = 4.9
	salsol = 0.65
	ts = TSOL - 273.157
	IF (ts <= -0.5) THEN
	  ! SOIL   =   FROZEN SOIL	(Ulaby et al., 1986, p2101) 	
		epx = 5.
		epy =  - 0.5
		epsol = CMPLX(epx,epy)
	ELSEIF (xmv < .02 .AND. SAND >= .9) THEN
	  ! DRY SAND (ref Matzler, 1998)
		epsi = 2.53
		epss = 2.79
		xa = .002
		f0 = .27e9
		epsol = epsi + (epss - epsi)/CMPLX(1, - fr/f0) + CMPLX(0,xa)
	ELSE
		!	SALINE WATER
		! ref: saline water (Ulaby/M/F p2024)
		!
		! E.24 ulaby 2024
		xa = 1 + 1.613E-5*ts*salsol - 3.656E-3*salsol + 3.21E-5*salsol*salsol - 4.232E-7*(salsol**3)
		!	E.23 dependence of e saline water
		epw0 = xa*(87.134 - 1.949E-1*ts - 1.276E-2*ts*ts + 2.491E-4*(ts**3)) 
		! E.26: correction to the free water relaxation time due to salinity
		xb = 1 + 2.282E-5*ts*salsol - 7.638E-4*salsol - 7.76E-6*salsol*salsol + 1.105E-8*(salsol**3)
		! E.25, E.17(2*PI*relaxation time of free water)
		toPI = xb*(1.1109E-10 - 3.824E-12*ts + 6.938E-14*ts*ts - 5.096E-16*(ts**3))
		! Dobson E.32 corrected, effective conductivity, function of soil texture
		seff =  - 1.645 + 1.939*ROB - 2.256*SAND + 1.594*CLAY
		! Dobson E.29
		cx = (epw0 - epwi)/CMPLX(1.,toPI*fr)
		! modified Debye equation, Dobson E.29		 
		epfw = epwi + cx - CMPLX(0,seff*(ros - ROB))/(2*PI*fr*permit*ros*xmv)
		! Dobson E.30 (related to real(e))
		bet1 = 1.275 - 0.519*SAND - 0.152*CLAY
		! Dobson E.31 (related to imag(e))
		bet2 = 1.338 - 0.603*SAND - 0.166*CLAY
		! Dobson E.22 soil permit
		eps = (1.01 + 0.44*ros)**2 - 0.062
		! Ulaby E.111, Dobson E.28
		x = 1 + ROB*(eps**alp - 1)/ros + (xmv**bet1)*(REAL(epfw)**alp) - xmv
		y = (xmv**bet2)*(ABS(AIMAG(epfw))**alp)
		epx = x**(1/alp)
		epy = y**(1/alp)
		! frozen soil
		epsol = CMPLX(epx, -epy)
	END IF
  eg = epsol

END SUBROUTINE DOBSOIL

END MODULE MOD_PERMITTIVITY