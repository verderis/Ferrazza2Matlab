! *************************************************************************************
! MOD_VEGETATION_FUNCTIONS contains the algorithms that compute the matrix scattering
! and extinction of dielectric bodies.
!
! RAYLEIGH_GANS 
! PHYSICAL_OPTICS
! CURVED_SHEET
! INFINITE_LENGTH
! CORNER
! *************************************************************************************
MODULE MOD_VEGETATION_FUNCTIONS

USE MOD_GLOBAL_PARAMETERS
USE MOD_PERMITTIVITY

IMPLICIT NONE

! Public declaration of member functions defined inside the module
PUBLIC	RAYLEIGH_GANS				,&
				PHYSICAL_OPTICS			,&
				CURVED_SHEET				,&
				INFINITE_LENGTH			,&
				HOLLOW_CYLINDER			,&
				TRUNKS_ABS					,&
				CORNER							
CONTAINS
!
! THIS SUBROUTINE CONTAINS OTHER SUBROUTINES  -  Osborn, Disc, Needle, Integr and FFATT
! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE RAYLEIGH_GANS
! CALLED BY: MAIN
! Computes the scatter matrix and the extinction vector uSINg the
! Rayleigh - Gans approximation for:
! 1) An ensemble of circular discs (a > l);
! 2) An ensemble of needles (l > a).
!
SUBROUTINE RAYLEIGH_GANS(a,b,l,amoi,dsdw,s,ke)
!Dummy varaiables declaration
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
REAL, DIMENSION(:,:), INTENT(OUT)			    :: ke
REAL, INTENT(IN)							            :: a,b,l,amoi,dsdw
! Local variabled declaration
COMPLEX, DIMENSION(3,3)		:: ua,u2,uaa,u1,ucom
COMPLEX, DIMENSION(3,1)		:: cv,ch,av,ah
COMPLEX, DIMENSION(NM)		:: xx

COMPLEX						        :: cscatc,ecom,fvv,fvh,fhv,fhh

REAL, DIMENSION(NM,2,2)		:: xff
REAL, DIMENSION(NN/2,2)		:: ks,ka
REAL, DIMENSION(3,3)		  :: u
REAL, DIMENSION(3,1)		  :: v,q,q2,h
REAL, DIMENSION(2,2)		  :: xdirez
REAL, DIMENSION(3)			  :: hs,vs
REAL, DIMENSION(2)			  :: jcos,icos,kstoro

REAL						:: lambda,k,kk,j1aqt,lqt,kav,kah,kavs,kahs,kdirez1,kdirez2,cscatr,cscatj,ccabs,thetaj,thetais,theta
REAL						:: stheta,ctheta,absctheta,thetas,sthetas,cthetas,abscthetas,o0,cdirez,phis,sphis,cphis
REAL						:: alpha1,beta1,gamma1,dalpha,dbeta,dgamma,delthet,delphis,c,er,ej,alpha,salpha,calpha,beta,sbeta,cbeta
REAL						:: gamma,sgamma,cgamma,qt,aqt,vscat,fvvfvv,fhhfhh,fhvfhv,fvhfvh,rl,rm,rn,BESJ1

INTEGER					:: nalpha,nbeta,ngamma,nabg,i,j,icode,k1,k2,k3,k3a,ialpha,ibeta,igamma,istoki,istoks

DATA ua/1.,0.,0.,0.,1.,0.,0.,0.,1./

DATA jCOS/ -1, -1/
DATA iCOS/1, -1/


IF(a > l) THEN
  ! Discs selection 
  alpha1 = ALPHA1DIS
  nalpha = NALPHADIS
  dalpha = DALPHADIS
  beta1 = BETA1DIS
  nbeta = NBETADIS
  dbeta = DBETADIS
  gamma1 = GAMMA1DIS
  ngamma = NGAMMADIS
  dgamma = DGAMMADIS
ELSE
  ! Petioles or ears selection
  alpha1 = ALPHA1PET
  nalpha = NALPHAPET
  dalpha = DALPHAPET
  beta1 = BETA1PET
  nbeta = NBETAPET
  dbeta = DBETAPET
  gamma1 = GAMMA1PET
  ngamma = NGAMMAPET
  dgamma = DGAMMAPET
ENDIF

nabg = nalpha*nbeta*ngamma
delthet = PI/(2*NIJ)	! interval amplitude in theta and theta_s
delphis = 2.*PI/nm		! interval amplitude in phi_s - phi

lambda = 30./F				! wavelength, in cm
k = 2.*PI/lambda			! wavenumber, in cm^ - 1
kk = k*k
c = 3.*l/4.					  ! ellipsoid small semi - axis (discs), in cm ellipsoid large semi - axis (needles), in cm

CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)  
                                    ! Permittivity computation
ecom = CMPLX(er,ej)

CALL OSBORN(a,b,c,rl,rm,rn)         ! Demagnetization factors 2.23 of [3]
ua(1,1) = 1./(1. + rl*(ecom - 1))   ! Diagonal elements in (24) of [4] and (2) of [6]
ua(2,2) = 1./(1. + rm*(ecom - 1))      
ua(3,3) = 1./(1. + rn*(ecom - 1))
!
cscatr = kk*a*b*l*(er - 1)          ! factors in (28) of [4], (5) of [6]
cscatj = kk*a*b*l*ej
cscatc = CMPLX(cscatr,cscatj)

ccabs = k*PI*a*b*l*ej			          ! fattore coeff. assorbimento Ref.1 appendice 

jLoop: DO j = 1,NIJ                 ! j = incidence angular interval index

  thetaj = (j - 1)*delthet          ! lower limit of angular interval in theta

  ks(j,:) = 0. ; ka(j,:) = 0.

  icodeLoop: DO icode = 1,2
    ! icode = 1  : upper half - space scattering
    ! icode = 2  : lower half - space scattering

    iLoop: DO i = 1,NIJ             ! i = scattering angular interval index
      
      thetais = (i - 1)*delthet     ! lower limit of angular interval in theta_s
      !
      kstoro(:) = 0.
      xff(:,:,:) = 0.

      k1Loop: DO k1 = 1,NS              ! integration within jth theta interval using Gauss - Legendre technique

        theta = thetaj + .5*delthet*(1. + chi(k1))
        stheta = SIN(theta)
        absctheta = COS(theta)
        ctheta = jCOS(icode)*absctheta

        k2Loop: DO k2 = 1,NS		        ! integration within ith theta_s interval using Gauss - Legendre technique

          thetas = thetais + .5*delthet*(1. + chi(k2))
          sthetas = SIN(thetas)
          abscthetas = COS(thetas)
          cthetas = iCOS(icode)*abscthetas
          o0 = om(k1)*om(k2)/4.         ! to be used in the double integration in theta and theta_s
          !
          ! The "s" functions returned to the main program include the
          ! (delthet*stheta/abscthetas) factor, to be used in formulas (3)
          ! and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0

          ! Computations for 0 < phi_s - phi < 180
          k3Loop: DO k3 = 1,NM1

            phis = (k3 - 1)*delphis   ! phi_s - phi
            sphis = SIN(phis)
            cphis = COS(phis)
            !
            !  Computations at page 193 of [4], 141 of [6]. Incidence and scattering polarization vectors.
            v(1,1) = ctheta           ! V pol, incidence
            v(2,1) = 0
            v(3,1) =  - stheta
            vs(1) = cthetas*cphis     ! V pol, scattering
            vs(2) = cthetas*sphis
            vs(3) =  - sthetas
            h(1,1) = 0                ! H pol, incidence
            h(2,1) = 1.
            h(3,1) = 0
            hs(1) =  - sphis          ! H pol, scattering
            hs(2) = cphis
            hs(3) = 0
            q(1,1) = sthetas*cphis - stheta       ! wavenumber difference vector (page 192 of [4])
            q(2,1) = sthetas*sphis
            q(3,1) = cthetas - ctheta
            !
            xdirez(:,:) = 0.
            IF(k2 == 1 .AND. k3 == 1 .AND. icode == 1 .AND. i == 1) THEN
              kavs = 0.
              kahs = 0.
            END IF
            !                  (averaging over scatterer orientation)
            alphaLoop: DO ialpha = 1,nalpha
              alpha = alpha1 + (ialpha - 1)*dalpha
              salpha = SIND(alpha)
              calpha = COSD(alpha)

              betaLoop: DO ibeta = 1,nbeta
                beta = beta1 + (ibeta - 1)*dbeta
                sbeta = SIND(beta)
                cbeta = COSD(beta)

                gammaLoop: DO igamma = 1,ngamma
                  gamma = gamma1 + (igamma - 1)*dgamma
                  sgamma = SIND(gamma)
                  cgamma = COSD(gamma)
                  ! Reference system transformation matrix, (22) of [4]
                  u(1,1) = cbeta*calpha
                  u(1,2) = salpha*cbeta
                  u(1,3) = sbeta
                  u(2,1) =  - calpha*sbeta*sgamma - salpha*cgamma
                  u(2,2) =  - salpha*sbeta*sgamma + calpha*cgamma
                  u(2,3) = cbeta*sgamma
                  u(3,1) =  - calpha*sbeta*cgamma + salpha*sgamma
                  u(3,2) =  - salpha*sbeta*cgamma - calpha*sgamma
                  u(3,3) = cbeta*cgamma
                  ! (24) of [4]
                  ucom(:,:) = u(:,:)
                  u1(:,:) = TRANSPOSE(u(:,:))            ! u^ - 1
                  u2(:,:) = MATMUL(u1(:,:),ua(:,:))      ! u^ - 1[ - ]
                  uaa(:,:) = MATMUL(u2(:,:),ucom(:,:))   ! u^ - 1[ - ]u
                  q2(:,:) = MATMUL(u(:,:),q(:,:))

                  IF(a > c)THEN
                    ! Disc CASE (form factor in (28) of [4])
                    qt = k*SQRT((q2(1,1)*q2(1,1)*a*a) + (q2(2,1)*q2(2,1)*b*b))
                    aqt = qt
                    IF(qt == 0.) THEN
                      vscat = .5/2.
                    ELSE					  
                      j1aqt = BESJ1(aqt)					  
                      vscat = j1aqt/(2.*qt)
                    END IF
                    !
                    ! Needle CASE (form factor in (5) of [6])
                    ELSE
                      lqt = k*l*q2(3,1)/2.
                      IF(q2(3,1) == 0.)THEN
                        vscat = 1./4.
                      ELSE
                        vscat = SIN(lqt)/(4.*lqt)
                    END IF
                  END IF
                  !
                  cv(:,:) = v(:,:)
                  ch(:,:) = h(:,:)
                  ! Transformation to principal reference frame ((28) of [4], (5) of [6])
                  av(:,:) = MATMUL(uaa(:,:),cv(:,:)) 
                  ah(:,:) = MATMUL(uaa(:,:),ch(:,:))
                  ! Elements of scattering amplitude matrix
                  fvv = cscatc*vscat*(av(1,1)*vs(1) + av(2,1)*vs(2) + av(3,1)*vs(3))
                  fhv = cscatc*vscat*(av(1,1)*hs(1) + av(2,1)*hs(2) + av(3,1)*hs(3))
                  fhh = cscatc*vscat*(ah(1,1)*hs(1) + ah(2,1)*hs(2) + ah(3,1)*hs(3))
                  fvh = cscatc*vscat*(ah(1,1)*vs(1) + ah(2,1)*vs(2) + ah(3,1)*vs(3))

                  IF(k2 == 1 .AND. k3 == 1 .AND. icode ==  1 .AND. i == 1) THEN
                    KAV = AV(1,1)*CONJG(AV(1,1)) + AV(2,1)*CONJG(AV(2,1)) + AV(3,1)*CONJG(AV(3,1))
                    KAH = AH(1,1)*CONJG(AH(1,1)) + AH(2,1)*CONJG(AH(2,1)) + AH(3,1)*CONJG(AH(3,1))
                    kavs = kavs + KAV
                    kahs = kahs + KAH
                  END IF		

                  fvvfvv = fvv*CONJG(fvv)
                  fvhfvh = fvh*CONJG(fvh)
                  fhvfhv = fhv*CONJG(fhv)
                  fhhfhh = fhh*CONJG(fhh)

                  xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nabg
                  xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nabg
                  xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nabg
                  xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nabg

                  xdirez(1,1) = xdirez(1,1) + fvvfvv
                  xdirez(1,2) = xdirez(1,2) + fvhfvh
                  xdirez(2,1) = xdirez(2,1) + fhvfhv
                  xdirez(2,2) = xdirez(2,2) + fhhfhh

                END DO gammaLoop

              END DO betaLoop

            END DO alphaLoop

            kdirez1 = O0*sthetas*(xdirez(1,1) + xdirez(2,1))
            kdirez2 = O0*sthetas*(xdirez(1,2) + xdirez(2,2))
            !     vedi appendice
            IF(K3 == 1 .OR. K3 == NM1) THEN
            kdirez1 = kdirez1/2.
            kdirez2 = kdirez2/2.
            END IF
            kstoro(1) = kstoro(1) + kdirez1/(nabg*absctheta)
            kstoro(2) = kstoro(2) + kdirez2/(nabg*absctheta)

          END DO k3Loop

        END DO k2Loop

        IF(icode == 1 .AND. i == 1) THEN
          ka(j,1) = ka(j,1) + kavs*ccabs*om(k1)/(2*nabg*absctheta)   
          ka(j,2) = ka(j,2) + kahs*ccabs*om(k1)/(2*nabg*absctheta)
        END IF

      END DO K1Loop
      kstoro(:) = kstoro(:)*2.
      ks(j,:) = ks(j,:) + kstoro(:)*delthet*delphis  ! somma contributi dei tori di scattering 
      !
      !  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      DO k3 = nm/2 + 2,nm
        k3a = 2*NM1 - k3
        xff(k3,:,:) = xff(k3a,:,:)
      END DO  ! k3
      !
      !      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
      !	
      DO istoki = 1,2
        DO istoks = 1,2
          xx(:) = CMPLX(xff(:,istoks,istoki), 0.)
          CALL CFSTFT(NF,xx)
          ! normalization
          xx(1) = xx(1)/NM
          xx(2:) = xx(2:)*2/NM
          s(i,j,icode,1:NM1,istoks,istoki) = REAL(xx(1:NM1))	
        END DO   
      END DO    

    END DO iLoop

  END DO icodeLoop

END DO jLoop

! Extinction cross senction
ke(:,:) = ks(:,:) + ka(:,:)

CONTAINS

	!******************************************************************
	SUBROUTINE OSBORN(a,b,c,rl,rm,rn)
	! CALLed by RAYLEIGH_GANS
	! 2.23 - 2.25 of [3] (demagnetization factors)

	! Dumyy variables declaration
	REAL, INTENT(IN)    :: a,b,c
	REAL, INTENT(OUT)   :: rl,rm,rn

	! Local variables declaration
	REAL				        :: rapba,rapca,e,e2,e1,ellk,elle
  
	IF((a == b) .AND. (a > c)) THEN
		CALL DISC(a,c,rl,rm,rn)
		RETURN
	ELSEIF((a == b) .AND. (c > a)) THEN
		CALL NEEDLE(a,c,rl,rm,rn)
		RETURN
	ELSE
		rapba = b/a
		rapca = c/a
		e = SQRT(1. - (rapba*rapba))
		e2 = e*e
		e1 = 1. - e2
		CALL INTEGR(e2,1,ellk)
		CALL INTEGR(e2,-1,elle)
		rl = rapca*SQRT(e1)*(ellk - elle)/(e2)
		rm = rapca*(elle-e1*ellk)/(e2*SQRT(e1))
		rn = 1. - rapca*elle/SQRT(e1)
	END IF

	END SUBROUTINE OSBORN
	!******************************************************************
	SUBROUTINE DISC(a,c,rl,rm,rn)
	! CALLed by OSBORN

	! Dummy variables declaration
	REAL, INTENT(IN)  :: a,c
	REAL, INTENT(OUT) :: rm,rl,rn

	! Local variables declaration
	REAL			        :: aa,cc,diff,rad,atg,ac3

	aa = a*a
	cc = c*c
	diff = aa - cc
	rad = SQRT(diff)
	atg = ATAN(SQRT(cc/diff))
	ac3 = (diff)**1.5
	rl = ( -c*rad/aa + PI/2. - atg)/ac3
	rl = aa*c/2.*rl
	rm = rl
	rn = 2.*(rad/c - PI/2. + atg)/ac3
	rn = aa*c/2.*rn

	END SUBROUTINE DISC
	!******************************************************************
	SUBROUTINE NEEDLE(a,c,rl,rm,rn)
	! CALLed by OSBORN

	! Dummy variables declaration
	REAL, INTENT(IN)  ::a,c
	REAL, INTENT(OUT) :: rl,rm,rn

	! Local variables declaration
	REAL			        :: aa,cc,ac1,ac2,ac3,aa0,aa1

	aa = a*a
	cc = c*c
	ac1 = SQRT(cc - aa)
	ac2 = ALOG((c - ac1)/(c + ac1))
	ac3 = 1/(cc - aa)
	aa0 = ac3*(c/aa + (0.5/ac1)*ac2)
	aa1 = -ac3*((1/ac1)*ac2 + 2./a)
	rl = aa*c/2.*aa0
	rm = rl
	rn = aa*c/2.*aa1

	END SUBROUTINE NEEDLE
	!******************************************************************
	SUBROUTINE INTEGR(rm,iflag,sum)
	! CALLed by OSBORN
	!
	!  Elliptic integrals computation [9]
	!  iflag = 1 : first kind
	!  iflag =  - 1: second kind

	! Dummy variables declaration
	REAL, INTENT(IN)    :: rm
	REAL, INTENT(OUT)   :: sum
	INTEGER, INTENT(IN) :: iflag

	! Local variables declaration
	REAL				  :: denom,coeff,param
	INTEGER				:: n,i,par,disp

	sum = PI/2.
	denom = 1.
	n = 10
	DO i = 1,n
		disp = 2*i - 1
		par = 2*i
		IF(iflag == - 1) denom = disp
		coeff = (ffatt(disp)/ffatt(par))*(ffatt(disp)/ffatt(par))
		param = rm**i
		sum = sum + PI/2.*(iflag*coeff*param/denom)    ! eq. 17.3.11 and 17.3.12 of [9]
	END DO

	END SUBROUTINE INTEGR

	!******************************************************************
	REAL FUNCTION FFATT(rk)
	! called by integr
	! double factorial function
	! Dummy variables declaration
	INTEGER, INTENT(IN)   :: rk
	! Local variables declaration
	INTEGER			          :: i,nk

	ffatt = 1.
	nk = rk
	DO i = nk,1, -2
	ffatt = ffatt*i
	END DO

	END FUNCTION FFATT

END SUBROUTINE RAYLEIGH_GANS
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE PHYSICAL_OPTICS
! CALLED BY: MAIN
! Computes the scatter matrix and the extinction vector using the
! Physical - Optics approximation for an ensemble of circular discs.
!
SUBROUTINE PHYSICAL_OPTICS(a,b,l,amoi,dsdw,s,ke)
! Dummy variables declaration
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
REAL, DIMENSION(:,:), INTENT(OUT)			    :: ke
REAL, INTENT(IN)							            :: a,b,l,amoi,dsdw
! Local declaration variables
COMPLEX, DIMENSION(NM)	:: xx
COMPLEX, DIMENSION(3,3)	:: v2p,v2m,h2,h2ehp,v2evp,h2ehm,v2evm,alev,epiu,emen,epsin,emsin
COMPLEX, DIMENSION(3,1) :: uvepsp,uvepsm,fv,fh
COMPLEX, DIMENSION(2,2) :: eq
COMPLEX, DIMENSION(2)   :: r,t,sig,sinc,tet

COMPLEX					        :: ubetz,betmen,betpiu,cscatc,alfexp,psi,denom,ecom,radeps,fvv,fhv,fvh,fhh

REAL, DIMENSION(NM,2,2) :: xff
REAL, DIMENSION(3,1)    :: h,v,av,ah
REAL, DIMENSION(1,3)	  :: ht,vt
REAL, DIMENSION(3)		  :: vs,hs,q,o,ri,rn,hvn,rx,ry
REAL, DIMENSION(2)		  :: icos,jcos,q2

REAL					:: delthet,delphis,er,ej,cscatr,cscatj,BESJ1
REAL					:: lambda,k,kk,j1qt,gamma,beta,alpha,theta,thetas,thetaj,thetais,stheta,ctheta,absctheta,sthetas,abscthetas
REAL					:: cthetas,o0,cdirez,phis,sphis,cphis,salpha,calpha,sbeta,cbeta,sgamma,cgamma,cthetau,tvi,thi,tt,sthetau
REAL					:: on,qt,aqt,vscat,cke,fvvfvv,fvhfvh,fhvfhv,fhhfhh,alpha1,beta1,gamma1,dalpha,dbeta,dgamma

INTEGER				:: ialpha,ibeta,igamma,nalpha,nbeta,ngamma,i,j,k1,k2,k3,icode,istoki,istoks,kj,nabg,k3a

DATA jcos/ -1, -1/
DATA icos/1, -1/

IF(a > l) THEN
  ! Discs selection 
  alpha1 = ALPHA1DIS
  nalpha = NALPHADIS
  dalpha = DALPHADIS
  beta1 = BETA1DIS
  nbeta = NBETADIS
  dbeta = DBETADIS
  gamma1 = GAMMA1DIS
  ngamma = NGAMMADIS
  dgamma = DGAMMADIS
ELSE
  ! Petioles or ears selection
  alpha1 = ALPHA1PET
  nalpha = NALPHAPET
  dalpha = DALPHAPET
  beta1 = BETA1PET
  nbeta = NBETAPET
  dbeta = DBETAPET
  gamma1 = GAMMA1PET
  ngamma = NGAMMAPET
  dgamma = DGAMMAPET
ENDIF

nabg = nalpha*nbeta*ngamma

delthet = PI/(2*NIJ)     ! interval amplitude in theta and theta_s
delphis = 2.*PI/nm       ! interval amplitude in phi_s - phi
!
lambda = 30./F           ! wavelength, in cm
k = 2.*PI/lambda         ! wavenumber, in cm^ - 1
kk = k*k
! 
CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)    
                               ! permittivity computation
ecom = CMPLX(er,ej)
!
cscatr = kk*a*b*l*(er - 1)     ! (22) of [5]
cscatj = kk*a*b*l*ej
cscatc = CMPLX(cscatr,cscatj)
!
sig(1) = 1./ecom               ! (20c) of [5]
sig(2) = (1.0,0.0)
radeps = CSQRT(ecom)

jLoop: DO j = 1,NIJ            ! j = incidence angular interval index
  thetaj = (j - 1)*delthet     ! lower limit of angular interval in theta
  ke(j,:) = 0.                 ! initialization to compute extinction coefficients
  !
  icodeLoop: DO icode = 1,2
    ! icode = 1  : upper half - space scattering
    ! icode = 2  : lower half - space scattering
    !

    iLoop: DO i = 1,NIJ          ! i = scattering angular interval index
      thetais = (i - 1)*delthet  ! lower limit of angular interval in theta_s
      xff(1:NM1,:,:) = 0.

      k1Loop: DO k1 = 1,2        ! integration within jth theta interval using Gauss - Legendre technique
        theta = thetaj + .5*delthet*(1. + chi(k1))
        stheta = SIN(theta)
        absctheta = COS(theta)
        ctheta = jCOS(icode)*absctheta
        ! Incidence direction vector
        ri(1) = stheta
        ri(2) = 0.
        ri(3) = ctheta
        ! Incidence vertical polarization vector
        v(1,1)  = ctheta
        v(2,1) = 0.0
        v(3,1) =  -stheta
        ! Incidence horizontal polarization vector
        h(1,1) = 0.0
        h(2,1)  = 1.0
        h(3,1)  = 0.0

        k2Loop: DO k2 = 1,2     ! integration within ith theta_s interval using Gauss - Legendre technique
          thetas = thetais + .5*delthet*(1. + chi(k2))
          sthetas = SIN(thetas)
          abscthetas = COS(thetas)
          cthetas = iCOS(icode)*abscthetas
          o0 = om(k1)*om(k2)/4.  ! to be used in the double integration in theta and theta_s
          ! The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor, 
          ! to be used in formulas (3) and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0

          ! Computations for 0 < phi_s - phi < 180
          k3Loop: DO k3 = 1,NM1
            phis = (k3 - 1)*delphis    ! phi_s - phi
            sphis = SIN(phis)
            cphis = COS(phis)
            !
            ! Scattering direction vector
            o(1) = sthetas*cphis
            o(2) = sthetas*sphis
            o(3) = cthetas
            ! Scattering vertical polarization vector
            vs(1) = cthetas*cphis
            vs(2) = cthetas*sphis
            vs(3) =  -sthetas
            ! Scattering horizontal polarization vector
            hs(1) =  -sphis
            hs(2) = cphis
            hs(3) = 0.

            alphaLoop: DO ialpha = 1,nalpha
              alpha = alpha1+(ialpha - 1)*dalpha
              salpha = SIND(alpha)
              calpha = COSD(alpha)

                betaLoop: DO ibeta = 1,nbeta
                  beta = beta1+(ibeta - 1)*dbeta
                  sbeta = SIND(beta)
                  cbeta = COSD(beta)

                  gammaLoop: DO igamma = 1,ngamma
                    gamma = gamma1+(igamma - 1)*dgamma
                    sgamma = SIND(gamma)
                    cgamma = COSD(gamma)
                    ! Unit vector perpendicular to disc
                    rn(1) = calpha*sbeta*cgamma + salpha*sgamma
                    rn(2) = salpha*sbeta*cgamma - calpha*sgamma
                    rn(3) = cbeta*cgamma
                    ! Direction of x axis in local reference frame
                    rx(1) = cbeta*calpha
                    rx(2) = cbeta*salpha
                    rx(3) =  - sbeta
                    ! Direction of y axis in local reference frame
                    ry(1) = calpha*sbeta*sgamma - salpha*cgamma
                    ry(2) = calpha*cgamma + salpha*sbeta*sgamma
                    ry(3) = cbeta*sgamma
                    ! Cosine of incidence direction in local reference frame
                    cthetau = ri(1)*rn(1) + ri(3)*rn(3)
                    IF(cthetau > 0.) THEN
                      cthetau =  - cthetau
                      rn(:) =  - rn(:)
                      rx(:) =  - rx(:)
                      ry(:) =  - ry(:)
                    END IF
                    sthetau = SQRT(1. - cthetau*cthetau)
                    ubetz = CSQRT(ecom - sthetau*sthetau) ! page 1257 of [5]
                    alfexp = CMPLX(0.,k*l/2.)
                    psi = 4*alfexp*ubetz
                    betmen = cthetau - ubetz
                    betpiu = cthetau + ubetz

                    tvi = (salpha*sgamma + calpha*sbeta*cgamma)*ctheta - cbeta*cgamma*stheta
                    thi = salpha*sbeta*cgamma - calpha*sgamma    
                    tt = SQRT(tvi*tvi + thi*thi)
                    ! Vertical and horizontal polarization vectors in local reference frame
                    IF(beta == 0. .AND. gamma == 0.)THEN
                      av = v(:,:)
                      ah = h(:,:)
                    ELSE
                      av(:,:) =  - (tvi*v(:,:) + thi*h(:,:))/tt
                      ah(:,:) = (thi*v(:,:) - tvi*h(:,:))/tt
                    END IF
                    ht(:,:) = TRANSPOSE(ah(:,:))
                    vt(:,:) = TRANSPOSE(av(:,:))
                    CALL PRODV(ah,rn,hvn)									        ! 18a of [5]
                    uvepsp(:,1) = (ubetz - cthetau)*hvn(:)
                    uvepsp(:,1) = (av(:,1) + uvepsp(:,1))/radeps	! 18b of [5]
                    uvepsm(:,1) =  - (ubetz + cthetau)*hvn(:)
                    uvepsm(:,1) = (av(:,1) + uvepsm(:,1))/radeps	! 18b of [5]

                    DO kj = 1,2
                      denom = cthetau - sig(kj)*ubetz
                      r(kj) = (cthetau + sig(kj)*ubetz)/denom			! 20a of [5]
                      t(kj) = 2.*CSQRT(sig(kj))*cthetau/denom			! 20b of [5]
                    END DO
                    DO kj = 1,2
                      denom = 1. - r(kj)*r(kj)*CEXP(psi)
                      eq(kj,1) =  - t(kj)*r(kj)*CEXP(psi)*CEXP(alfexp*betmen)/denom   ! 19a of [5]
                      eq(kj,2) = t(kj)*CEXP(alfexp*betpiu)/denom					            ! 19b of [5]
                    END DO
                    h2(:,:) = MATMUL(ah(:,:),ht(:,:))
                    v2p(:,:) = MATMUL(uvepsp(:,:),vt(:,:))
                    v2m(:,:) = MATMUL(uvepsm(:,:),vt(:,:))
                    h2ehp(:,:) = h2(:,:)*eq(2,1)
                    v2evp(:,:) = v2p(:,:)*eq(1,1)
                    h2ehm(:,:) = h2(:,:)*eq(2,2)
                    v2evm(:,:) = v2m(:,:)*eq(1,2)
                    epiu(:,:) = h2ehp(:,:) + v2evp(:,:)								    ! 23a of [5]
                    emen(:,:) = h2ehm(:,:) + v2evm(:,:)
                    on = DOT_PRODUCT(o(:),rn(:))
                    tet(1) = k*l*(on + ubetz)/2.						! 23c of [5]
                    tet(2) =  - k*l*(on - ubetz)/2.
                    sinc(:) = CSIN(tet(:))/tet(:)
                    epsin(:,:) = epiu(:,:)*sinc(2)
                    emsin(:,:) = emen(:,:)*sinc(1)
                    ! Scattering amplitude matrix in local reference frame
                    alev(:,:) = epsin(:,:) + emsin(:,:)
                    ! Wavenumber difference vector
                    q(:) = o(:) - ri(:)
                    ! Transformation to local reference frame
                    q2(1) = DOT_PRODUCT(q(:),rx(:))
                    q2(2) = DOT_PRODUCT(q(:),ry(:))
                    qt = k*SQRT(a*a*q2(1)*q2(1) + b*b*q2(2)*q2(2))
                    ! Form factor of circular disc ((23d) + (23e) + (23f) of [5])
                    aqt = qt
                    IF(qt == 0.) THEN
                      vscat = .5/2.
                    ELSE
                      j1qt = BESJ1(aqt)					 
                      vscat = j1qt / (2. * aqt)
                    END IF
                    ! Transformation to principal reference frame  ((28)of [4])
                    fv(:,:) = MATMUL(alev(:,:),v(:,:))
                    fh(:,:) = MATMUL(alev(:,:),h(:,:))
                    ! Elements of the scattering amplitude matrix
                    fvv = vs(1)*fv(1,1) + vs(2)*fv(2,1) + vs(3)*fv(3,1)
                    fhv = hs(1)*fv(1,1) + hs(2)*fv(2,1) + hs(3)*fv(3,1)
                    fvh = vs(1)*fh(1,1) + vs(2)*fh(2,1) + vs(3)*fh(3,1)
                    fhh = hs(1)*fh(1,1) + hs(2)*fh(2,1) + hs(3)*fh(3,1)
                    fvv = fvv*vscat*cscatc
                    fhv = fhv*vscat*cscatc
                    fvh = fvh*vscat*cscatc
                    fhh = fhh*vscat*cscatc
                    ! Extinction cross - section computation using the forward scattering theorem
                    IF(icode == 2 .AND. i == j .AND. k2 == k1 .AND. k3 == 1) THEN
                      ! ke values RETURNed to the main program are divided by absctheta as requested in formula (8) of [1]
                      cke = om(k1)/(2*absctheta*nabg)
                      ke(j,1) = ke(j,1) + cke*ABS(4*PI/k*AIMAG(fvv))
                      ke(j,2) = ke(j,2) + cke*ABS(4*PI/k*AIMAG(fhh))
                    END IF
                    !
                    fvvfvv = fvv*CONJG(fvv)
                    fvhfvh = fvh*CONJG(fvh)
                    fhvfhv = fhv*CONJG(fhv)
                    fhhfhh = fhh*CONJG(fhh)
                    !
                    xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nabg
                    xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nabg
                    xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nabg
                    xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nabg

                  END DO gammaLoop
                END DO betaLoop
            END DO alphaLoop

          END DO k3Loop

        END DO k2Loop

      END DO k1Loop
      !
      !  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      !
      DO k3 = nm/2 + 2,nm
        k3a = 2*NM1 - k3
        xff(k3,:,:) = xff(k3a,:,:)
      END DO
      !
      !      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
      DO istoki = 1,2
        DO istoks = 1,2
          xx(:) = CMPLX(xff(:,istoks,istoki), 0.)
          CALL CFSTFT(NF,xx)
          ! normalization
          xx(1) = xx(1)/NM
          xx(2:) = xx(2:)*2/NM
          s(i,j,icode,1:NM1,istoks,istoki) = REAL(xx(1:NM1))	      
        END DO
      END DO    

    END DO iLoop

  END DO icodeLoop

END DO jLoop

RETURN
CONTAINS
!*********************************************************
SUBROUTINE PRODV(a,b,c)
! called by PHYSICAL_OPTICS
! vector product;
REAL, INTENT(IN)  :: a(3,1),b(3)
REAL, INTENT(OUT) :: c(3)

c(1) = a(2,1)*b(3) - a(3,1)*b(2)
c(2) = a(3,1)*b(1) - a(1,1)*b(3)
c(3) = a(1,1)*b(2) - a(2,1)*b(1)

RETURN 
END SUBROUTINE PRODV

END SUBROUTINE PHYSICAL_OPTICS
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE CurvedSheet
! CALLED BY: MAIN
! Computes the scatter matrix and the extinction vector using the
! curved sheet approximation for a whole leaf.
!
SUBROUTINE CURVED_SHEET(ac,bc,lc,amoi,dsdw,s,ke)
! Dummy variables declaration
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
REAL, DIMENSION(:,:),INTENT(OUT)			    :: ke
REAL,INTENT(IN)														:: ac,bc,lc,amoi,dsdw
! Local variables declaration
COMPLEX, DIMENSION(3,2)		:: SS,Sl,Sm
COMPLEX, DIMENSION(NM)		:: xx
COMPLEX										:: fvv,fvh,fhv,fhh,molt,RR

REAL, DIMENSION(NM,2,2)		:: xff
REAL, DIMENSION(3,2)			:: es

REAL		:: delta_theta,delta_phi,thetamin,phi_l,theta,thetas,phis,theta_f,phi_f,stheta
REAL		:: thetaj,thetais,delthet,delphis,abscthetas,cdirez,ctheta,sthetas,cthetas,sphis,cphis
REAL		:: alpharad,den_r,den_i,thetamax,o0,cke,fvvfvv,fhhfhh,fvhfvh,fhvfhv,r,k,z,soglia,er,ei,lambda,a,b,l

INTEGER	:: nphi,dphi,l2,m2_1,m2,ck,test,i,j,k1,k2,k3,istoki,istoks,k3a,icode,i_phi,i_elev,i_azim

!***************************************************************************************
! ANGULAR CURVATURE OF DIELECTRIC SHEET
!***************************************************************************************
INTEGER, PARAMETER  :: ANGLE_CURVATURE = 90			! Curvature angle of dielectric curved sheet


fvv = (0.,0.)
fvh = (0.,0.)
fhv = (0.,0.)
fhh = (0.,0.)

! conversion from centimeter to meter of sheet dimension
a = ac/100.                
b = bc/100.
l = lc/100.

delthet = PI/(2*NIJ)      ! interval amplitude in theta and theta_s
delphis = 2.*PI/Nm        ! interval amplitude in phi_s - phi

z = SQRT(14400.*PI*PI)    ! Vacuum impedance  377 ohm	
alpharad = ANGLE_CURVATURE*PI/180. 
r = b/alpharad            ! sphere radius 
lambda = 0.3/F		        
k = 2.*PI/lambda		      

! Integration step
soglia = a/5.	

! Permittivity constant routine
CALL VEGEDIEL(dsdw,amoi,er,ei,PERMITTIVITY)			    
! Definition of constant
molt = (0.,1.)*((k*r*r)/(2.*PI))    ! simplification in k in from SS to Fpq
den_r = k*l*(er - 1)                ! (1)  pag.655
den_i = k*l*ei						          ! (1)  pag.655 
RR = CMPLX(0.,z)/CMPLX(den_r,den_i) ! (1)  pag.655
!
! Determination of the integration point with trapezium method, al least two point must be selected
!
l2 = CEILING(b/soglia)
m2_1 = CEILING(a/soglia)

delta_theta = alpharad/l2        ! Angular amplitude along the meridian
thetamin = ASIN(a/(2*PI*r))      ! Theta value where the sheet begin
thetamax = thetamin + alpharad   ! Finishing theta

dphi = 10                        ! Angular step amplitude for the rotation of the leaf
nphi = 360/dphi		               ! Number	of leaf position in Azimuth

jLoop: DO j = 1,NIJ			         ! j = incidence angular interval index

  thetaj = (j - 1)*delthet		   ! lower limit of angular interval in theta   
  ke(j,:) = 0.				           ! initialization to compute extinction coefficients

  icodLoop: DO icode = 1,2
    ! icode = 1  : upper half - space scattering
    ! icode = 2  : lower half - space scattering

    iLoop: DO i = 1,NIJ          ! i = scattering angular interval index
      thetais = (i - 1)*delthet  ! lower limit of angular interval in theta_s
      ! Initialization to compute scatter matrix
      xff(1:NM1,:,:) = 0.	   

      k1Loop: DO k1 = 1,2	       ! integration within jth theta interval using Gauss - Legendre technique
        theta = thetaj + .5*delthet*(1. + chi(k1))
        stheta = SIN(theta)
        ctheta = COS(theta)
        k2Loop: DO k2 = 1,2	    ! integration within ith theta_s interval using Gauss - Legendre technique

          thetas = thetais + .5*delthet*(1. + chi(k2))
          IF(icode==2) thetas = PI - thetas ! changing of reference systems [1] e [4]
          sthetas = SIN(thetas)
          cthetas = COS(thetas)
          abscthetas = ABS(cthetas)
          o0 = om(k1)*om(k2)/4. ! to be used in the double integration in theta and theta_s
          !
          ! The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor, 
          ! to be used in formulas (3) and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0
          !
          ! Computations for 0 < phi_s - phi < 180
          k3Loop: DO k3 = 1,NM1

            phis = (k3 - 1)*delphis     ! phi_s - phi	
            phis = phis + PI			      ! changing of reference systems  [1] e [4]
            sphis = SIN(phis)
            cphis = COS(phis)
            !
            ! Initialization of scattering polarization vector
            !
            es(1,1) = -cthetas*cphis
            es(2,1) = -cthetas*sphis
            es(3,1) = sthetas
            es(1,2) = -sphis
            es(2,2) = cphis
            es(3,2) = 0.
            !	
            ! I suppose 36 position in azimuth plane
					  !
            iphiLoop: DO i_phi = 1,nphi

              phi_l = (i_phi - 1)*dphi*PI/180.	    ! Centre azimuth of the leaf
              theta_f = thetamin - delta_theta      ! Minimum value in elevation

              ! S Initialization  
              SS(:,:) = (0.,0.)
              !
              ! Integration respect to theta
              !        
              elevation: DO i_elev = 0,l2

                theta_f = theta_f + delta_theta    

                m2 = CEILING(SIN(thetamax)/SIN(theta_f))
                IF(m2 < m2_1) m2 = m2_1

                ! Determination of angular amplitude of title along the parallel
                delta_phi = a/(m2*r*SIN(theta_f))	
                phi_f = phi_l - a/(2.*r*SIN(theta_f)) - delta_phi 

                ! Sl (Azimuthal integral) and Sm (point value) initialization
                Sl(:,:) = (0.,0.)
                Sm(:,:) = (0.,0.)	

                azimuth: DO i_azim = 0,m2   

                  ! Integration respect to phi
                  phi_f = phi_f + delta_phi
                  !
                  ! The following procedure has the task to test if the part of considered sheet is in shadow
									! or it is illuminated by the field
                  !
                  !	Test = 1  This is a singularity of the function (Q=0 at denominator), I take the contribution of the next
									!						tile as valid contribution
									!
                  !	Test = 2  The tile is not considered, it is in the shadow part
                  !
                  CALL VERIFICA(theta,theta_f,phi_f,test)

                  IF(test == 1) THEN
                    Sl(:,:) = Sl(:,:) + ck*Sm(:,:)
                    CYCLE
                  ELSEIF(test == 2) THEN
                    CYCLE
                  END IF
                  !   Estimation of S matrix
                  !	
                  CALL TRAPEZIO(k,r,RR,z,theta_f,phi_f,stheta,ctheta,sthetas,cthetas,sphis,cphis,Sm)
                  !	
                  ck = 2
                  IF((i_azim == 0.).OR.(i_azim == m2))  ck = 1
                  Sl(:,:) = Sl(:,:) + ck*Sm(:,:)
                END DO azimuth  

                ! Here is concluded the integration for row, or with fixed theta and phi variable
                Sl(:,:) = .5*delta_phi*Sl(:,:)
                ! Integration in theta
                ck = 2
                IF((i_elev==0).or.(i_elev==l2))  ck = 1

                SS(:,:) = SS(:,:) + ck*Sl(:,:)
              END DO elevation

              ! Constant moltiplication factor (18) di [1]
              SS(:,:) = .5*delta_theta*SS(:,:)
              ! The amplitude scattering function are in [cm]

              fvv = 100.*molt*(SS(1,1)*es(1,1) + SS(2,1)*es(2,1) + SS(3,1)*es(3,1))
              fhv = 100.*molt*(SS(1,1)*es(1,2) + SS(2,1)*es(2,2) + SS(3,1)*es(3,2))
              fvh = 100.*molt*(SS(1,2)*es(1,1) + SS(2,2)*es(2,1) + SS(3,2)*es(3,1))
              fhh = 100.*molt*(SS(1,2)*es(1,2) + SS(2,2)*es(2,2) + SS(3,2)*es(3,2))

              ! Extinction cross - section computation using the forward scattering theorem
              IF(icode==2.and.i==j.and.k2==k1.and.k3==1) THEN
                ! ke values returned to the main program are divided by ctheta as requested in formula (8) of [1]
                cke = om(k1)/(2*ctheta*nphi)
                ke(j,1) = ke(j,1) + cke*ABS(4*PI/k*AIMAG(fvv))*100   ! k in cm^ - 1
                ke(j,2) = ke(j,2) + cke*ABS(4*PI/k*AIMAG(fhh))*100   ! k in cm^ - 1
              END IF
              !
              fvvfvv = fvv*CONJG(fvv)
              fvhfvh = fvh*CONJG(fvh)
              fhvfhv = fhv*CONJG(fhv)
              fhhfhh = fhh*CONJG(fhh)
              !
              xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nphi
              xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nphi
              xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nphi
              xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nphi

            END DO iphiLoop

          END DO k3Loop

        END DO k2Loop

      END DO  k1Loop
      !
      !  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      !
      DO k3 = nm/2 + 2,nm
      k3a = 2*NM1 - k3
      xff(k3,:,:) = xff(k3a,:,:)
      END DO   ! k3
      !
      ! Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
      DO istoki = 1,2
				DO istoks = 1,2
					xx(:) = CMPLX(xff(:,istoks,istoki), 0.)
					CALL CFSTFT(NF,xx)
					! normalization
					xx(1) = xx(1)/NM
					xx(2:) = xx(2:)*2/NM
					s(i,j,icode,1:NM1,istoks,istoki) = REAL(xx(1:NM1))		
				END DO     ! istoks
      END DO     ! istoki

    END DO iLoop

  END DO icodLoop

END DO jLoop

END SUBROUTINE CURVED_SHEET
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE VERIFICA
! CALLED BY: FOGLIA_CURVA
!
SUBROUTINE VERIFICA(theta,theta_f,phi_f,test)
! Dummy variables declaration
REAL, INTENT(IN)				:: theta,theta_f,phi_f
INTEGER, INTENT(OUT)		:: test
! Local variables declaration
REAL                    :: x,QQ,st,ct,stf,ctf,spf,cpf,cp1

st = SIN(theta)
ct = COS(theta)
stf = SIN(theta_f)
ctf = COS(theta_f)
spf = SIN(phi_f)
cpf = COS(phi_f)

! cos(phi_1) is internal to RE and RH
cp1 = st*stf*cpf + ct*ctf

! verification about the singularity of the function
QQ = stf*stf*spf*spf + (st*ctf - stf*cpf*ct)**2

IF((QQ == 0.).or.(cp1 == 0))THEN
  test = 1
  RETURN
END IF

! Verification about the presence of the considered tile in the shadow part of the sphere 
x = stf*cpf*st + ctf*ct

IF(x <= 0.) THEN
  test = 2
  RETURN
END IF

test = 0.

END SUBROUTINE VERIFICA
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE TRAPEZIO
! CALLED BY: FOGLIA_CURVA
!
!
! Tecnica di integrazione numerica a due punti
!
!
SUBROUTINE TRAPEZIO(k,r,RR,z,theta_f,phi_f,st,ct,sts,cts,sps,cps,Sm)

! Dummy variable declaration
COMPLEX, DIMENSION(3,2), INTENT(OUT)	:: Sm
COMPLEX, INTENT(IN)										:: RR
REAL, INTENT(IN)											:: theta_f,phi_f,st,ct,sts,cts,sps,cps
REAL, INTENT(IN)											:: k,r,z

! Local variable declaration
COMPLEX					      :: RH,RE,c_arg,comune

REAL, DIMENSION(3,2)	:: PPI
REAL, DIMENSION(2)		:: A,B

REAL					        :: stf,ctf,spf,cpf,QQ,arg,cp1
INTEGER					      :: pol

! Estimation of sinusoidal function in theta_f and phi_f
stf = SIN(theta_f)
ctf = COS(theta_f)
spf = SIN(phi_f)
cpf = COS(phi_f)
cp1 = st*stf*cpf + ct*ctf

RH = 1/(1 + (2.*RR)/(z*cp1))
RE = 1/(1 + (2.*RR*cp1)/z)

QQ = stf*stf*spf*spf + (st*ctf - stf*cpf*ct)*(st*ctf - stf*cpf*ct)

A(1) = ctf*st - stf*cpf*ct
A(2) = stf*spf
B(1) = stf*stf*cpf*spf*st + stf*ctf*ct*spf
B(2) =  - stf*cpf*ctf + st*ct*(stf*stf*cpf*cpf - ctf*ctf) + 2.*stf*ctf*ct*ct*cpf

! Estimation od PI component for V Polarization
PPI(1,1) = A(1)*RH*(stf*ctf*cpf*ct - ctf*ctf*st - st*stf*stf*spf*spf) - B(1)*RE*ct*stf*spf
PPI(2,1) = A(1)*RH*(stf*stf*spf*cpf*st + stf*ctf*spf*ct) + B(1)*RE*(stf*cpf*ct - st*ctf)
PPI(3,1) = A(1)*RH*(st*stf*ctf*cpf - stf*stf*ct) + B(1)*RE*st*stf*spf

! Estimation od PI component for V Polarization
PPI(1,2) = A(2)*RH*(stf*ctf*cpf*ct - ctf*ctf*st - st*stf*stf*spf*spf) - B(2)*RE*ct*stf*spf
PPI(2,2) = A(2)*RH*(stf*stf*spf*cpf*st + stf*ctf*spf*ct) + B(2)*RE*(stf*cpf*ct - st*ctf)
PPI(3,2) = A(2)*RH*(st*stf*ctf*cpf - stf*stf*ct) + B(2)*RE*st*stf*spf

DO pol = 1,2
  Sm(1,pol) = (1 - sts*sts*cps*cps)*PPI(1,pol) - sts*sts*sps*cps*PPI(2,pol) - sts*cts*cps*PPI(3,pol)
  Sm(2,pol) =  - sts*sts*sps*cps*PPI(1,pol) + (1 - sts*sts*sps*sps)*PPI(2,pol) - sts*cts*sps*PPI(3,pol)
  Sm(3,pol) =  - sts*cts*cps*PPI(1,pol) - sts*cts*sps*PPI(2,pol) + sts*sts*PPI(3,pol)
END DO  

arg =  -stf*cpf*st - ctf*ct - sts*cps*stf*cpf - sts*sps*stf*spf - cts*ctf
c_arg = k*r*CMPLX(0.,arg)
comune = stf*EXP(c_arg)/QQ

Sm(:,:) = Sm(:,:)*comune

END SUBROUTINE TRAPEZIO
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE INFINITE_LENGTH
! CALLED BY: MAIN
! Computes the scatter matrix and the extinction vector using the Infinite-length approximation for an ensemble of cylinders.
! Requires double precision.
!
SUBROUTINE INFINITE_LENGTH(acyl,lcyl,amoi,dsdw,s,ke,select)
! Dummy variables declaration
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)   :: s
REAL, DIMENSION(:,:), INTENT(OUT)		        :: ke
REAL, INTENT(IN)						                :: acyl,lcyl,amoi,dsdw
CHARACTER, INTENT(IN)					              :: select
! Lopcal variables declaration

COMPLEX(8), DIMENSION(0:NMAX0+1)	:: z
COMPLEX(8), DIMENSION(0:NMAX0)		:: hni,hni1,ju,ju1,r,ev,eh,hv,hh,a,b
COMPLEX(8)							          :: li,u,uu,uunn,rr1,rr2,rr3,rr4,rr5,zmin,ff,ec,svv,svh,shv,shh
COMPLEX(8)							          :: fv1,fv2,fh1,fh2,fvv,fvh,fhv,fhh,fvfs,fhfs,ff1vv,ff1vh,ff1hv,ff1hh,fv1fs,fh1fs

COMPLEX, DIMENSION(NM)				    :: xx

REAL(8), DIMENSION(NM,7,2,2)		  :: xff						 
REAL(8), DIMENSION(0:NMAX0)				:: js

REAL, DIMENSION(2)				        :: icos,jcos

REAL(8)								:: lambda,k,kk,ni,nini,mi,ls,conv,delthet,delphis,eta,theta,thetaj,thetais,alpha,beta,gamma
REAL(8)								:: salpha,calpha,sbeta,cbeta,sgamma,cgamma,stheta,ctheta,absctheta,tvi,thi,cthetil,ccthetil
REAL(8)								:: sthetil,sphil,tt,cphil,phil,ure,uim,thetas,sthetas,cthetas,abscthetas,o0,cdirez
REAL(8)								:: phis,sphis,cphis,ths,cthetsl,ccthetsl,sthetsl,tts,argmi,phisg,tvs,cphisl,sphisl
REAL(8)								:: phisl,arg2,dd,fvvfvv,fhhfhh,fhvfhv,fvhfvh,cke,rap1,rap2

REAL							  	:: lcyl2,er,ej,alpha1,dalpha,beta1,dbeta,gamma1,dgamma,aww,sww,ww

INTEGER								:: nabg,i,j,k1,k2,k3,k3a,iflag,imax,ialpha,ibeta,igamma,icode,nmax,n,istoki,istoks
INTEGER								:: nalpha,nbeta,ngamma,icosn

! Variables for library functions
COMPLEX(8), DIMENSION(0:NMAX0 + 1)	:: JB,YB,dJB,dYB,SIG
COMPLEX(8)													:: ETAA,ZLMIN,cni

INTEGER															:: KFN,MODE,JFAIL,JPR

!
DATA jcos/ -1, -1/
DATA icos/1, -1/

! Initiation of input parameter for WCLBES function 
ETAA = (0.d0, 0.d0)
ZLMIN = (0.d0, 0.d0)
KFN = 2
MODE = 1

icosn = 1
aww = 1.


SELECT CASE (select)
CASE ('S')
  ! Stem or Trunks selection
  alpha1 = ALPHA1STE
  nalpha = NALPHASTE
  dalpha = DALPHASTE
  beta1 = BETA1STE
  nbeta = NBETASTE
  dbeta = DBETASTE
  gamma1 = GAMMA1STE
  ngamma = NGAMMASTE
  dgamma = DGAMMASTE
CASE('P','N','R')
  ! Petioles, needle, secondary branches
  alpha1 = ALPHA1PET
  nalpha = NALPHAPET
  dalpha = DALPHAPET
  beta1 = BETA1PET
  nbeta = NBETAPET
  dbeta = DBETAPET
  gamma1 = GAMMA1PET
  ngamma = NGAMMAPET
  dgamma = DGAMMAPET
CASE('E')
  ! Ears
  alpha1 = ALPHA1EAR
  nalpha = NALPHAEAR
  dalpha = DALPHAEAR
  beta1 = BETA1EAR
  nbeta = NBETAEAR
  dbeta = DBETAEAR
  gamma1 = GAMMA1EAR
  ngamma = NGAMMAEAR
  dgamma = DGAMMAEAR
CASE('B')
  ! Branches selection, in this case a weight function is used
  alpha1 = ALPHA1BRA
  nalpha = NALPHABRA
  dalpha = DALPHABRA
  beta1 = BETA1BRA
  nbeta = NBETABRA
  dbeta = DBETABRA
  gamma1 = GAMMA1BRA
  ngamma = NGAMMABRA
  dgamma = DGAMMABRA
  icosn=2
  sww=0.
  DO ibeta = 1,nbeta
	  beta = beta1 + (ibeta - 1)*dbeta
	  sww = sww + COSD(beta - 90)**2
  END DO
  aww = nbeta/sww
END SELECT

nabg = nalpha*nbeta*ngamma
conv = 180.d0/DPI
delthet = DPI/(2*NIJ)	    	! interval amplitude in theta and theta_s
delphis = 2.d0*DPI/nm	    	! interval amplitude in phi_s - phi
!
lambda = 30.d0/F				    ! wavelength, in cm
k = 2.d0*DPI/lambda			    ! wavenumber, in cm^ - 1
kk = k*k
eta = 1.d0					        ! free-space impedance
lcyl2 = lcyl/2.
!
CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)
   	                        ! permittivity computation
ec = CMPLX(er, -ej,8)

jLoop: DO j = 1,NIJ			  	! j = incidence angular interval index
  thetaj = (j - 1)*delthet	! lower limit of angular interval in theta
  !
  ke(j,:) = 0.				    	! initialization to compute extinction coefficients   
  !
  icodeLoop: DO icode = 1,2
    ! icode = 1  : upper half - space scattering
    ! icode = 2  : lower half - space scattering

    iLoop: DO i = 1,NIJ         ! i = scattering angular interval index
    thetais = (i - 1)*delthet		! lower limit of angular interval in theta_s
    xff(:,:,:,:) = 0.d0						      ! scattering function initialization
    iflag = 0						        ! Controls order of Bessel functions series
    imax = 1

    alphaLoop: DO ialpha = 1,nalpha
      alpha = alpha1 + (ialpha - 1)*dalpha
      salpha = DSIND(alpha)
      calpha = DCOSD(alpha)

      betaLoop: DO ibeta = 1,nbeta
        beta = beta1 + (ibeta - 1)*dbeta
        sbeta = DSIND(beta)
        cbeta = DCOSD(beta)

			  IF(icosn == 2) THEN
			    ww = aww*COSD(beta - 90)**2
			  ELSE
			    ww = 1.
			  ENDIF

        gammaLoop: DO igamma = 1,ngamma
          gamma = gamma1 + (igamma - 1)*dgamma
          sgamma = DSIND(gamma)
          cgamma = DCOSD(gamma)

          k1Loop: DO k1 = 1,NS     ! integration within jth theta interval
            ! using Gauss - Legendre technique
            theta = thetaj + .5d0*delthet*(1. + chi(k1))
            stheta = DSIN(theta)
            absctheta = DCOS(theta)
            ctheta = jCOS(icode)*absctheta
            !
            ! computation of enp,hnp (15) of [7].(dependent on incidence angle only)
            !
            tvi = (salpha*sgamma + calpha*sbeta*cgamma)*ctheta - cbeta*cgamma*stheta                    ! (7) of [7]
            thi = salpha*sbeta*cgamma - calpha*sgamma													! (7) of [7]
            cthetil = (salpha*sgamma + calpha*sbeta*cgamma)*stheta + cbeta*cgamma*ctheta                ! (6) of [7]
            ccthetil = cthetil*cthetil
            sthetil = 0.d0
            IF(ccthetil <= 1.d0) sthetil = DSQRT(1.d0 - ccthetil)
            tt = DSQRT(tvi*tvi + thi*thi)
            !
            IF(beta == 0.d0 .AND. gamma == 0.d0)THEN
              cphil = calpha
              sphil =  - salpha
            ELSE
              cphil = (cbeta*calpha*stheta - sbeta*ctheta)/tt    ! (6) of [7]
              sphil = ((sbeta*sgamma*calpha - cgamma*salpha)*stheta + cbeta*sgamma*ctheta)/tt
            END IF
            phil = 0.d0
            IF(cphil < 0.d0) phil = DPI
            IF(DABS(cphil) <= 1.d0) phil = DACOS(cphil)
            IF(sphil < 0.d0)phil = 2.d0*DPI - phil
            li = k*CDSQRT(ec - ccthetil)                      ! (15) of [7]
            u = acyl*li                                       ! (16) of [7]
            ure = REAL(u,8)
            uim = DIMAG(u)
            uu = u*u
            ni = k*acyl*sthetil                               ! (16) of [7]
            nini = ni*ni
            uunn = 1.d0/nini - 1.d0/uu
            !
            ! The number of Bessel functions is increased until  relative error is  <  1/1000 (see control later)
            !
       200  IF(iflag == 0)THEN
              imax = imax + 1
              nmax = 2**imax
              IF(imax == 7) THEN
                nmax = 100
                IFLAG = 1
                WRITE(*,*) "The Bessel functions cannot be computed with further accuracy, the order is 100 in INFINITE_LENGTH"
              ENDIF
            END IF
           
						! First kind of Bessel functions
            CALL WBSJA(u, 0.d0, nmax, 3, ju)

						! First, Second and Third kinds Bessel functions
						cni = (1.d0, 0.d0) * ni
						CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
 					  IF(JFAIL /= 0) STOP "Error in WCLBES function in INFINITE_LENGTH"
						! Hankel second kind functions computation [-NH,NH]
						hni(0:) = JB(0:NMAX0) - (0.d0, 1.d0) * YB(0:NMAX0)
						! Derivates of Hankel functions 
						hni1(:) = dJB(0:NMAX0) - (0.d0, 1.d0) * dYB(0:NMAX0)
            ! Derivates of Bessel functions
            ju1(0) = -ju(1)                         ! first derivatives (9.1.28 of [9])
            DO n = 1,nmax
              ju1(n) = ju(n - 1) - (n/u)*ju(n)      ! first derivatives (9.1.27 of [9])
            END DO

            DO n = 0,nmax
              rr1 = ju1(n)/(u*ju(n))
              rr2 = hni1(n)/(ni*hni(n))
              rr3 = DPI*nini*hni(n)/2.d0
              rr4 = (rr2 - rr1)*(rr2 - ec*rr1)
              rr5 = n*n*ccthetil*uunn*uunn
              r(n) = rr3*(rr4 - rr5)
              ev(n) = (0.d0,1.d0)*sthetil*(rr2 - rr1)/(r(n)*ju(n))           ! (15) of [7]
              hv(n) = (sthetil/eta)*uunn*n*( - cthetil)/(r(n)*ju(n))         ! (15) of [7]
              hh(n) = (0.d0,1.d0)*(sthetil/eta)*(rr2 - ec*rr1)/(r(n)*ju(n))  ! (23) of [7]
              eh(n) =  - sthetil*uunn*n*( - cthetil)/(r(n)*ju(n))            ! (23) of [7]
            END DO

            k2Loop: DO k2 = 1,NS     ! integration within ith theta_s interval using Gauss - Legendre technique
              thetas = thetais + .5d0*delthet*(1. + chi(k2))
              sthetas = DSIN(thetas)
              abscthetas = DCOS(thetas)
              cthetas = iCOS(icode)*abscthetas
              o0 = om(k1)*om(k2)/4.d0  ! to be used in the double integration in theta and theta_s
              ! The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor,
              ! to be used in formulas (3) and (4) of [7].
              cdirez = (delthet*stheta/abscthetas)*o0
              ! Computations for 0 < phi_s - phi < 180

              k3Loop: DO k3 = 1,NM1
                phis = (k3 - 1)*delphis    ! phi_s - phi
                sphis = SIN(phis)
                cphis = COS(phis)
                !
                phisg = phis*conv
                ! (6) and (7) of [7] in the scattering direction
                tvs = (DSIND(alpha - phisg)*sgamma + DCOSD(alpha - phisg)*sbeta*cgamma)*cthetas - cbeta*cgamma*sthetas
                ths = DSIND(alpha - phisg)*sbeta*cgamma - DCOSD(alpha - phisg)*sgamma
                cthetsl = (DSIND(alpha - phisg)*sgamma + DCOSD(alpha - phisg)*sbeta*cgamma)*sthetas + cbeta*cgamma*cthetas
                ccthetsl = cthetsl*cthetsl
                sthetsl = 0.d0
                IF(ccthetsl <= 1.d0) sthetsl = DSQRT(1.d0 - ccthetsl)
                tts = DSQRT(tvs*tvs + ths*ths)

                argmi = k*lcyl2*(cthetsl - cthetil)              ! (28 of [7])
                IF(argmi == 0.d0)THEN
                  mi = 1.d0
                ELSE
                 mi = DSIN(argmi)/argmi
                END IF

                IF(beta == 0.d0 .AND. gamma == 0.d0) THEN
                  cphisl = DCOSD(alpha - phisg)
                  sphisl =  - DSIND(alpha - phisg)
                ELSE
                  cphisl = (cbeta*sthetas*DCOSD(alpha - phisg) - sbeta*cthetas)/tts
                  sphisl = ((sbeta*sgamma*DCOSD(alpha - phisg) - cgamma*DSIND(alpha - phisg))*sthetas + cbeta*sgamma*cthetas)/tts
                END IF
       
                phisl = 0.d0
                IF(cphisl < 0.d0) phisl = DPI
                IF(DABS(cphisl) <= 1.d0) phisl = DACOS(cphisl)
                IF(sphisl  <  0.d0)phisl = 2.d0*DPI - phisl
                !
                ls = k*sthetsl                        ! (31 of [7])
                arg2 = ls*acyl
                CALL DBSJA(arg2, 0.d0, nmax, 3, js)
                z(0) = (acyl/(li*li - ls*ls))*(li*js(0)*ju(1) - ls*ju(0)*js(1))  ! (a10 of [7])
                z(1) = (acyl/(li*li - ls*ls))*(li*js(1)*ju(2) - ls*ju(1)*js(2))  ! (a11 of [7])
                zmin = z(1)
                IF(ls == 0.d0)THEN
                  z(2) = z(0) - 2*ju(1)/li*acyl*0.5
                ELSE
                 z(2) = z(0) - 2*ju(1)*js(1)/(ls*li)  ! (a13 of [7])
                END IF
                DO n = 3,nmax + 1
                  IF(ls == 0.d0)THEN
                    z(n) = z(n - 2)
                  ELSE
                    z(n) = z(n - 2) - 2*(n - 1)*ju(n - 1)*js(n - 1)/(ls*li)  ! (a13 of [7])
                  END IF
                END DO
                a(0) = (k/(2.d0*li))*(zmin - z(1))    !  (38 of [7])
                b(0) = (k/(2.d0*li))*(zmin + z(1))
                !
                ! Computations in 42 of [7]
                !
                ff = kk*lcyl2*mi*(ec - 1.d0)    ! to be used in (42) of [7]
                ff1vv = ev(0)*(b(0)*cthetsl*( - cthetil) - sthetsl*z(0))
                ff1vh = 0.d0
                ff1hv = 0.d0
                ff1hh = b(0)*eta*hh(0)
                DO n = 1,nmax
                  a(n) = (k/(2*li))*(z(n - 1) - z(n + 1))      ! (38 of [7])
                  b(n) = (k/(2*li))*(z(n - 1) + z(n + 1))      ! (38 of [7])
                  svv = 2.d0*((ev(n)*( - cthetil)*b(n) - (0.d0,1.d0)*eta*hv(n)*a(n))*cthetsl - sthetsl*ev(n)*z(n))
                  svh = ((eh(n)*( - cthetil)*b(n) - (0.d0,1.d0)*eta*hh(n)*a(n))*cthetsl - sthetsl*eh(n)*z(n))
                  shv = (eta*hv(n)*b(n) + (0.d0,1.d0)*ev(n)*( - cthetil)*a(n))
                  shh = 2.d0*(eta*hh(n)*b(n) + (0.d0,1.d0)*eh(n)*( - cthetil)*a(n))
                  ff1vv = ff1vv + svv*DCOS(n*(phisl - phil))
                  ff1vh = ff1vh + svh*DSIN(n*(phisl - phil))
                  ff1hv = ff1hv + shv*DSIN(n*(phisl - phil))
                  ff1hh = ff1hh + shh*DCOS(n*(phisl - phil))
                END DO
                ! (42) of [7]
                ff1vv = ff1vv*ff
                ff1vh = 2.d0*(0.d0,1.d0)*ff1vh*ff
                ff1hv = 2.d0*(0.d0,1.d0)*ff1hv*ff
                ff1hh = ff1hh*ff
                !
                ! From local frame to absolute reference frame
                !
                IF(beta==0.d0.and.gamma==0.d0)THEN
                  fvv = ff1vv
                  fvh = ff1vh
                  fhv = ff1hv
                  fhh = ff1hh
                ELSE
                  fv1 =  - ff1vv*tvs + ff1hv*ths       ! (49 of [7])
                  fv2 =  - ff1vh*tvs + ff1hh*ths
                  fh1 =  - ff1vv*ths - ff1hv*tvs
                  fh2 =  - ff1vh*ths - ff1hh*tvs
                  dd = tt*tts
                  fvv = ( - fv1*tvi + fv2*thi)/dd      ! (48 of [7])
                  fvh = ( - fv2*tvi - fv1*thi)/dd
                  fhv = ( - fh1*tvi + fh2*thi)/dd
                  fhh = ( - fh1*thi - fh2*tvi)/dd
                END IF
                !
                ! Computations of scattering functions and averaging over Eulerian angles and angular intervals of theta and theta_s
                !
                fvvfvv = fvv*DCONJG(fvv)
                fvhfvh = fvh*DCONJG(fvh)
                fhvfhv = fhv*DCONJG(fhv)
                fhhfhh = fhh*DCONJG(fhh)
                xff(k3,imax,1,1) = xff(k3,imax,1,1) + ww*cdirez*fvvfvv/nabg
                xff(k3,imax,1,2) = xff(k3,imax,1,2) + ww*cdirez*fvhfvh/nabg
                xff(k3,imax,2,1) = xff(k3,imax,2,1) + ww*cdirez*fhvfhv/nabg
                xff(k3,imax,2,2) = xff(k3,imax,2,2) + ww*cdirez*fhhfhh/nabg
                !
                ! Extinction cross sections (51 of [7])
                !
                IF(icode == 2 .AND. i == j .AND. k2 == k1 .AND. k3 == 1) THEN ! forward direction
                  fv1fs = DCMPLX(REAL(ff1vv,8),DABS(DIMAG(ff1vv)))
                  fh1fs = DCMPLX(REAL(ff1hh,8),DABS(DIMAG(ff1hh)))
                  IF(beta == 0.d0 .AND. gamma == 0.d0) THEN
                    fvfs = fv1fs
                    fhfs = fh1fs
                  ELSE
                    fvfs = (fv1fs*tvs*tvi + fh1fs*ths*thi)/dd          ! (49 of [7])
                    fhfs = (fv1fs*ths*thi + fh1fs*tvs*tvi)/dd          ! (49 of [7])
                  END IF
                  !
                  ! ke values RETURNed to the main program are divided by absctheta
                  ! as requested in formula (8) of [1]
                  cke = om(k1)/(2*absctheta*nabg)
                  ke(j,1) = ke(j,1) + ww*cke*DABS(4.d0*DPI/k*DIMAG(fvfs))
                  ke(j,2) = ke(j,2) + ww*cke*DABS(4.d0*DPI/k*DIMAG(fhfs))
                END IF
                IF(iflag == 0)THEN
                  rap1 = (xff(k3,imax,1,1) - xff(k3,imax - 1,1,1))/xff(k3,imax,1,1)
                  rap2 = (xff(k3,imax,2,2) - xff(k3,imax - 1,2,2))/xff(k3,imax,2,2)
                  IF(DABS(rap1) < 1.d-3 .AND. DABS(rap2) < 1.d-3) THEN
                    iflag = 1
                    imax = imax - 1
                    nmax = 2**imax
                  ELSE
                    IF(icode == 2 .AND. i == j .AND. k2 == k1 .AND. k3 == 1) ke(j,:) = 0.d0
                    IF(imax < 8) GOTO 200
                    GOTO 200
                  END IF
                END IF

              END DO k3Loop
 
            END DO k2Loop

          END DO k1Loop

        END DO gammaLoop

      END DO betaLoop

    END DO alphaLoop
    !  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
    DO k3 = nm/2 + 2,nm
      k3a = 2*NM1 - k3
      xff(k3,imax,:,:) = xff(k3a,imax,:,:)
    END DO
    !
    !      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
    !
    DO istoki = 1,2
      DO istoks = 1,2
        xx(:) = CMPLX(xff(:,imax,istoks,istoki),0.)
        CALL CFSTFT(NF,xx)
        ! normalization
        xx(1) = xx(1) / NM
        xx(2:) = xx(2:)*2/NM
        s(i,j,icode,1:NM1,istoks,istoki) = REAL(xx(1:NM1))	   
      END DO
    END DO

    END DO iLoop

  END DO icodeLoop

END DO jLoop

END SUBROUTINE INFINITE_LENGTH
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE HOLLOW_CYLINDER
! CALLED BY: MAIN
! Computes backscatter coefficient and extinction of an ensemble of "vertical" cylinders (main stems) over the soil.
! Requires double precision.
!
SUBROUTINE HOLLOW_CYLINDER(acili,acile,lstem,amoi,dsdw,s,ke,select)
! Dummy variables declaration
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
REAL, DIMENSION(:,:), INTENT(OUT)					:: ke
REAL, INTENT(IN)											    :: acili,acile,amoi,lstem,dsdw
CHARACTER, INTENT(IN)                     :: select
! Local variables declaration
INTEGER,PARAMETER		:: NH = 2
INTEGER,PARAMETER		:: NI = NH - 1

! Coefficients (12.a) di [1] 
COMPLEX(8), DIMENSION(4,-NI:NI,4,4)	:: Xn
COMPLEX(8), DIMENSION(-NI:NI,4,4)		:: Dn,Fn
COMPLEX(8), DIMENSION(4,4)			    :: Inv2,Inv4,InvFn

! Coefficients (15.b) di [1]
COMPLEX(8), DIMENSION(-NI:NI,4)			:: CoeffV,CoeffH,CoeffVS,CoeffHS
COMPLEX(8), DIMENSION(-NI:NI,2)			:: Ezn,Epn,Hzn,Hpn

! Array contenente i valori delle funzioni di Bessel 
COMPLEX(8), DIMENSION(-NH-1:NH+1)		:: J,Y,dJ,dY,H1,H2,dH1,dH2
COMPLEX(8), DIMENSION(-NI:NI)				:: U1n,U2n,U3n

! Incident polarization vector
REAL(8), DIMENSION(4), PARAMETER		:: PoliVV = (/1.d0, 1.d0, 0.d0, 0.d0/)
REAL(8), DIMENSION(4), PARAMETER		:: PoliHH = (/0.d0, 0.d0,-1.d0,-1.d0/)
REAL(8), DIMENSION(4)								:: PoliV,PoliH

! Propagation dimensions
COMPLEX(8), DIMENSION(2)	:: II,KK
COMPLEX(8)								:: ec,ek,A,B,C,k,kr,kz,fvv,fhv,fvh,fhh
REAL(8)										:: Z0,k0,lambda,V,sincV,y0,cb
REAL											:: er,ej

! Coordinate System Reference variables
REAL(8), DIMENSION(3)	:: ki,ks,hs,vs,x1,y1,z1,zeta
REAL(8)								:: prx,pry

! Angular dimensions
REAL(8) ::	thetaj,thetas,sthetas,cthetas,delthet,Phip,phis,sphis,cphis,delphis,thetais,theta,stheta
REAL(8) ::	ctheta,abscthetas,absctheta,cphit,sphit,phit,cthetap,sthetap
INTEGER	::  nalpha,nbeta,ngamma
REAL    ::  alpha1,dalpha,beta1,dbeta,gamma1,dgamma

! Generic variables
REAL(8), DIMENSION(NM,2,2)	:: Xff

COMPLEX, DIMENSION(NM)			:: XX

REAL(8)		::	o0,cdirez,cke,fvvfvv,fhhfhh,fhvfhv,fvhfvh,rho
INTEGER		::  i,ij,n,i0,i1,i2,icode,k1,k2,k3,k3a,istoki,istoks
LOGICAL		::  test

! Variables for the average values of Eularian angles
REAL(8)		::  alpha,beta,salpha,sbeta,calpha,cbeta
INTEGER	  ::  nabg,ialpha,ibeta

! Variables for CERNILB library functions
COMPLEX(8), DIMENSION(0:NH + 1) :: SIG
COMPLEX(8)											:: ETA,ZLMIN,cni
INTEGER													:: KFN,MODE,JFAIL,JPR

!	Variables for LAPACK library functions
! Dimension is 64*4, where 4x4 is the dimension of matrix to invert 
INTEGER, PARAMETER									:: LWORK = 256
COMPLEX(8), DIMENSION(LWORK)				:: WORK
COMPLEX(8), DIMENSION(4)						:: IPIV
INTEGER															:: INFO

! Initiation of input parameter for WCLBES function 
ETA = (0.d0, 0.d0)
ZLMIN = (0.d0, 0.d0)
KFN = 2
MODE = 1

SELECT CASE (select)
CASE ('S')
  ! Stem selection
  alpha1 = ALPHA1STE
  nalpha = NALPHASTE
  dalpha = DALPHASTE
  beta1 = BETA1STE
  nbeta = NBETASTE
  dbeta = DBETASTE
  gamma1 = GAMMA1STE
  ngamma = NGAMMASTE
  dgamma = DGAMMASTE
CASE ('E')
  ! Ears
  alpha1 = ALPHA1EAR
  nalpha = NALPHAEAR
  dalpha = DALPHAEAR
  beta1 = BETA1EAR
  nbeta = NBETAEAR
  dbeta = DBETAEAR
  gamma1 = GAMMA1EAR
  ngamma = NGAMMAEAR
  dgamma = DGAMMAEAR
END SELECT

nabg = nalpha*nbeta*ngamma

! In order that cylinder is hollow, the internal radius must be great than 0.5 cm
IF (acili .GT. .05) THEN
	test = .TRUE.		! Hollow
	i0 = 1
ELSE
	test = .FALSE.	! Full
	i0 = 3
END IF

Z0 = DSQRT(14400*DPI*DPI)    ! free - space impedance 		
PoliV(:) = PoliVV(:)
PoliH(:) = PoliHH(:)/Z0		   ! Relative magnetic field amplitude
lambda = 30./F						   ! free - space wavelength in m
k0 = 2.d0*DPI/lambda		 		 ! free - space wavenumber in m^ - 1				  
delthet = DPI/NN					   ! interval amplitude in theta and theta_s
delphis = 2.d0*DPI/NM				 ! interval amplitude in phi_s - phi

CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)	! Permittivity constant
ec = DCMPLX(er,ej)

jLoop: DO ij = 1,NIJ				! ij = incidence angular interval index
	thetaj = (ij - 1)*delthet	! lower limit of angular interval in theta   
	ke(ij,:) = 0.d0					    ! initialization to compute extinction coefficients
 
	icodLoop: DO icode = 1,2
		! icode = 1  : upper half - space scattering
		! icode = 2  : lower half - space scattering
		!
		iLoop: DO i = 1,NIJ						! i = scattering angular interval index
			thetais = (i - 1)*delthet		! lower limit of angular interval in theta_s
			! Initialization to compute scatter matrix
			xff(1:NM1,:,:) = 0.d0
      !                  (averaging over scatterer orientation)
			Lalpha: DO ialpha = 1,nalpha       
				alpha = alpha1 + (ialpha - 1)*dalpha
				salpha = DSIND(alpha)
				calpha = DCOSD(alpha)

				Lbeta: DO ibeta = 1,nbeta		   
					beta = beta1 + (ibeta - 1)*dbeta
					sbeta = DSIND(beta)
					cbeta = DCOSD(beta)
					!
					! Unit normal vector into the locar reference system
					z1(1) = calpha*sbeta
					z1(2) = salpha*sbeta
					z1(3) = cbeta		

					k1Loop: DO k1 = 1,NS			! integration within jth theta interval using Gauss - Legendre technique

						! theta angle into the the generic reference system
						theta = DPI - (thetaj + .5d0*delthet*(1.d0 + chi(k1)))
						stheta = DSIN(theta)
						ctheta = DCOS(theta)
						absctheta = DABS(ctheta)
						!
						! incidence vector
						ki(1) =  - stheta
						ki(2) = 0.d0
						ki(3) = ctheta
						!
						! Unit normal vector into the locar reference system
						y1(:) = CROSS(ki,z1)
						!	
						! Unit normal vector into the locar reference system
						x1(:) = CROSS(y1,z1)
						!
						!	Phip angle into the local reference system
						!
						prx = DOT_PRODUCT(ki,x1)
						pry = DOT_PRODUCT(ki,y1)
						Phip = DATAN(pry/prx) 
						IF(prx < 0.d0) THEN 
							Phip = DPI + Phip
						ELSEIF(pry < 0.d0) THEN
							Phip = 2.d0*DPI + Phip
						END IF
						!
						!	Theta angle into the local reference system
						!
						cthetap = DOT_PRODUCT(z1,ki)
						sthetap = DSQRT(1 - cthetap*cthetap)
						!
						! If the cyclynder is full the cycle start from 2, else 1
						NInterfaccia: DO i1 = i0,4		

							SELECT CASE (i1)
								CASE (1)
									rho = acili
									ek = (1.d0,0.d0)			
								CASE (2)
									rho = acili
									ek = ec
								CASE (3)
									rho = acile
									ek = ec
								CASE DEFAULT
									rho = acile
									ek = (1.d0,0.d0)
							END SELECT

							k = k0*CDSQRT(ek)

							kz = k*cthetap
							kr = k*sthetap

							! First and second types of Bessel functions are computed for positive order
							cni = (1.d0, 0.d0)*kr*rho
							CALL WCLBES(cni,ETA,ZLMIN,NH,J(0:NH+1),Y(0:NH+1),dJ(0:NH+1),dY(0:NH+1),SIG,KFN,MODE,JFAIL,JPR)
							IF(JFAIL /= 0) STOP "Error in WCLBES function in HOLLOW_CYLINDER"
							! Computation of negative orders [-NH,-1] rif. [2]
							DO n = NH,1,-1
								J(-n) = J(n)*(-1.d0,0.d0)**n
								Y(-n) = Y(n)*(-1.d0,0.d0)**n
								dJ(-n) = dJ(n)*(-1.d0,0.d0)**n
								dY(-n) = dY(n)*(-1.d0,0.d0)**n
							END DO

							! Hankel functions computation [-NH,NH]
							H1(:) = J(:) + (0.d0, 1.d0) * Y(:)
							H2(:) = J(:) - (0.d0, 1.d0) * Y(:)
							dH1(:) = dJ(:) + (0.d0, 1.d0) * dY(:)
							dH2(:) = dJ(:) - (0.d0, 1.d0) * dY(:)

							! Multiplicative constants of matrix elements
							B =  -(0.d0,1.d0)*k0*Z0/kr
							C = (0.d0,1.d0)*k0*ek/(kr*Z0)

							Ordine: DO i2 = -NI,NI
								A = -(i2*kz)/(rho*kr*kr)

								! First row
								Xn(i1,i2,1,1) = H1(i2)
								Xn(i1,i2,1,2) = H2(i2)
								Xn(i1,i2,1,3) = (0.d0,0.d0)
								Xn(i1,i2,1,4) = (0.d0,0.d0)

								! Second row
								Xn(i1,i2,2,1) = (0.d0,0.d0)
								Xn(i1,i2,2,2) = (0.d0,0.d0)
								Xn(i1,i2,2,3) = H1(i2)
								Xn(i1,i2,2,4) = H2(i2)

								! Third row
								Xn(i1,i2,3,1) = A*H1(i2)
								Xn(i1,i2,3,2) = A*H2(i2)
								Xn(i1,i2,3,3) = B*dH1(i2)
								Xn(i1,i2,3,4) = B*dH2(i2)

								! Fourth row 
								Xn(i1,i2,4,1) = C*dH1(i2)
								Xn(i1,i2,4,2) = C*dH2(i2)
								Xn(i1,i2,4,3) = A*H1(i2)
								Xn(i1,i2,4,4) = A*H2(i2)
							END DO Ordine
						END DO NInterfaccia

						Ordine2: DO i1 = -NI,NI

							! Inverse matrix of Xn32
							Inv4(:,:) = Xn(4,i1,:,:)
							! LU factorization
							CALL ZGETRF(4, 4, Inv4, 4, IPIV, INFO)
							IF(INFO /= 0) STOP "Error during Xn32 matrix factorization in HOLLOW function"
							! Matrix invertion
							CALL ZGETRI(4, Inv4, 4, IPIV, WORK, LWORK, INFO)
							IF(INFO /= 0) STOP "Error during Xn32 matrix inversion in HOLLOW function"
						
							IF (Test .EQ. .FALSE.) THEN
								! Full cylinder
								! Computation of matrix Dn (15.a) in [1]  (Xn32^ - 1*Xn22)
								Dn(i1,:,:) = MATMUL(Inv4,Xn(3,i1,:,:))
							ELSE	
								! Hollow cylinder	
								! Inverse matrix of Xn21			
								Inv2(:,:) = Xn(2,i1,:,:)		
								! LU factorization
								CALL ZGETRF(4, 4, Inv2, 4, IPIV, INFO)
								IF(INFO /= 0) STOP "Error during Xn21 matrix factorization in HOLLOW function"
								! Matrix invertion
								CALL ZGETRI(4, Inv2, 4, IPIV, WORK, LWORK, INFO)
								IF(INFO /= 0) STOP "Error during Xn21 matrix inversion in HOLLOW function"
								! Computation of matrix Dn (15.a) in [1]  (Xn32^ - 1*Xn22)*(Xn21^ - 1*Xn11)
								Dn(i1,:,:) = MATMUL(MATMUL(Inv4,Xn(3,i1,:,:)),MATMUL(Inv2,Xn(1,i1,:,:)))
							END IF

							! alpha1 and alpha2
							Fn(i1,1,1) = ( - 1.d0,0.d0)
							Fn(i1,1,2) = Dn(i1,1,1) + Dn(i1,1,2)
							Fn(i1,1,3) = Dn(i1,1,3) + Dn(i1,1,4)
							Fn(i1,1,4) = (0.d0,0.d0)

							! alpha3 and alpha4
							Fn(i1,2,1) = (0.d0,0.d0)
							Fn(i1,2,2) = Dn(i1,2,1) + Dn(i1,2,2)
							Fn(i1,2,3) = Dn(i1,2,3) + Dn(i1,2,4)
							Fn(i1,2,4) = (0.d0,0.d0)

							! alpha5 and alpha6
							Fn(i1,3,1) = (0.d0,0.d0)
							Fn(i1,3,2) = Dn(i1,3,1) + Dn(i1,3,2)
							Fn(i1,3,3) = Dn(i1,3,3) + Dn(i1,3,4)
							Fn(i1,3,4) = ( - 1.d0,0.d0)

							! alpha7 and alpha8
							Fn(i1,4,1) = (0.d0,0.d0)
							Fn(i1,4,2) = Dn(i1,4,1) + Dn(i1,4,2)
							Fn(i1,4,3) = Dn(i1,4,3) + Dn(i1,4,4)
							Fn(i1,4,4) = (0.d0,0.d0)

							! Computation of inverse Fn
							InvFn(:,:) = Fn(i1,:,:)		
							! LU factorization
							CALL ZGETRF(4, 4, InvFn, 4, IPIV, INFO)
							IF(INFO /= 0) STOP "Error during Fn matrix factorization in HOLLOW function"
							! Matrix invertion
							CALL ZGETRI(4, InvFn, 4, IPIV, WORK, LWORK, INFO)
							IF(INFO /= 0) STOP "Error during Fn matrix inversion in HOLLOW function"

							A =  - .5d0*stheta*(0.d0,1.d0)**i1

							! Computation of scattering coefficients Ans, Cns, An1 and Cn1 (15.b) [1]
							CoeffVS(i1,:) = A*MATMUL(InvFn,PoliV)
							CoeffHS(i1,:) = A*MATMUL(InvFn,PoliH)

							! The conventiom adopted in the first vector at second member in (10) of [1] 
							CoeffVS(i1,3) = CoeffVS(i1,4)
							CoeffVS(i1,2) = (0.d0,0.d0)
							CoeffVS(i1,4) = (0.d0,0.d0)

							CoeffHS(i1,3) = CoeffHS(i1,4)
							CoeffHS(i1,2) = (0.d0,0.d0)
							CoeffHS(i1,4) = (0.d0,0.d0)

							! Computation of coefficients Ank,Bnk,Cnk e Dnk, linked to the exteral interface field (10) of [1]
							CoeffV(i1,:) = CoeffVS(i1,:) + A*PoliV(:)
							CoeffH(i1,:) = CoeffHS(i1,:) + A*PoliH(:)

							B = CDEXP(-(0.d0,1.d0)*i1*Phip)   

							CoeffV(i1,:) = B*MATMUL(Xn(4,i1,:,:),CoeffV(i1,:)) 	  ! Vertical polarization (11) of [1]
							Ezn(i1,1) = CoeffV(i1,1)
							Hzn(i1,1) = Z0*CoeffV(i1,2)
							Epn(i1,1) = CoeffV(i1,3)
							Hpn(i1,1) = Z0*CoeffV(i1,4)

							CoeffH(i1,:) = B*MATMUL(Xn(4,i1,:,:),CoeffH(i1,:))	  ! Horizontal polarization (11) of [1]
							Ezn(i1,2) = CoeffH(i1,1)
							Hzn(i1,2) = Z0*CoeffH(i1,2)
							Epn(i1,2) = CoeffH(i1,3)
							Hpn(i1,2) = Z0*CoeffH(i1,4)

						END DO Ordine2

						k2Loop: DO k2 = 1,NS			! integration within ith theta_s interval using Gauss - Legendre technique
							thetas = thetais + .5d0*delthet*(1.d0 + chi(k2))
							IF(icode.EQ.2) thetas = DPI - thetas
							sthetas = DSIN(thetas)
							cthetas = DCOS(thetas)
							abscthetas = DABS(cthetas)
							o0 = om(k1)*om(k2)/4.		! to be used in the double integration in theta and theta_s
							!
							! The "s" functions returned to the main program include the (delthet*stheta/abscthetas) factor, 
							! to be used in formulas (3) and (4) of [1].
							cdirez = (delthet*stheta/abscthetas)*o0
							!
							! Computations for phij = DPI and DPI<phis<2*DPI
							!
							k3Loop: DO k3 = 1,NM1  
								phis = DPI + (k3 - 1)*delphis			
								sphis = DSIN(phis)
								cphis = DCOS(phis)
								!
								! Unit vector Ks				  
								ks(1) = sthetas*cphis
								ks(2) = sthetas*sphis
								ks(3) = cthetas				  
								!
								! Unit vector Hs
								zeta(:) = (/0.d0,0.d0,1.d0/)
								hs(:) = CROSS(zeta,ks)
								!
								! Unit vector Vs
								vs(:) = CROSS(hs,ks)
								!					
								prx = DOT_PRODUCT(ks,x1)
								pry = DOT_PRODUCT(ks,y1)
								!
								! From this point the formulation used is tha same adopted in [3] 
								CB = DSQRT(prx*prx + pry*pry)
								!
								phit = DATAN(pry/prx) 
								IF(prx < 0.d0) THEN 
									phit = DPI + phit
								ELSEIF(pry < 0.d0) THEN
									phit = 2.d0*DPI + phit
								END IF
								cphit = DCOS(phit)
								sphit = DSIN(phit)
								y0 = k0*acile*CB

								! First and second types of Bessel functions are computed for positive order
								cni = (1.d0, 0.d0)*k0*acile*CB
								CALL WCLBES(cni,ETA,ZLMIN,NH,J(0:NH+1),Y(0:NH+1),dJ(0:NH+1),dY(0:NH+1),SIG,KFN,MODE,JFAIL,JPR)
								IF(JFAIL /= 0) STOP "Error in WCLBES function in HOLLOW_CYLINDER"
								! Computation of negative orders [-NH,-1] rif. [2]
								DO n = NH,1,-1
									J(-n) = J(n)*(-1.d0,0.d0)**n
									dJ(-n) = dJ(n)*(-1.d0,0.d0)**n
								END DO

								! Multiplicative coefficients (29 - 30) of [3]
								II(:) = (0.d0,0.d0)
								KK(:) = (0.d0,0.d0)				

								DO n = -NI,NI
									! Multiplicative coefficients (31) of [3]
									A = 2.d0*DPI*CDEXP((0.d0,1.d0)*n*phit)*(0.d0, - 1.d0)**n					  
									U1n(n) = A*J(n)
									U2n(n) = A*((0.d0,1.d0)*cphit*dJ(n) + sphit*J(n)*n/y0)
									U3n(n) = A*((0.d0,1.d0)*sphit*dJ(n) - cphit*J(n)*n/y0)

									! Coefficients (29 - 30) of [3]
									II(:) = II(:) + (Hpn(n,:)*U1n(n) + (DOT_PRODUCT(ks,z1)/CB)*Hzn(n,:)*(sphit*U2n(n) - cphit*U3n(n)) - (Ezn(n,:)/CB)*(sphit*U3n(n) + cphit*U2n(n)))
									KK(:) = KK(:) - (Epn(n,:)*U1n(n) + (DOT_PRODUCT(ks,z1)/CB)*Ezn(n,:)*(sphit*U2n(n) - cphit*U3n(n)) + (Hzn(n,:)/CB)*(sphit*U3n(n) + cphit*U2n(n)))
								END DO

								V = .5*k0*lstem*DOT_PRODUCT((ki(:) - ks(:)),z1)
								IF(V == 0.d0) THEN
									sincV = 1.d0
								ELSE
									sincV = DSIN(V)/V
								END IF
					
								A = (0.d0,1.d0)*acile*lstem*sincV*k0/(4.*DPI)
								fvv = A*(DOT_PRODUCT(vs,z1)*II(1) + DOT_PRODUCT(hs,z1)*KK(1))
								fvh = A*(DOT_PRODUCT(vs,z1)*II(2) + DOT_PRODUCT(hs,z1)*KK(2))
								fhv = A*(DOT_PRODUCT(hs,z1)*II(1) - DOT_PRODUCT(vs,z1)*KK(1))
								fhh = A*(DOT_PRODUCT(hs,z1)*II(2) - DOT_PRODUCT(vs,z1)*KK(2))		

								! Extinction cross - section computation using the forward scattering theorem
								IF(icode == 2 .and. i == ij .and. k2 == k1 .and. k3 == 1) THEN
									! ke values returned to the main program are divided by ctheta as requested in formula (8) of [1]
									cke = om(k1)/(2*absctheta)
									ke(ij,1) = ke(ij,1) + cke*4.*DPI*DIMAG(fvv)/(k0*nabg)  
									ke(ij,2) = ke(ij,2) + cke*4.*DPI*DIMAG(fhh)/(k0*nabg)   
								END IF

								fvvfvv = fvv*DCONJG(fvv)
								fvhfvh = fvh*DCONJG(fvh)
								fhvfhv = fhv*DCONJG(fhv)
								fhhfhh = fhh*DCONJG(fhh)

								XFF(k3,1,1) = XFF(k3,1,1) + cdirez*fvvfvv/nabg
								XFF(k3,1,2) = XFF(k3,1,2) + cdirez*fvhfvh/nabg
								XFF(k3,2,1) = XFF(k3,2,1) + cdirez*fhvfhv/nabg
								XFF(k3,2,2) = XFF(k3,2,2) + cdirez*fhhfhh/nabg

							END DO k3Loop
						END DO k2Loop
					END DO  k1Loop
				END DO Lbeta      
			END DO Lalpha        
			!
			! 180<phi_s - phi<360 values (azimuthal symmetry is assumed)
			!
			DO k3 = nm/2 + 2,nm
				k3a = 2*NM1 - k3
				XFF(k3,:,:) = XFF(k3a,:,:)
			END DO   ! k3
			!
			! Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
			DO istoki = 1,2
				DO istoks = 1,2
					xx(:) = CMPLX(xff(:,istoks,istoki),0.d0)
					CALL CFSTFT(NF,xx)
					xx(1) = xx(1) / NM
					xx(2:) = xx(2:)*2/NM
					S(i,ij,icode,1:NM1,istoks,istoki) = REAL(xx(1:NM1))		
				END DO     ! istoks
			END DO     ! istoki

		END DO iLoop
	END DO icodLoop
END DO jLoop

CONTAINS
	!
	! The function compute the vectorial cross product
	!
	FUNCTION CROSS(A,B)
	! Function declaration
	REAL(8), DIMENSION(3) 						 :: CROSS

	!Dummy arguments
	REAL(8), DIMENSION(:), INTENT(IN) :: A,B
	!Local variables declaration
	REAL(8) mod

		CROSS(1) = A(2)*B(3) - A(3)*B(2)
		CROSS(2) = A(3)*B(1) - A(1)*B(3)
		CROSS(3) = A(1)*B(2) - A(2)*B(1)

		mod = DSQRT(CROSS(1)*CROSS(1) + CROSS(2)*CROSS(2) + CROSS(3)*CROSS(3))

		CROSS(:) = CROSS(:)/mod

	END FUNCTION CROSS
END SUBROUTINE HOLLOW_CYLINDER
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE TRUNKS_ABS
! CALLED BY: MAIN
! Computes the absorption coefficient ensemble of "vertical" cylinders (main stems) over the soil.
! Requires double precision.
!
SUBROUTINE TRUNKS_ABS(acil,length,amoist,dsdw,Ka,Kes)
!Dummy arguments declaration
REAL,INTENT(IN)										        :: length,acil,amoist,dsdw
REAL,DIMENSION(:,:),INTENT(OUT)		        :: Ka
REAL,DIMENSION(:,:),INTENT(OUT), OPTIONAL	:: Kes
!Local variables declaration
COMPLEX(8),DIMENSION(-NMAX0:NMAX0)		:: R,Ev,Eh,Hv,Hh,cv,ch,dv,dh
COMPLEX(8),DIMENSION(0:NMAX0)					:: A,B

COMPLEX(8),DIMENSION(0:NMAX0+1)				:: z

COMPLEX(8)							:: li,u,uu,uunn,rr1,rr2,rr3,rr4,rr5,zmen,ff,ec,svv,svh,shv,shh,ff1vv
COMPLEX(8)							:: ff1vh,ff1hv,ff1hh,fv1,fv2,fh1,fh2,fvv,fvh,fhv,fhh

REAL(8),DIMENSION(9,2,2)						:: Xtoro
REAL(8),DIMENSION(-NMAX0-1:NMAX0+1)	:: Zlic
REAL(8),DIMENSION(0:NMAX0+1)				:: Js
REAL(8),DIMENSION(NN,2)							:: Ke
REAL(8),DIMENSION(2)								:: Jcos,Icos

REAL(8)		:: dalfa,alfa1,dbeta,beta1,gamma1,dgamma,conv,sqrt2,c1,eta
REAL(8)		:: ccabs,tetaj,alfa,salfa,calfa,beta,sbeta,cbeta,gamma,sgamma,cgamma,teta,steta,abscteta,cteta
REAL(8)		:: tvi,thi,ctetil,cctetil,stetil,tt,cfil,sfil,fil,ure,uim,tetas,stetas,ctetas,absctetas,o0,cdirez
REAL(8)		:: fis,sfis,cfis,fisg,tvs,ths,ctetsl,cctetsl,stetsl,tts,argmi,cfisl,sfisl,fisl,arg2,dd
REAL(8)		:: fvvfvv,fhvfhv,fvhfvh,fhhfhh,rap1,rap2,cka,ni,deltet,lcil
REAL(8)		:: kav1,kav2,kav3,kah1,kah2,kah3,kav,kah,lambda,k,kk,nini,mi,ls

REAL			::	er, ej

INTEGER		:: nalfa,nbeta,ngamma,nabg,icode,iflag,imax,ialfa,ibeta,igamma,k1,nmax,n,k3,j

! Variables for library functions
COMPLEX(8), DIMENSION(-NMAX0-1 : NMAX0 + 1)	:: JU,JB,YB,dJU,dJB,dYB,Hn,dHn,SIG
COMPLEX(8)																	:: ETAA,ZLMIN,cni

INTEGER																			:: KFN,MODE,JFAIL,JPR

! Initiation of input parameter for WCLBES function 
ETAA = (0.d0, 0.d0)
ZLMIN = (0.d0, 0.d0)
KFN = 2
MODE = 1

DATA Jcos/-1,-1/
DATA Icos/1,-1/ 

alfa1 = ALPHA1STE
nalfa = NALPHASTE
dalfa = DALPHASTE
beta1 = BETA1STE
nbeta = NBETASTE
dbeta = DBETASTE
gamma1 = GAMMA1STE
ngamma = NGAMMASTE
dgamma = DGAMMASTE
	
nabg=nalfa*nbeta*ngamma

conv=180.d0/DPI
deltet=DPI/(2*nij)
sqrt2=DSQRT(2.d0)
c1=2.d0/DSQRT(DPI)
lambda=30.d0/f
k=2.d0*DPI/lambda
kk=k*k
eta=1.d0   ! e' ininfluente

CALL VEGEDIEL(dsdw,amoist,er,ej,PERMITTIVITY)
ec = CMPLX(er, -ej, 8)

lcil = length/2
ccabs=4.*k*DPI*lcil*ej

toroloop: DO j=1,nij       ! toro di incidenza
	tetaj=(j-1)*deltet
	ke(j,:)=0.d0
	ka(j,:)=0.d0
	kav=0.; kav1=0.; kav2=0.; kav3=0.
	kah=0.; kah1=0.; kah2=0.; kah3=0.
	!<--------------------------------------------------loop di icode

	icodeLoop: DO icode=2,2
		! icode=1  : scattering semispazio superiore
		! icode=2  : scattering semispazio inferiore
		!     
		!<--------------------------------------------------loop di i
		xtoro(3:,:,:)=0.d0
		iflag=0               !indice controllo nmax
		imax=3
		! --------------------------------- loop alfa,beta,gamma

		alfaLoop: DO ialfa=1,nalfa
			alfa=alfa1+(ialfa-1)*dalfa
			salfa=DSIND(alfa)
			calfa=DCOSD(alfa)

			betaLoop: DO ibeta=1,nbeta
				beta=beta1+(ibeta-1)*dbeta
				sbeta=DSIND(beta)
				cbeta=DCOSD(beta)

				gammaLoop: DO igamma=1,ngamma
					gamma=gamma1+(igamma-1)*dgamma
					sgamma=DSIND(gamma)
					cgamma=DCOSD(gamma)

					tetaLoop: DO k1=1,2     ! integrale (gauss) in teta
						teta=.5d0*(deltet*chi(k1)+2.d0*tetaj+deltet)
						!teta=tetaj
						steta=DSIN(teta)
						abscteta=dCOS(teta)
						cteta=jCOS(icode)*abscteta
						! calcolo enp,hnp (15),(16). Dipendono solo dall'angolo di incidenza; (fi=0).
						tvi=(salfa*sgamma+calfa*sbeta*cgamma)*cteta-cbeta*cgamma*steta      ! (7)
						thi=salfa*sbeta*cgamma-calfa*sgamma                                 ! (7)
						ctetil=(salfa*sgamma+calfa*sbeta*cgamma)*steta+cbeta*cgamma*cteta   ! (6)
						cctetil=ctetil*ctetil
						stetil=0.d0
						IF(cctetil.le.1.d0) stetil=DSQRT(1.d0-cctetil)
						tt=DSQRT(tvi*tvi+thi*thi)
						! steta e' sempre diverso da 0!!
						IF(beta.eq.0.d0.and.gamma.eq.0.d0)THEN
							cfil=calfa
							sfil=-salfa
						ELSE
							cfil=(cbeta*calfa*steta-sbeta*cteta)/tt             ! (6)
							sfil=((sbeta*sgamma*calfa-cgamma*salfa)*steta+cbeta*sgamma*cteta)/tt
						ENDIF
						fil=0.d0
						IF(cfil.lt.0.d0)fil=DPI
						IF(DABS(cfil).le.1.d0) fil=DACOS(cfil)
						IF(sfil.lt.0.d0) fil=2.d0*DPI-fil 
						li=k*CDSQRT(ec-cctetil)
						u=acil*li
						ure=DREAL(u)
						uim=DIMAG(u)
						uu=u*u
						ni=k*acil*stetil
						nini=ni*ni
						uunn=1.d0/nini-1.d0/uu

       200  IF(iflag == 0)THEN
              imax = imax + 1
              nmax = 2**imax
              IF(imax == 7) THEN
                nmax = 100
                IFLAG = 1
                WRITE(*,*) "The Bessel functions cannot be computed with further accuracy, the order is 100 in TRUNKS_ABS"
              ENDIF
            END IF

						CALL WCLBES(u,ETAA,ZLMIN,nmax,JU(0:nmax+1),YB(0:nmax+1),dJU(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
						IF(JFAIL /= 0) STOP "Error in WCLBES function in TRUNK_ABS "
						! Computation of negative orders [-nmax,-1] rif. [2]

          	! First and second types of Bessel functions are computed for positive order
						cni = (1.d0, 0.d0) * ni
						CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
						IF(JFAIL /= 0) STOP "Error in WCLBES function in TRUNK_ABS "
						! Computation of negative orders [-nmax,-1] rif. [2]
						DO n = nmax,1,-1
							JU(-n) = JU(n)*(-1.d0,0.d0)**n
							JB(-n) = JB(n)*(-1.d0,0.d0)**n
							YB(-n) = YB(n)*(-1.d0,0.d0)**n
							dJU(-n) = dJU(n)*(-1.d0,0.d0)**n
							dJB(-n) = dJB(n)*(-1.d0,0.d0)**n
							dYB(-n) = dYB(n)*(-1.d0,0.d0)**n
						END DO
						! Hankel functions computation [-nmax,nmax]
						Hn(:) = JB(:) - (0.d0, 1.d0) * YB(:)
						dHn(:) = dJB(:) - (0.d0, 1.d0) * dYB(:)

						DO n=-nmax,nmax
							rr1=dJU(n)/(u*JU(n))
							rr2=dHn(n)/(ni*Hn(n))
							rr3=DPI*nini*Hn(n)/2.d0
							rr4=(rr2-rr1)*(rr2-ec*rr1)
							rr5=n*n*cctetil*uunn*uunn
							r(n)=rr3*(rr4-rr5)
							ev(n)=(0.d0,1.d0)*stetil*(rr2-rr1)/(r(n)*JU(n))
							hv(n)=(stetil/eta)*uunn*n*(-ctetil)/(r(n)*JU(n))
							hh(n)=(0.d0,1.d0)*(stetil/eta)*(rr2-ec*rr1)/(r(n)*JU(n))
							eh(n)=-stetil*uunn*n*(-ctetil)/(r(n)*JU(n))
						END DO  
						tetas=teta
						stetas=DSIN(tetas)
						absctetas=dCOS(tetas)
						ctetas=iCOS(icode)*absctetas
						o0=1.d0
						!cdirez=(deltet*steta/absctetas)*o0 
						cdirez=1.
						!
						! cdirez e' il coefficiente per il calcolo del contributo di una coppia di direzioni alla matrice 
						! di scatter di una coppia di direzioni alla matrice di scatter di una coppia di tori (vd. (2) di appunto).                   
						! L'integrale e' calcolato con la quadratura di gauss, il fattore deltet**2 si elide perche' sono 
						! valori medi. Cdirezg e' per la sigma del terreno (v. appunto emissivita')
						!                                                                 
						fisLoop: DO k3=1,1      ! calcoli per fis=0                    
							fis=0.
							sfis=SIN(fis)     
							cfis=COS(fis)     
							fisg=fis*conv
							! espressioni analoghe alle (6) e (7) per direzione di scattering.
							tvs=(DSIND(alfa-fisg)*sgamma+DCOSD(alfa-fisg)*sbeta*cgamma)*ctetas-cbeta*cgamma*stetas
							ths=DSIND(alfa-fisg)*sbeta*cgamma-DCOSD(alfa-fisg)*sgamma
							ctetsl=(DSIND(alfa-fisg)*sgamma+DCOSD(alfa-fisg)*sbeta*cgamma)*stetas+cbeta*cgamma*ctetas
							cctetsl=ctetsl*ctetsl
							stetsl=0.d0
							IF(cctetsl.le.1.d0)  stetsl=DSQRT(1.d0-cctetsl)
							tts=DSQRT(tvs*tvs+ths*ths)
							argmi=k*lcil*(ctetsl-ctetil)              ! (28)

							IF(argmi.eq.0.d0)THEN
								mi=1.d0
							ELSE
								mi=DSIN(argmi)/argmi
							ENDIF

							IF(beta.eq.0.d0.and.gamma.eq.0.d0)THEN
								cfisl=DCOSD(alfa-fisg)
								sfisl=-DSIND(alfa-fisg)
							ELSE
								cfisl=(cbeta*stetas*DCOSD(alfa-fisg)-sbeta*ctetas)/tts
								sfisl=((sbeta*sgamma*calfa-cgamma*DSIND(alfa-fisg))*stetas+cbeta*sgamma*ctetas)/tts
							ENDIF
							fisl=0.d0
							IF(cfisl.lt.0.d0) fisl=DPI
							IF(DABS(cfisl).le.1.d0) fisl=DACOS(cfisl)
							IF(sfisl.lt.0.d0) fisl=2.d0*DPI-fisl
							ls=k*stetsl                               ! (31)
							arg2=ls*acil

							! First kind of Bessel functions
							CALL DBSJA(arg2, 0.d0, nmax, 3, js)

							z(0)=(acil/(li*li-ls*ls))*(li*js(0)*JU(1)-ls*JU(0)*js(1))  ! (a10)
							z(1)=(acil/(li*li-ls*ls))*(li*js(1)*JU(2)-ls*JU(1)*js(2))  ! (a11)
							!zmen=(acil/(li*li-ls*ls))*(-li*js(1)*JU(0)+ls*JU(1)*js(0))  ! (a11)+abram.,9.1.5
							zmen=z(1)
							IF(ls.eq.0.d0)THEN
								z(2)=z(0)-2*JU(1)/li*acil*0.5
							ELSE
								z(2)=z(0)-2*JU(1)*js(1)/(ls*li)  ! (a13)
							ENDIF
							DO n=3,nmax+1
								IF(ls.eq.0.d0)THEN
									z(n)=z(n-2)
								ELSE
									z(n)=z(n-2)-2*(n-1)*JU(n-1)*js(n-1)/(ls*li)  ! (a13)
								ENDIF
							END DO 
							a(0)=(k/(2.d0*li))*(zmen-z(1))  !  (38)
							b(0)=(k/(2.d0*li))*(zmen+z(1))
							!                                  calcolo (42)
							ff=kk*lcil*mi*(ec-1.d0)               
							ff1vv=ev(0)*(b(0)*ctetsl*(-ctetil)-stetsl*z(0))
							ff1vh=0.d0
							ff1hv=0.d0
							ff1hh=b(0)*eta*hh(0)
							DO n=1,nmax
								a(n)=(k/(2*li))*(z(n-1)-z(n+1))      ! (38)
								b(n)=(k/(2*li))*(z(n-1)+z(n+1))      ! (38)
								svv=2.d0*((ev(n)*(-ctetil)*b(n)-(0.d0,1.d0)*eta*hv(n)*a(n))*ctetsl-stetsl*ev(n)*z(n))
								svh=((eh(n)*(-ctetil)*b(n)-(0.d0,1.d0)*eta*hh(n)*a(n))*ctetsl-stetsl*eh(n)*z(n))
								shv=(eta*hv(n)*b(n)+(0.d0,1.d0)*ev(n)*(-ctetil)*a(n))
								shh=2.d0*(eta*hh(n)*b(n)+(0.d0,1.d0)*eh(n)*(-ctetil)*a(n))
								ff1vv=ff1vv+svv*DCOS(n*(fisl-fil))
								ff1vh=ff1vh+svh*DSIN(n*(fisl-fil))
								ff1hv=ff1hv+shv*DSIN(n*(fisl-fil))
								ff1hh=ff1hh+shh*DCOS(n*(fisl-fil))
							END DO  
							ff1vv=ff1vv*ff
							ff1vh=2.d0*(0.d0,1.d0)*ff1vh*ff
							ff1hv=2.d0*(0.d0,1.d0)*ff1hv*ff
							ff1hh=ff1hh*ff

							IF(beta.eq.0.d0.and.gamma.eq.0.d0)THEN
								fvv=ff1vv
								fvh=ff1vh
								fhv=ff1hv
								fhh=ff1hh
							ELSE
								!  Cambiamento di riferimento
								fv1=-ff1vv*tvs+ff1hv*ths       ! (49)
								fv2=-ff1vh*tvs+ff1hh*ths
								fh1=-ff1vv*ths-ff1hv*tvs
								fh2=-ff1vh*ths-ff1hh*tvs
								dd=tt*tts    
								fvv=(-fv1*tvi+fv2*thi)/dd      ! (48)
								fvh=(-fv2*tvi-fv1*thi)/dd
								fhv=(-fh1*tvi+fh2*thi)/dd
								fhh=(-fh1*thi-fh2*tvi)/dd
							ENDIF
							fvvfvv=fvv*DCONJG(fvv)
							fvhfvh=fvh*DCONJG(fvh)
							fhvfhv=fhv*DCONJG(fhv)
							fhhfhh=fhh*DCONJG(fhh)
							xtoro(imax,1,1)=xtoro(imax,1,1)+cdirez*fvvfvv/nabg
							xtoro(imax,1,2)=xtoro(imax,1,2)+cdirez*fvhfvh/nabg
							xtoro(imax,2,1)=xtoro(imax,2,1)+cdirez*fhvfhv/nabg
							xtoro(imax,2,2)=xtoro(imax,2,2)+cdirez*fhhfhh/nabg
							ke(j,1)=ke(j,1)+DABS(DIMAG(fvv))*om(k1)/(2.d0*abscteta)
							ke(j,2)=ke(j,2)+DABS(DIMAG(fhh))*om(k1)/(2.d0*abscteta)
							DO n=-nmax+1,nmax-1
								zlic(n)=DIMAG(li*JU(n+1)*CONJG(JU(n)))/(2.*DREAL(li)*DIMAG(li))
								zlic(n)=acil*zlic(n)
							END DO
							DO n=-nmax+2,nmax-2
								cv(n)=k/(2*li)*((0.,1)*ctetil*ev(n)+hv(n))
								ch(n)=k/(2*li)*((0.,1)*ctetil*eh(n)+hh(n))
								dv(n)=k/(2*li)*(-(0.,1)*ctetil*ev(n)+hv(n))
								dh(n)=k/(2*li)*(-(0.,1)*ctetil*eh(n)+hh(n))
							END DO
							kav1=0.
							kav2=0.
							kav3=0.
							kah1=0.
							kah2=0.
							kah3=0.
							DO n=-nmax+2,nmax-2
								kav1=kav1+ev(n)*CONJG(ev(n))*zlic(n)
								kav2=kav2+cv(n)*CONJG(cv(n))*zlic(n+1)
								kav3=kav3+dv(n)*CONJG(dv(n))*zlic(n-1)
								kah1=kah1+eh(n)*CONJG(eh(n))*zlic(n)
								kah2=kah2+ch(n)*CONJG(ch(n))*zlic(n+1)
								kah3=kah3+dh(n)*CONJG(dh(n))*zlic(n-1)
							END DO
							kav=kav1+2.*kav2+2.*kav3
							kah=kah1+2.*kah2+2.*kah3

			        IF(iflag.eq.0)THEN
			          rap1=(xtoro(imax,1,1)-xtoro(imax-1,1,1))/xtoro(imax,1,1)
			          rap2=(xtoro(imax,2,2)-xtoro(imax-1,2,2))/xtoro(imax,2,2)
			          IF(DABS(rap1).lt.1.d-3.and.DABS(rap2).lt.1.d-3)THEN
			            iflag=1
			            imax=imax-1
			            nmax=2**imax
			  	      ELSE
				          ke(j,:)=0.d0
				          kav1=0.; kav2=0.; kav3=0.
				          kah1=0.; kah2=0.; kah3=0.
                  IF(imax < 8) GOTO 200
				        ENDIF
			        ENDIF 

						END DO fisLoop
						cka=om(k1)/(2.*abscteta)
						kav=kav*ccabs*cka/nabg
						kah=kah*ccabs*cka/nabg 
						ka(j,1)=ka(j,1)+(tvi*tvi*kav+thi*thi*kah)/(tt*tt)
						ka(j,2)=ka(j,2)+(thi*thi*kav+tvi*tvi*kah)/(tt*tt)

					END DO tetaLoop

				END DO  gammaLoop
			
			END DO betaLoop

		END DO  alfaLoop
	
	END DO icodeLoop

	ke(j,1)=ke(j,1)*4*PI/(k*nabg)
	ke(j,2)=ke(j,2)*4*PI/(k*nabg)

END DO toroLoop
IF(PRESENT(Kes(:,:))) kes = ke (:,:)
RETURN

END SUBROUTINE TRUNKS_ABS
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUBROUTINE CORNER
! CALLED BY: MAIN
! Computes scattering coefficients in forward scattering, backscattering and in
! specular direction
!
SUBROUTINE CORNER(avert,hvert,amoi,dsdw,sigmae,sigmad,sigmau)
! Dummey arguments declaration
REAL, DIMENSION(:,:,:), INTENT(out)	:: sigmad,sigmau
REAL, DIMENSION(:,:), INTENT(out)		:: sigmae
REAL, INTENT(IN)										:: amoi,avert,hvert,dsdw
! Local variables declaration
REAL, DIMENSION(7,2,2)							:: ssigmad, ssigmau
REAL, DIMENSION(7,2)								:: ssigmae

COMPLEX(8), DIMENSION(0:NMAX0 + 2)	:: ju
COMPLEX(8), DIMENSION(0:NMAX0 + 1)	:: z
COMPLEX(8), DIMENSION(0:NMAX0)			:: hni,hni1,ju1,r,ev,eh,hv,hh,a,b

COMPLEX(8)  :: fv1fs,fh1fs,fvfs,fhfs,li,u,uu,uunn,rr1,rr2,rr3,rr4,rr5,zmeno,ff,ff0
COMPLEX(8)  :: ff1vv0,ff1vh0,ff1hv0,ff1hh0,svv,svh,shv,shh,ff1vv,ff1vh,ff1hv,ff1hh
COMPLEX(8)  :: fv1,fv2,fh1,fh2,fvv,fvh,fhv,fhh,ec

REAL(8), DIMENSION(0:NMAX0 + 2) :: js

REAL(8)     :: lcil,lambda,k,kk,ni,nini,mi,ls,teta,alfa,beta,gamma,salfa,calfa,sbeta,cbeta,sgamma,cgamma
REAL(8)			:: conv,deltet,delfis,eta,roh2o,ak1,ak1k1,tetagraz,abscteta,steta,tetas,ak1steta,cteta,tvi,thi
REAL(8)			:: ctetil,cctetil,stetil,tt,cfil,sfil,fil,ure,uim,stetas,ctetas,absctetas,fis,sfis,cfis
REAL(8)			:: fisg,tvs,ths,ctetsl,cctetsl,stetsl,tts,cfisl,sfisl,fisl,arg2,arg,dd,cks,cka,ccv,cch
REAL(8)			:: fvvfvv,fhvfhv,fvhfvh,fhhfhh,fvvj,fhhj,alfa1,dalfa,beta1,dbeta,gamma1,dgamma

REAL, DIMENSION(2)	:: jcos,icos

REAL				:: er,ej

INTEGER			:: ialfa,ibeta,igamma,nalfa,nbeta,ngamma,nabg,j,imax,nmax,n,icode,k3

LOGICAL			:: test

! Variables for library functions
COMPLEX(8), DIMENSION(0:NMAX0 + 1)	:: JB,YB,dJB,dYB,SIG
COMPLEX(8)													:: ETAA,ZLMIN,cni

INTEGER															:: KFN,MODE,JFAIL,JPR

! Initiation of input parameter for WCLBES function 
ETA = (0.d0, 0.d0)
ZLMIN = (0.d0, 0.d0)
KFN = 2
MODE = 1

DATA jcos/ -1, -1/
DATA icos/1, -1/

alfa1 = ALPHA1STE
nalfa = NALPHASTE
dalfa = DALPHASTE
beta1 = BETA1STE
nbeta = NBETASTE
dbeta = DBETASTE
gamma1 = GAMMA1STE
ngamma = NGAMMASTE
dgamma = DGAMMASTE
	
nabg=nalfa*nbeta*ngamma
conv = 180.d0/DPI
deltet = DPI/(2*nij)
delfis = 2.d0*DPI/nm
lambda = 30.d0/f
k = 2.d0*DPI/lambda
kk = k*k
eta = 1.d0  
roh2o = 1.d0 
lcil = hvert/2.

CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)     ! cylinder permittivity
ec = DCMPLX(er, -ej)
ff0 = kk*lcil*(ec - 1.d0)
ak1 = (2.*DPI)/lambda           
ak1k1 = ak1*ak1

jloop: DO j = 1,nij
	teta = deltet*(j - .5)
	tetagraz = DPI/2. - teta
	abscteta = DCOS(teta)
	steta = DSIN(teta)
	tetas = teta
	ak1steta = ak1*steta
	cteta = jcos(1)*abscteta
  test = .FALSE.
  imax = 2
	bessel: DO WHILE(imax <=  7 .and. test  ==  .false.)
	  imax = imax + 1
	  
    IF(imax == 7) THEN
      nmax = 100
      WRITE(*,*) "The Bessel functions cannot be computed with further accuracy, the order is 100 in CORNER"
    ELSE
			nmax = 2**imax
		ENDIF
		ssigmae(imax,:) = 0.d0
		ssigmad(imax,:,:) = 0.d0
		ssigmau(imax,:,:) = 0.d0

		alfaloop: DO ialfa = 1,nalfa
			alfa = alfa1 + (ialfa - 1)*dalfa
			salfa = DSIND(alfa)
			calfa = DCOSD(alfa)

			betaloop: DO ibeta = 1,nbeta
				beta = beta1 + (ibeta - 1)*dbeta
				sbeta = DSIND(beta)
				cbeta = DCOSD(beta)

				gammaloop: DO igamma = 1,ngamma
					gamma = gamma1 + (igamma - 1)*dgamma
					sgamma = DSIND(gamma)
					cgamma = DCOSD(gamma)
					! evaluation of enp,hnp (15), (16).
					! tey depend only from incidence angle
					tvi = (salfa*sgamma + calfa*sbeta*cgamma)*cteta - cbeta*cgamma*steta
					thi = salfa*sbeta*cgamma - calfa*sgamma
					ctetil = (salfa*sgamma + calfa*sbeta*cgamma)*steta + cbeta*cgamma*cteta
					cctetil = ctetil*ctetil
					stetil = DSQRT(1.d0 - cctetil)
					tt = DSQRT(tvi*tvi + thi*thi)
					IF(beta == 0.d0.and.gamma == 0.d0)THEN
						cfil = calfa
						sfil =  - salfa
					ELSE
						cfil = (cbeta*calfa*steta - sbeta*cteta)/tt
						sfil = ((sbeta*sgamma*calfa - cgamma*salfa)*steta + cbeta*sgamma*cteta)/tt
					ENDIF
					fil = DACOS(cfil)
					IF(sfil < 0.d0)fil = 2.d0*DPI - fil
					li = k*CDSQRT(ec - cctetil)
					u = avert*li
					ure = DREAL(u)
					uim = DIMAG(u)
					uu = u*u
					ni = k*avert*stetil
					nini = ni*ni
					uunn = 1.d0/nini - 1.d0/uu
					! First kind of Bessel functions
					CALL WBSJA(u, 0.d0, nmax, 3, ju)
					! First, Second and Third kinds Bessel functions

					cni = (1.d0, 0.d0) * ni
					CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
					IF(JFAIL /= 0) STOP "Error in WCLBES function in CORNER"
					! Hankel second kind functions computation [0,NMAX0]
					hni(:) = JB(0:NMAX0) - (0.d0, 1.d0) * YB(0:NMAX0)
					! Derivates of Hankel functions 
					hni1(:) = dJB(0:NMAX0) - (0.d0, 1.d0) * dYB(0:NMAX0)
					! Derivates of Bessel functions 
					ju1(0) = -ju(1)                    ! first derivatives (9.1.28 of [9])
					DO n = 1,nmax
						ju1(n) = ju(n - 1) - (n/u)*ju(n) ! first derivatives (9.1.27 of [9])
					END DO

					DO n = 0,nmax
						rr1 = ju1(n)/(u*ju(n))
						rr2 = hni1(n)/(ni*hni(n))
						rr3 = DPI*nini*hni(n)/2.d0
						rr4 = (rr2 - rr1)*(rr2 - ec*rr1)
						rr5 = n*n*cctetil*uunn*uunn
						r(n) = rr3*(rr4 - rr5)
						ev(n) = (0.d0,1.d0)*stetil*(rr2 - rr1)/(r(n)*ju(n))
						hv(n) = (stetil/eta)*uunn*n*( -ctetil)/(r(n)*ju(n))
						hh(n) = (0.d0,1.d0)*(stetil/eta)*(rr2 - ec*rr1)/(r(n)*ju(n))
						eh(n) =  -stetil*uunn*n*( -ctetil)/(r(n)*ju(n))
					END DO 
					! icode = 1  : upward semispace scattering
					! icode = 2  : downward semispace scattering
					icodeloop: DO icode = 2,1, - 1
						!     
						stetas = DSIN(tetas)
						absctetas = DCOS(tetas)
						ctetas = icos(icode)*absctetas

						fisloop: DO k3 = 1,2         !nm1,nm1 - 1      ! calcoli per 0<fis<180
							fis = (k3 - 1)*DPI                 
							sfis = DSIN(fis)
							cfis = DCOS(fis)
							fisg = fis*conv
							tvs = (DSIND(alfa - fisg)*sgamma + DCOSD(alfa - fisg)*sbeta*cgamma)*ctetas - cbeta*cgamma*stetas
							ths = DSIND(alfa - fisg)*sbeta*cgamma - DCOSD(alfa - fisg)*sgamma
							ctetsl = (DSIND(alfa - fisg)*sgamma + DCOSD(alfa - fisg)*sbeta*cgamma)*stetas + cbeta*cgamma*ctetas
							cctetsl = ctetsl*ctetsl
							stetsl = DSQRT(1.d0 - cctetsl)
							tts = DSQRT(tvs*tvs + ths*ths)

							IF(beta == 0.d0.and.gamma == 0.d0)THEN
								cfisl = DCOSD(alfa - fisg)
								sfisl =  - DSIND(alfa - fisg)
							ELSE
								cfisl = (cbeta*stetas*DCOSD(alfa - fisg) - sbeta*ctetas)/tts
								sfisl = ((sbeta*sgamma*calfa - cgamma*DSIND(alfa - fisg))*stetas + cbeta*sgamma*ctetas)/tts
							ENDIF
							fisl = DACOS(cfisl)
							IF(sfisl < 0.d0)fisl = 2.d0*DPI - fisl

							ls = k*stetsl                               ! (31)
							arg2 = ls*avert
	            CALL DBSJA(arg2, 0.d0, nmax, 3, js)

							z(0) = (avert/(li*li - ls*ls))*(li*js(0)*ju(1) - ls*ju(0)*js(1))  ! (a10)
							z(1) = (avert/(li*li - ls*ls))*(li*js(1)*ju(2) - ls*ju(1)*js(2))  ! (a11)
							zmeno = z(1)
							IF(ls == 0.d0)THEN
								z(2) = z(0) - 2*ju(1)/li*avert*0.5
							ELSE
								z(2) = z(0) - 2*ju(1)*js(1)/(ls*li)  ! (a13)
							ENDIF
					
							DO n = 3,nmax + 1
					
								IF(ls == 0.d0)THEN
									z(n) = z(n - 2)
								ELSE
									z(n) = z(n - 2) - 2*(n - 1)*ju(n - 1)*js(n - 1)/(ls*li)  ! (a13)
								ENDIF
							END DO 
							a(0) = (k/(2.d0*li))*(zmeno - z(1))  !  (38)
							b(0) = (k/(2.d0*li))*(zmeno + z(1))
							!                                  calcolo (42)
							ff1vv0 = ev(0)*(b(0)*ctetsl*( - ctetil) - stetsl*z(0))    
							ff1vh0 = 0.d0
							ff1hv0 = 0.d0
							ff1hh0 = b(0)*eta*hh(0)
							DO n = 1,nmax
								a(n) = (k/(2*li))*(z(n - 1) - z(n + 1))      ! (38)
								b(n) = (k/(2*li))*(z(n - 1) + z(n + 1))      ! (38)
								svv = 2.d0*((ev(n)*( - ctetil)*b(n) - (0.d0,1.d0)*eta*hv(n)*a(n))*ctetsl - stetsl*ev(n)*z(n)) 
								svh = ((eh(n)*( - ctetil)*b(n) - (0.d0,1.d0)*eta*hh(n)*a(n))*ctetsl - stetsl*eh(n)*z(n))
								shv = (eta*hv(n)*b(n) + (0.d0,1.d0)*ev(n)*( - ctetil)*a(n))
								shh = 2.d0*(eta*hh(n)*b(n) + (0.d0,1.d0)*eh(n)*( - ctetil)*a(n))
								ff1vv0 = ff1vv0 + svv*DCOS(n*(fisl - fil))
								ff1vh0 = ff1vh0 + svh*DSIN(n*(fisl - fil))
								ff1hv0 = ff1hv0 + shv*DSIN(n*(fisl - fil))
								ff1hh0 = ff1hh0 + shh*DCOS(n*(fisl - fil))
							END DO ! n
							arg = k*lcil*(ctetsl - ctetil)              ! (28)

							IF(arg == 0.d0)THEN
								mi = 1.d0
							ELSE
								mi = DSIN(arg)/arg
							ENDIF
							ff = ff0*mi
							ff1vv = ff1vv0*ff
							ff1vh = 2.d0*(0.d0,1.d0)*ff1vh0*ff
							ff1hv = 2.d0*(0.d0,1.d0)*ff1hv0*ff
							ff1hh = ff1hh0*ff
							!
							fv1 =  -ff1vv*tvs + ff1hv*ths
							fv2 =  -ff1vh*tvs + ff1hh*ths
							fh1 =  -ff1vv*ths - ff1hv*tvs
							fh2 =  -ff1vh*ths - ff1hh*tvs
							dd = tt*tts
							IF(beta == 0.d0.and.gamma == 0.d0)THEN
								fvv = ff1vv
								fvh = ff1vh
								fhv = ff1hv
								fhh = ff1hh
							ELSE
								fvv = ( -fv1*tvi + fv2*thi)/dd
								fvh = ( -fv2*tvi - fv1*thi)/dd
								fhv = ( -fh1*tvi + fh2*thi)/dd
								fhh = ( -fh1*thi - fh2*tvi)/dd
							ENDIF
							!
							fvvfvv = fvv*DCONJG(fvv)
							fvhfvh = fvh*DCONJG(fvh)
							fhvfhv = fhv*DCONJG(fhv)
							fhhfhh = fhh*DCONJG(fhh)
							IF(k3  ==  2) THEN
								IF(icode  ==  2) THEN
									ssigmad(imax,1,1) = ssigmad(imax,1,1) + 4*DPI*fvvfvv/nabg
									ssigmad(imax,2,1) = ssigmad(imax,2,1) + 4*DPI*fhvfhv/nabg
									ssigmad(imax,1,2) = ssigmad(imax,1,2) + 4*DPI*fvhfvh/nabg
									ssigmad(imax,2,2) = ssigmad(imax,2,2) + 4*DPI*fhhfhh/nabg
								ELSE
									cks = 4*DPI/(nabg*abscteta)
									ssigmau(imax,1,1) = ssigmau(imax,1,1) + cks*fvvfvv
									ssigmau(imax,2,1) = ssigmau(imax,2,1) + cks*fhvfhv
									ssigmau(imax,1,2) = ssigmau(imax,1,2) + cks*fvhfvh
									ssigmau(imax,2,2) = ssigmau(imax,2,2) + cks*fhhfhh
								END IF
							END IF

							IF(k3 == 1 .and. icode == 2) THEN
								fv1fs = DCMPLX(DREAL(ff1vv),DABS(DIMAG(ff1vv)))
								fh1fs = DCMPLX(DREAL(ff1hh),DABS(DIMAG(ff1hh)))
								IF(beta == 0. .and. gamma == 0.) THEN
									fvfs = fv1fs
									fhfs = fh1fs
								ELSE
									fvfs = (fv1fs*tvs*tvi + fh1fs*ths*thi)/dd          ! (49)
									fhfs = (fv1fs*ths*thi + fh1fs*tvs*tvi)/dd          ! (49)
								ENDIF
								!
								fvvj = DIMAG(fvfs)
								fhhj = DIMAG(fhfs)
								cka = (4*DPI/k)/(abscteta*nabg)
								ssigmae(imax,1) = ssigmae(imax,1) + cka*DABS(fvvj)
								ssigmae(imax,2) = ssigmae(imax,2) + cka*DABS(fhhj)
							END IF
						END DO fisloop
					END DO icodeloop
				ENDDO gammaloop
			END DO betaloop
		ENDDO alfaloop
		!
		!
		IF(imax > 3) THEN
			ccv = (ssigmae(imax,1) - ssigmae(imax - 1,1))/ssigmae(imax,1)
			cch = (ssigmae(imax,2) - ssigmae(imax - 1,2))/ssigmae(imax,2)
			IF(ccv < .01 .and. cch < .01) test = .TRUE.
		END IF
	END DO bessel

	sigmae(j,:) = ssigmae(imax,:)
	sigmad(j,:,:) = ssigmad(imax,:,:)
	sigmau(j,:,:) = ssigmau(imax,:,:)

END DO jloop

RETURN

END SUBROUTINE CORNER

END MODULE MOD_VEGETATION_FUNCTIONS