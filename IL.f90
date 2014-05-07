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