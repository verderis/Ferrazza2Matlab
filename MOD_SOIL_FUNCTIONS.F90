! *************************************************************************************
! MOD_SOIL_FUNCTIONS contains the algorithms that compute the matrix scattering
! and extinction of soil.
!
! IEM
! *************************************************************************************
MODULE MOD_SOIL_FUNCTIONS

USE MOD_GLOBAL_PARAMETERS
USE MOD_PERMITTIVITY

IMPLICIT NONE

PUBLIC IEM               ,&
			 GEOMETRICALOPTICS ,&
       LITTER 
CONTAINS
!
!******************************************************************************
!
SUBROUTINE IEM(sigmaz,lcorr,vmoist,Sgg,rr,litter_gmoist,biomass)
! Dummy arguments declaration 
REAL, INTENT(IN)														  :: sigmaz,vmoist,lcorr
REAL, DIMENSION(:,:,:,:,:),INTENT(OUT)			  :: Sgg
REAL, OPTIONAL, DIMENSION(:,:),INTENT(OUT)  	:: RR
REAL, OPTIONAL, INTENT(IN)                    :: litter_gmoist, biomass
! Local variable delaration
COMPLEX, DIMENSION(NM)	:: GX

COMPLEX					:: eta2,ff,rhh0,rvv0,egc,fvv,fhh,f1vv,f1hh,fhv,fvh,f1hv,f1vh,f1hhs,f1vvs,f1hvs,f1vhs
COMPLEX					:: tn,tnm,th,thm,tp,tm,r,sq,sqs,ihhn,ihvn,ivhn,ivvn,rpa,ror,ihh2,ivv2,ihh1,ivv1
COMPLEX					:: r21,r22,rt,rii,ctet2,ivh2,ihv2,ihh3,ivv3,ivh3,ihv3,ivh1,ihv1

REAL(8)					:: Sommahh,Sommavv,Sommahv,Sommavh

REAL, DIMENSION(NM,2,2)	:: Gtoro
REAL, DIMENSION(2,2)		:: Gdirez
REAL, DIMENSION(64)			:: FFss

REAL		:: kx,ky,kz,ksx,ksy,ksz,kdx,kdy,wn,lambda,k,kk,llcorr,fi,NM1,conv,delfis,sqrt2,ak1k1,eta1
REAL		:: alfag,betag,ppg,qqg1,cgg,tetaj,tetais,teta,steta
REAL		:: cteta,abscteta,steta2,cteta2,ssteta,ccteta,tetagraz,ak1steta,qqg,tet2,c1,ak1,stet2,r12,beckspiz1
REAL		:: beckspiz2,beckspiz,tteta,tteta2,r11,tetas,stetas,ctetas,sstetas,cctetas,absctetas,o0,cdirez,cdirezg
REAL		:: fis,sfis,cfis,s,ss,cs,css,sf,sf1,csf,c2,c1s,c2s,sigma1,a,b,sigmahv,sigmavh,sigmavv,sigmahh,WN1,WN2
REAL		:: deltet

INTEGER	:: nsup,k1,k2,k3,n,i,j,istoki,istoks

LOGICAL :: RR_present,Litter_present

IF(PRESENT(RR)) THEN
	RR_present = .TRUE.
ELSE
	RR_present = .FALSE.
ENDIF

Litter_present = .FALSE.
IF(PRESENT(litter_gmoist) .AND. PRESENT(biomass)) THEN
  IF(biomass > 0.) Litter_present = .TRUE.
END IF

fi = 0
NM1 = NM/2 + 1
conv = 180./PI
deltet = PI/(2*NIJ)
delfis = 2.*PI/NM
sqrt2 = SQRT(2.)
c1 = 2./SQRT(PI)
nsup = 32
lambda = 30./f
k = 2.*PI/lambda
kk = k*k
ak1 = (2.*PI)/lambda 
ak1k1 = ak1*ak1
eta1 = 1.

IF(Litter_present) THEN
  CALL LITTER(vmoist,litter_gmoist,biomass,egc) 
ELSE
  CALL DOBSOIL(vmoist,egc) 
END IF
eta2 = CSQRT(1/egc)
alfag = k*AIMAG(CSQRT(egc))
betag = k*REAL(CSQRT(egc))
ppg = 2.*alfag*betag
qqg1 = betag*betag - alfag*alfag
cgg = CABS((eta1 - eta2)/(eta1 + eta2))
cgg = cgg*cgg

llcorr = lcorr*lcorr
toroi: DO j = 1,NIJ           ! Incidence loop
	tetaj = (j - 1)*deltet

	IF(RR_present) RR(j,:) = 0.

	toros: DO i = 1,NIJ         ! Scattering loop
		tetais = (i - 1)*deltet
		Gtoro(:,:,:) = 0.

		tetaloop: DO k1 = 1,2     ! integrale (gauss) in teta
			teta = .5*(deltet*Chi(k1) + 2*tetaj + deltet)
			steta = SIN(teta)
			cteta = COS(teta)
			abscteta = COS(teta)
			steta2 = steta*steta
			cteta2 = cteta*cteta
			ssteta = steta2
			ccteta = cteta2
			IF(i .eq. 1) THEN
				tetagraz = PI/2. - teta
				ak1steta = ak1*steta
				qqg = qqg1 - ak1steta*ak1steta
				tet2 = ATAN(ak1steta*sqrt2/SQRT(SQRT(ppg*ppg + qqg*qqg) + qqg))
				ctet2 = COS(tet2)
				ctet2 = CSQRT(1 - (steta*steta)/egc)
				stet2 = SIN(tet2)
				r11 = eta1*abscteta
				r12 = eta1*ctet2
				r21 = eta2*abscteta
				r22 = eta2*ctet2
				rt = (r21 - r12)/(r21 + r12)
				rii = (r11 - r22)/(r11 + r22)
				beckspiz1 = 4*PI*sigmaz*SIN(tetagraz)/lambda
				beckspiz2 =  beckspiz1*beckspiz1
				beckspiz = EXP( - beckspiz2)
				IF(RR_present) THEN
					rr(j,1) = rr(j,1) + (OM(k1)/2.)*beckspiz*rii*CONJG(rii)
					rr(j,2) = rr(j,2) + (OM(k1)/2.)*beckspiz*rt*CONJG(rt)
				ENDIF
				tteta = steta/abscteta
				tteta2 = tteta*tteta
			END IF
			ff = CSQRT(egc - ssteta)

			tetasloop: DO k2 = 1,2       ! integrale (gauss) in tetas
				tetas = .5*(deltet*Chi(k2) + 2*tetais + deltet)
				stetas = SIN(tetas)
				ctetas = COS(tetas)
				sstetas = stetas*stetas
				cctetas = ctetas*ctetas
				absctetas = COS(tetas)
				o0 = OM(k1)*OM(k2)/4.
				cdirez = (deltet*steta/absctetas)*o0
				cdirezg = cdirez/(4.*PI)

				fisloop: DO k3 = 1,NM1      ! calcoli per 0<fis<180
					fis = (k3 - 1)*delfis
					sfis = SIN(fis)
					cfis = COS(fis)
					! ***********************************************
					! Soil bistatic scattering coefficient
					! ***********************************************
					kx = k*SIN(teta)*COS(fi)
					ky = k*SIN(teta)*SIN(fi)
					kz = k*COS(teta)
					ksx = k*SIN(tetas)*COS(fis)
					ksy = k*SIN(tetas)*SIN(fis)
					ksz = k*COS(tetas)
					kdx = ksx - kx
					kdy = ksy - ky
					! ************ pag 206,211 fung**********
					s = SIN(teta)
					ss = SIN(tetas)
					cs = COS(teta)
					css = COS(tetas)
					sf = SIN(fis - fi)
					sf1 = SIN(fi - fis)
					csf = COS(fis - fi)
					sq = CSQRT(egc - (s**2))
					c1 = (csf - s*ss)/(sq*css)
					c2 = s*(ss - s*csf)/css
					sqs = CSQRT(egc - (ss**2))
					c1s = (csf - s*ss)/(sqs*cs)
					c2s = ss*(s - ss*csf)/cs
					! ********** coeff. di fresnel *********
					rhh0 = (cs - sq)/(cs + sq)
					rvv0 = (egc*cs - sq)/(egc*cs + sq)
					r = (rvv0 - rhh0)/2
					tn = 1 + rvv0
					tnm = 1 - rvv0
					th = 1 + rhh0
					thm = 1 - rhh0
					tp = 1 + r
					tm = 1 - r
					! ******** kirchhoff field coefficent ********
					rpa = rvv0
					ror = rhh0
					fvv = ((2*rpa)/(cs + css))*(s*ss - (1 + cs*css)*csf)
					fhh =  - ((2*ror)/(cs + css))*(s*ss - (1 + cs*css)*csf)
					fhv = 2*r*sf
					fvh = 2*r*sf1
					! ******** complementary field coefficent ********
					f1hhs = (css*thm - sqs*th)*(th*csf + thm*c1s) - (th**2 - (th*thm*css)/sqs)*c2s
					f1hh = (cs*thm - sq*th)*(th*csf + thm*c1) - (thm**2 - (th*thm*cs)/sq)*c2
					f1vv =  - (cs*tnm - sq*tn/egc)*(tn*csf + tnm*egc*c1) + (tnm**2 - (tn*tnm*cs)/sq)*c2
					f1vvs =  - (css*tnm - sqs*tn/egc)*(tn*csf + tnm*egc*c1s) + (tn**2 - (tn*tnm*css)/sqs)*c2s
					f1vhs =  - (css*tm - sqs*tp/egc)*((tp/cs) + tm*egc/sqs)*sf - (tp**2 - (tp*tm*css)/sqs)*(ss**2)*sf
					f1hvs =  - (css*tp - sqs*tm)*(tm/cs + tp/sqs)*sf - (tm**2 - (tp*tm*css)/sqs)*(ss**2)*sf
					f1vh = (cs*tp - sq*tm)*((tm/css) + tp/sq)*sf + ((tp**2) - (tp*tm*cs)/sq)*(s**2)*sf
					f1hv = (cs*tm - sq*tp/egc)*((tp/css) + (tm*egc)/sq)*sf + (tm**2 - tp*tm*cs/sq)*(s**2)*sf
					! ******* scattering coefficient computation ********

					sigma1 = 0.5*(k**2)*EXP( - (sigmaz**2)*(kz*kz + ksz*ksz))

					Sommahh = 0.
					Sommahv = 0.
					Sommavh = 0.
					Sommavv = 0.

					a = (ksz + kz)
					b = EXP( - (sigmaz*sigmaz)*kz*ksz)
					ihh1 = fhh*b
					ihv1 = fhv*b
					ivh1 = fvh*b
					ivv1 = fvv*b
					ihh2 = f1hh
					ihv2 = f1hv
					ivh2 = f1vh
					ivv2 = f1vv
					ihh3 = f1hhs
					ihv3 = f1hvs
					ivh3 = f1vhs
					ivv3 = f1vvs
					nloop: DO n = 1,nsup
						ihh1 = ihh1*a
						ihv1 = ihv1*a
						ivh1 = ivh1*a
						ivv1 = ivv1*a
						ihh2 = ihh2*ksz
						ihv2 = ihv2*ksz
						ivh2 = ivh2*ksz
						ivv2 = ivv2*ksz
						ihh3 = ihh3*kz
						ihv3 = ihv3*kz
						ivh3 = ivh3*kz
						ivv3 = ivv3*kz

						ihhn = ihh1 + 0.5*(ihh2 + ihh3)
						ihvn = ihv1 + 0.5*(ihv2 + ihv3)
						ivhn = ivh1 + 0.5*(ivh2 + ivh3)
						ivvn = ivv1 + 0.5*(ivv2 + ivv3)

						! Roughness spectrum
						IF(CORRELATION  ==  1) THEN
							! 1) Exponential
							wn1 = llcorr/(n*n)
							wn2 = (1. + (llcorr*(kdx*kdx + kdy*kdy))/(n*n))**( - 1.5)
							wn = wn1*wn2
						ELSE
							! 2) Gaussian
							wn = (llcorr/(2*n))*EXP( - llcorr*((kdx*kdx + kdy*kdy)/(4*n)))
						ENDIF

						IF(n  ==  1) THEN
							FFss(n) = sigmaz**2
						ELSE
							FFss(n) = FFss(n - 1)*sigmaz**2/n
						END IF
						!
						Sommahh = Sommahh + wn*CABS(ihhn)*FFss(n)*CABS(ihhn)
						Sommahv = Sommahv + wn*CABS(ihvn)*FFss(n)*CABS(ihvn)
						Sommavh = Sommavh + wn*CABS(ivhn)*FFss(n)*CABS(ivhn)
						Sommavv = Sommavv + wn*CABS(ivvn)*FFss(n)*CABS(ivvn)

					END DO nloop

					sigmahh = (sigma1*Sommahh)
					sigmahv = (sigma1*Sommahv)
					sigmavh = (sigma1*Sommavh)
					sigmavv = (sigma1*Sommavv)

					Gdirez(1,1) = sigmavv
					Gdirez(1,2) = sigmavh
					Gdirez(2,1) = sigmahv
					Gdirez(2,2) = sigmahh
					!
					!  ************************
					!  * Gdirez(1,1) = sigma0vv *
					!  * Gdirez(1,2) = sigma0vh *
					!  * Gdirez(2,1) = sigma0hv *
					!  * Gdirez(2,2) = sigma0hh *
					!  ************************
					!
					Gtoro(k3,:,:) = Gtoro(k3,:,:) + cdirezg*Gdirez(:,:)
					! It is enough to compute the averaged value in fis, either incident radiation and scattered one
					! are uniformly distribuited in the space
					!
				END DO fisloop
			END DO tetasloop
		END DO  tetaloop
		!
		! contributions for 180<fis<360
		Gtoro(NM/2 + 2:NM,:,:) = Gtoro(NM/2:2: - 1,:,:)
		!
		DO istoki = 1,2
			DO istoks = 1,2
				GX(:) = CMPLX(Gtoro(:,istoks,istoki),0.)
				CALL CFSTFT(NF,gx)
				! normalization
				Gx(1) = Gx(1) / NM
				Gx(2:) = Gx(2:)*2/NM
				Sgg(i,j,1:NM1,istoks,istoki) = REAL(GX(1:NM1))
			END DO     ! istoks
		END DO     ! istoki
	END DO toros
END DO toroi
RETURN
END SUBROUTINE IEM
!
!*****************************************************************************
!
SUBROUTINE GEOMETRICALOPTICS(sigmaz,lcorr,vmoist,Sgg,rr,litter_gmoist,biomass)
! Dummy arguments declaration
REAL, INTENT(IN)														  :: sigmaz,vmoist,lcorr
REAL, DIMENSION(:,:,:,:,:),INTENT(OUT)			  :: Sgg
REAL, OPTIONAL, DIMENSION(:,:),INTENT(OUT)  	:: RR
REAL, OPTIONAL, INTENT(IN)                    :: litter_gmoist, biomass
! Local variable declaration
COMPLEX, DIMENSION(NM)		:: gx

REAL, DIMENSION(NN/2,2,2)	:: sbacg0
REAL, DIMENSION(NM,2,2)		:: gtoro
REAL, DIMENSION(2,2)			:: gdirez

COMPLEX		:: egc,eta2,r21,r22,rt,rii,uvv,uvh,uhv,uhh,g2,g3,g4,uvvc,uhhc,uhvc,uvhc
COMPLEX		:: uvvuvvc,uhvuhvc,uvhuvhc,uhhuhhc

REAL			:: lambda,k,kk,nsv,nsh,nivs,nihs,ms,msms,conv,deltet,delfis,sqrt2,c1,ak1,ak1k1,eta1,alfaregr
REAL			:: ross,epsss,g1,betaregr,ff0,egr,egj,alfag,betag,ppg,qqg1,cgg,tetaj,tetais,ctetis,teta
REAL			:: steta,abscteta,tetagraz,ak1steta,qqg,tet2,ctet2,stet2,r11,r12,beckspiz,beckspiz1,beckspiz2
REAL			:: tteta2,expgg,cteta4,tetas,stetas,absctetas,o0,cdirez,cdirezg,fis,sfis,cfis,qx,qqx,qy,qqy,qz,qqz
REAL			:: a0g,tteta,e1g,qg,qgqg,bg,d0q,ctet1,stet1,ak1stet1,erfc,ccg

INTEGER		:: i,j,istoki,istoks,k1,k2,k3

LOGICAL :: RR_present,Litter_present

IF(PRESENT(RR)) THEN
	RR_present = .TRUE.
ELSE
	RR_present = .FALSE.
ENDIF

Litter_present = .FALSE.
IF(PRESENT(litter_gmoist) .AND. PRESENT(biomass)) THEN
  IF(biomass > 0.) Litter_present = .TRUE.
END IF

conv = 180./PI
deltet = PI/(2*NIJ)
delfis = 2.*PI/NM
sqrt2 = SQRT(2.)
c1 = 2./SQRT(PI)
lambda = 30./f
k = 2.*PI/lambda
kk = k*k
ak1 = (2.*PI)/lambda  ! Soil constant
ak1k1 = ak1*ak1
eta1 = 1.

IF(Litter_present) THEN
  CALL LITTER(vmoist,litter_gmoist,biomass,egc) 
ELSE
  CALL DOBSOIL(vmoist,egc) 
END IF

egr = REAL(egc)
egj = AIMAG(egc)
eta2 = CSQRT(1/egc)
alfag = k*AIMAG(CSQRT(egc))
betag = k*REAL(CSQRT(egc))
ppg = 2.*alfag*betag
qqg1 = betag*betag - alfag*alfag
cgg = CABS((eta1 - eta2)/(eta1 + eta2))
cgg = cgg*cgg
!
toroi: DO j = 1,NIJ       ! Incident planes
	tetaj = (j - 1)*deltet

	IF(RR_present) RR(j,:) = 0.
	sbacg0(j,:,:) = 0.

  toros: DO i = 1,NIJ     ! Scattering planes
		tetais = (i - 1)*deltet
		ctetis = COS(tetais)
		gtoro(:,:,:) = 0.
		ms = sqrt2*sigmaz/lcorr
		msms = ms*ms

		tetaloop: DO k1 = 1,2         ! Gauss integral in teta variables
			teta = .5*(deltet*chi(k1) + 2*tetaj + deltet)
			steta = SIN(teta)
			abscteta = COS(teta)
			IF(i == 1) THEN
				tetagraz = PI/2. - teta
				ak1steta = ak1*steta
				qqg = qqg1 - ak1steta*ak1steta
				tet2 = ATAN(ak1steta*sqrt2/SQRT(SQRT(ppg*ppg + qqg*qqg) + qqg))
				ctet2 = COS(tet2)
				stet2 = SIN(tet2)
				r11 = eta1*abscteta
				r12 = eta1*ctet2
				r21 = eta2*abscteta
				r22 = eta2*ctet2
				rt = (r21 - r12)/(r21 + r12)
				rii = (r11 - r22)/(r11 + r22)
				beckspiz1 = 4*PI*sigmaz*SIN(tetagraz)/lambda
				beckspiz2 =  beckspiz1*beckspiz1
				beckspiz = EXP( - beckspiz2)
				rr(j,1) = rr(j,1) + (om(k1)/2.)*beckspiz*rii*CONJG(rii)
				rr(j,2) = rr(j,2) + (om(k1)/2.)*beckspiz*rt*CONJG(rt)
				tteta = steta/abscteta
				tteta2 = tteta*tteta
				expgg = EXP( - tteta2/(2*msms))
				cteta4 = abscteta**4
				sbacg0(j,:,:) = sbacg0(j,:,:) + (om(k1)/2)*cgg*expgg/(2*msms*cteta4)
			ENDIF

			tetasloop: DO k2 = 1,2      ! Gauss integral in tetas variables
				tetas = .5*(deltet*chi(k2) + 2*tetais + deltet)
				stetas = SIN(tetas)
				absctetas = COS(tetas)
				o0 = om(k1)*om(k2)/4.
				cdirez = (deltet*steta/absctetas)*o0 
				cdirezg = cdirez/(4*PI)

				fisloop: DO k3 = 1,NM1    ! Computation for 0<fis<180
					fis = (k3 - 1)*delfis
					sfis = SIN(fis)
					cfis = COS(fis)
					!
					! Soil bistatic coefficient
					qx = ak1*((stetas*cfis) - steta)
					qqx = qx*qx
					qy = ak1*stetas*sfis
					qqy = qy*qy
					qz = ak1*(absctetas + abscteta)
					qqz = qz*qz
					a0g = 2*qqz*msms
					e1g = EXP( - (qqx + qqy)/a0g)
					qg = SQRT(qqx + qqy + qqz)
					qgqg = qg*qg
					bg = a0g*(qqz/qgqg)
					nsh = stetas*sfis
					nsv = abscteta*stetas*cfis + steta*absctetas
					nihs =  - steta*sfis
					nivs = steta*absctetas*cfis + abscteta*stetas
					d0q = nivs*nivs + nihs*nihs
					! Computation of local transmission incident and scattering angles
					! and of Fresnel coefficients
					ctet1 = qg*ABS(qz)/(2.*ak1*qz)  ! Ulaby, appendix 12a
					IF(ctet1 >= 1.) THEN
						ctet1 = 1.
						ctet2 = 1.
					ELSE
						stet1 = SQRT(1 - ctet1*ctet1)
						ak1stet1 = ak1*stet1
						qqg = qqg1 - ak1stet1*ak1stet1
						tet2 = ATAN(ak1stet1*sqrt2/SQRT(SQRT(ppg*ppg + qqg*qqg) + qqg))
						stet2 = SIN(tet2)
						ctet2 = SQRT(1 - stet2*stet2)
					ENDIF
					r11 = eta1*ctet1
					r12 = eta1*ctet2
					r21 = eta2*ctet1
					r22 = eta2*ctet2
					rt = (r21 - r12)/(r21 + r12)
					rii = (r11 - r22)/(r11 + r22)
					uvv = qg*((nsh*nihs*rt) + (nsv*nivs*rii))/(ak1*d0q)
					uvh = qg*( - (nsv*nihs*rt) + (nsh*nivs*rii))/(ak1*d0q)
					uhv = qg*( - (nsh*nivs*rt) + (nsv*nihs*rii))/(ak1*d0q)
					uhh = qg*((nsv*nivs*rt) + (nsh*nihs*rii))/(ak1*d0q)
					uvvc = CONJG(uvv)
					uvhc = CONJG(uvh)
					uhvc = CONJG(uhv)
					uhhc = CONJG(uhh)
					uvvuvvc = uvv*uvvc
					uvhuvhc = uvh*uvhc
					uhvuhvc = uhv*uhvc
					uhhuhhc = uhh*uhhc
					erfc = 1 - ERF(teta,c1,ms)
					! Computation of Stokes (gdirez(i,j)) matrix elements 
					ccg = SH(teta,PI,ms,erfc)*ak1k1*e1g/bg
					gdirez(1,1) = ccg*REAL(uvvuvvc)
					gdirez(1,2) = ccg*REAL(uvhuvhc)
					gdirez(2,1) = ccg*REAL(uhvuhvc)
					gdirez(2,2) = ccg*REAL(uhhuhhc)
					!
					gtoro(k3,:,:) = gtoro(k3,:,:) + cdirezg*gdirez(:,:)
				END DO fisloop
			END DO tetasloop
		END DO tetaloop
		!
		! contributions for 180<fis<360
		gtoro(NM/2 + 2:NM,:,:) = gtoro(NM/2:2: - 1,:,:)
		!
		! Computation of Fourier elements and storage of scattering matrix s
		DO istoki = 1,2
			DO istoks = 1,2
				gx(:) = CMPLX(gtoro(:,istoks,istoki),0.)
				CALL CFSTFT(nf,gx)
				! normalization
				gx(1) = gx(1) / NM
				gx(2:) = gx(2:)*2/NM
				sgg(i,j,1:NM1,istoks,istoki) = REAL(gx(1:NM1))
			END DO     ! istoks
		END DO     ! istoki
	END DO toros
END DO toroi
RETURN
END SUBROUTINE GEOMETRICALOPTICS


REAL FUNCTION ERF(teta,c1,ms)
! Dummy arguments declaration
REAL, INTENT(IN)	:: teta, c1, ms
! Local variable delaration
REAL		x,xx,x1,x2,z,d,sum,y
INTEGER i
IF(teta == 0) THEN 
	erf = 1.
	RETURN
END IF
x1 = (COS(teta)/SIN(teta))
x2 = SQRT(2.)*ms
z = x1/x2
d = z/100
sum = 0
DO i = 1,100    
	x = (2*i - 1)*(z/200)
	xx = x*x
	IF(xx > 180.) THEN
		y = 0.
	ELSE
		y = EXP(-xx)
	ENDIF
	sum = sum + y
END DO
ERF = c1*d*sum

END FUNCTION ERF
!
REAL FUNCTION SH(teta,PI,ms,erfc)
! Dummy arguments declaration
REAL, INTENT(IN)	:: teta, PI, ms, erfc 
! Local variable delaration
REAL a0,a1,aa1,a2,aa3,a3,f0,f1,ff,mm2,mm

IF(teta  ==  0.) THEN
	sh = 1
	RETURN
END IF
a0 = SQRT(2/PI)
a1 = COS(teta)/SIN(teta)
aa1 = a1*a1
a2 = ms/a1
mm = ms*ms
mm2 = mm*2
aa3 = aa1/mm2
IF(aa3 > 180.) THEN
	a3 = 0
ELSE
	a3 = EXP(-(aa1/mm2))
ENDIF
f0 = 0.5*(a0*a2*a3 - erfc)
f1 = 1 + f0
ff = 1/f1
SH = (1 - (0.5*erfc))*ff

END FUNCTION SH
!
!*****************************************************************************
!
!This subroutine contains a subroutine to calculate litter permittivity
!
SUBROUTINE LITTER(vmoist,lmoist,biomass,etc)
                  
REAL, INTENT(IN)      :: vmoist,lmoist
REAL, INTENT(IN)      :: biomass
COMPLEX, INTENT(OUT)  :: etc

INTEGER, PARAMETER	  :: BIOMAX = 201      ,&  ! bio = (0 - 2) [g/cm^2] => 20 [Kg/m^2]
                         NUM_SMC = 201     ,&  ! smc_eq = (0 - 0.4)
                         NUM_K_MAX = 201       ! k = (0 - 200)
!
! The maximum value that the exponential value can assume is limited to 200,
! higher values do not give further improvements to the function trends 
!
REAL, PARAMETER       :: DELTA_BIO = .01   ,&
                         DELTA_K = 1.      ,&
!
! The values of delta k is set to 1, because finer values do not give better
! sensitivity to the function trend
!
                         DELTA_SMC = .002  ,&
                         LITTER_SLOPE = 8.5    ! Slope of line which link the litter biomass [g]
                                               ! to litter thickness [cm]

   
REAL, DIMENSION(NN)            :: RS,RL
REAL, DIMENSION(NN,BIOMAX)     :: r0   
COMPLEX, DIMENSION(NN,BIOMAX)  :: RT0

! Variables used for the interpolation
REAL, DIMENSION(NN,BIOMAX)     :: fexp
REAL, DIMENSION(NUM_K_MAX)     :: k,RMS_K
REAL, DIMENSION(NN)            :: a,b
REAL, DIMENSION(BIOMAX)        :: deltax
INTEGER, DIMENSION(NN,1)       :: minimum
INTEGER                        :: h

! Variables used for the equivalent permittivity
COMPLEX, DIMENSION(NUM_SMC)    :: ec
REAL, DIMENSION(NUM_SMC)       :: RMS_SMC
REAL, DIMENSION(NN)            :: r, Req

! Generic variables
COMPLEX :: elc,k1z,k2z,R1,R2,esp,egc,egceq
REAL    :: elr,elj,beta21,app,dsdwl,Bio_tot,deltet,theta,smceq,Vv
REAL    :: deltah,lambda,sin2,alfa2,ka,L2,costheta2,beta2,omega,epsilon0,k02,k0z,mu0
INTEGER :: j,ja,istoki,ibio,i,ismceq

deltet = PI/NN
lambda = 30./f	
dsdwl = DSDWF

! Soil permittivity
CALL DOBSOIL(vmoist,egc)

! Vegetation permittivity 
CALL VEGEDIEL(dsdwl,lmoist,elr,elj,PERMITTIVITY)
elc = CMPLX(elr, -elj)

! Litter permittivity
Vv = (1 - lmoist*(1 - dsdwl))/ (dsdwl * LITTER_SLOPE)

elc = (1 + Vv*(CSQRT(elc) - 1))**2

! Estimation of litter reflectivity (coherent approach)
Biomassa_tot:  DO ibio = 1,BIOMAX
  Bio_tot = (ibio - 1)*DELTA_BIO    ! [g/cm^2]  
  deltah = Bio_tot * LITTER_SLOPE   ! [cm]

  DO istoki = 1,2  
    DO j = 1,NIJ
      ja = NIJ*(istoki - 1) + j
      theta = (j - .5)*deltet 
      sin2 = (SIN(theta))**2
      costheta2 = CSQRT(1 - (1/elc)*sin2)
      alfa2 = (2*PI/lambda)*ABS(IMAG(CSQRT(elc)))   
      ka = 2*alfa2                                
      L2 = EXP((ka*deltah)/costheta2)              
      beta2 = (2*PI/lambda)*REAL(CSQRT(elc))	 
      beta21 = beta2/costheta2
      app = 2*beta21*deltah
      esp = CMPLX(0., -app)
      omega = 2*PI*f*1.e9
      epsilon0 = 8.854e-12
      mu0 = 4*PI*1.e-7
      k02 = omega**2*epsilon0*mu0
      k0z = SQRT(k02 - k02*sin2)
      k1z = CSQRT(k02*elc - k02*sin2)
      k2z = CSQRT(k02*egc - k02*sin2)
      !
      ! R1 air - litter reflectivity
      ! R2 litter - soil reflectivity
      !
      IF(istoki.EQ.1)THEN
        RS(ja) = (CABS((egc*COS(theta) - CSQRT(egc - sin2))/(egc*COS(theta) + CSQRT(egc - sin2))))**2 
        RL(ja) = (CABS((elc*COS(theta) - CSQRT(elc - sin2))/(elc*COS(theta) + CSQRT(elc - sin2))))**2 
        R1 = (elc*k0z - k1z)/(elc*k0z + k1z)           
        R2 = (egc*k1z - elc*k2z)/(egc*k1z + elc*k2z)	 
      ELSE
        RS(ja) = (CABS((COS(theta) - CSQRT(egc - sin2))/(COS(theta) + CSQRT(egc - sin2))))**2
        RL(ja) = (CABS((COS(theta) - CSQRT(elc - sin2))/(COS(theta) + CSQRT(elc - sin2))))**2
        R1 = (k0z - k1z)/(k0z + k1z)
        R2 = (k1z - k2z)/(k1z + k2z)
      END IF
      RT0(ja,ibio) = (R1 + (R2/L2)*CEXP(esp))/(1 + ((R1*R2*CEXP(esp))/L2))
      r0(ja,ibio) = RT0(ja,ibio)*CONJG(RT0(ja,ibio))
    END DO
  END DO
END DO Biomassa_tot
!
! An exponential function has been used to reproduct the coherent reflectivity function
! f(x) = A e^(-kx) + B
! f(0) = A + B   is the reflectivity coefficient of soil
! f(inf) = B     is the reflectivity coefficient of litter
!
! Interpolation of coherent reflectivity function
DO istoki = 1,2   
  DO j = 1,NIJ
    ja = NIJ*(istoki - 1) + j
    b(ja) = RL(ja)             
    a(ja) = r0(ja,1) - b(ja)
    h = 1
    DO i = 1,NUM_K_MAX
      k(h) = (i - 1)*DELTA_K
      DO ibio = 1,BIOMAX
        Bio_tot = (ibio - 1)*DELTA_BIO                         
        fexp(ja,ibio) = a(ja)*exp(-k(h)*Bio_tot) + b(ja) 
        deltax(ibio) = (r0(ja,ibio) - fexp(ja,ibio))**2
      END DO
      RMS_K(h) = SQRT(SUM(deltax(:))/BIOMAX)
      h = h + 1
    END DO
    minimum(ja,:) = MINLOC(RMS_K(:))
  END DO
END DO
!
! Estimation of equivalent permittivity
i = 1
DO ismceq = 1,NUM_SMC
  smceq = 0.002 + (ismceq - 1)*DELTA_SMC
  CALL DOBSOIL(smceq,egceq)
  ec(i) = egceq
  DO istoki = 1,2   
    DO j = 1,NIJ 
      ja = NIJ*(istoki - 1) + j
      theta = (j - .5)*deltet    
      sin2 = SIN(theta)*SIN(theta)
      ! Reflectivity value for the input biomass value
      r(ja) = a(ja)*exp(-k(minimum(ja,1))*biomass) + b(ja)    
      ! Reflectivity of bare soil with equivalent SMC
      IF(istoki == 1) THEN
        Req(ja) = ((CABS((ec(i)*COS(theta) - CSQRT(ec(i) - sin2))/(ec(i)*COS(theta) + CSQRT(ec(i) - sin2))))**2) 
      ELSE
        Req(ja) = ((CABS((COS(theta) - CSQRT(ec(i) - sin2))/(COS(theta) + CSQRT(ec(i) - sin2))))**2) 
      END IF
    END DO
  END DO   
  RMS_SMC(i) = SQRT(SUM((Req(:) - r(:))**2)/NN) 
  i = i + 1
END DO
minimum(1,:) = MINLOC(RMS_SMC(:))
etc = ec(minimum(1,1))
RETURN

END SUBROUTINE LITTER
!
!******************************************************************************
END MODULE MOD_SOIL_FUNCTIONS
