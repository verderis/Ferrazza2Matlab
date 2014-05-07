% SUBROUTINE TRUNKS_ABS
% CALLED BY: MAIN
% Computes the absorption coefficient ensemble of "vertical" cylinders (main stems) over the soil.
% Requires double precision.
%
%SUBROUTINE TRUNKS_ABS(acil,length,amoist,dsdw,Ka,Kes)
function [Ka,Kes] = MOD_VEGERARION_FUNCTIONS_TRUNKS_ABS(acil,length,amoist,dsdw)

%Dummy arguments declaration
REAL,INTENT(IN)										        :: length,acil,amoist,dsdw
REAL,DIMENSION(:,:),INTENT(OUT)		        :: Ka
REAL,DIMENSION(:,:),INTENT(OUT), OPTIONAL	:: Kes

%Local variables declaration
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

% Variables for library functions
COMPLEX(8), DIMENSION(-NMAX0-1 : NMAX0 + 1)	:: JU,JB,YB,dJU,dJB,dYB,Hn,dHn,SIG
COMPLEX(8)																	:: ETAA,ZLMIN,cni

INTEGER																			:: KFN,MODE,JFAIL,JPR

% Initiation of input parameter for WCLBES function 
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
eta=1.d0   % e' ininfluente

CALL VEGEDIEL(dsdw,amoist,er,ej,PERMITTIVITY)
ec = CMPLX(er, -ej, 8)

lcil = length/2
ccabs=4.*k*DPI*lcil*ej

toroloop: DO j=1,nij       % toro di incidenza
	tetaj=(j-1)*deltet
	ke(j,:)=0.d0
	ka(j,:)=0.d0
	kav=0.; kav1=0.; kav2=0.; kav3=0.
	kah=0.; kah1=0.; kah2=0.; kah3=0.
	%<--------------------------------------------------loop di icode

	icodeLoop: DO icode=2,2
		% icode=1  : scattering semispazio superiore
		% icode=2  : scattering semispazio inferiore
		%     
		%<--------------------------------------------------loop di i
		xtoro(3:,:,:)=0.d0
		iflag=0               !indice controllo nmax
		imax=3
		% --------------------------------- loop alfa,beta,gamma

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

					tetaLoop: DO k1=1,2     % integrale (gauss) in teta
						teta=.5d0*(deltet*chi(k1)+2.d0*tetaj+deltet)
						!teta=tetaj
						steta=DSIN(teta)
						abscteta=dCOS(teta)
						cteta=jCOS(icode)*abscteta
						% calcolo enp,hnp (15),(16). Dipendono solo dall'angolo di incidenza; (fi=0).
						tvi=(salfa*sgamma+calfa*sbeta*cgamma)*cteta-cbeta*cgamma*steta      % (7)
						thi=salfa*sbeta*cgamma-calfa*sgamma                                 % (7)
						ctetil=(salfa*sgamma+calfa*sbeta*cgamma)*steta+cbeta*cgamma*cteta   % (6)
						cctetil=ctetil*ctetil
						stetil=0.d0
						IF(cctetil.le.1.d0) stetil=DSQRT(1.d0-cctetil)
						tt=DSQRT(tvi*tvi+thi*thi)
						% steta e' sempre diverso da 0!%
						IF(beta.eq.0.d0.and.gamma.eq.0.d0)THEN
							cfil=calfa
							sfil=-salfa
						ELSE
							cfil=(cbeta*calfa*steta-sbeta*cteta)/tt             % (6)
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
						% Computation of negative orders [-nmax,-1] rif. [2]

          	% First and second types of Bessel functions are computed for positive order
						cni = (1.d0, 0.d0) * ni
						CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
						IF(JFAIL /= 0) STOP "Error in WCLBES function in TRUNK_ABS "
						% Computation of negative orders [-nmax,-1] rif. [2]
						DO n = nmax,1,-1
							JU(-n) = JU(n)*(-1.d0,0.d0)**n
							JB(-n) = JB(n)*(-1.d0,0.d0)**n
							YB(-n) = YB(n)*(-1.d0,0.d0)**n
							dJU(-n) = dJU(n)*(-1.d0,0.d0)**n
							dJB(-n) = dJB(n)*(-1.d0,0.d0)**n
							dYB(-n) = dYB(n)*(-1.d0,0.d0)**n
						END DO
						% Hankel functions computation [-nmax,nmax]
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
						%
						% cdirez e' il coefficiente per il calcolo del contributo di una coppia di direzioni alla matrice 
						% di scatter di una coppia di direzioni alla matrice di scatter di una coppia di tori (vd. (2) di appunto).                   
						% L'integrale e' calcolato con la quadratura di gauss, il fattore deltet**2 si elide perche' sono 
						% valori medi. Cdirezg e' per la sigma del terreno (v. appunto emissivita')
						%                                                                 
						fisLoop: DO k3=1,1      % calcoli per fis=0                    
							fis=0.
							sfis=SIN(fis)     
							cfis=COS(fis)     
							fisg=fis*conv
							% espressioni analoghe alle (6) e (7) per direzione di scattering.
							tvs=(DSIND(alfa-fisg)*sgamma+DCOSD(alfa-fisg)*sbeta*cgamma)*ctetas-cbeta*cgamma*stetas
							ths=DSIND(alfa-fisg)*sbeta*cgamma-DCOSD(alfa-fisg)*sgamma
							ctetsl=(DSIND(alfa-fisg)*sgamma+DCOSD(alfa-fisg)*sbeta*cgamma)*stetas+cbeta*cgamma*ctetas
							cctetsl=ctetsl*ctetsl
							stetsl=0.d0
							IF(cctetsl.le.1.d0)  stetsl=DSQRT(1.d0-cctetsl)
							tts=DSQRT(tvs*tvs+ths*ths)
							argmi=k*lcil*(ctetsl-ctetil)              % (28)

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
							ls=k*stetsl                               % (31)
							arg2=ls*acil

							% First kind of Bessel functions
							CALL DBSJA(arg2, 0.d0, nmax, 3, js)

							z(0)=(acil/(li*li-ls*ls))*(li*js(0)*JU(1)-ls*JU(0)*js(1))  % (a10)
							z(1)=(acil/(li*li-ls*ls))*(li*js(1)*JU(2)-ls*JU(1)*js(2))  % (a11)
							!zmen=(acil/(li*li-ls*ls))*(-li*js(1)*JU(0)+ls*JU(1)*js(0))  % (a11)+abram.,9.1.5
							zmen=z(1)
							IF(ls.eq.0.d0)THEN
								z(2)=z(0)-2*JU(1)/li*acil*0.5
							ELSE
								z(2)=z(0)-2*JU(1)*js(1)/(ls*li)  % (a13)
							ENDIF
							DO n=3,nmax+1
								IF(ls.eq.0.d0)THEN
									z(n)=z(n-2)
								ELSE
									z(n)=z(n-2)-2*(n-1)*JU(n-1)*js(n-1)/(ls*li)  % (a13)
								ENDIF
							END DO 
							a(0)=(k/(2.d0*li))*(zmen-z(1))  %  (38)
							b(0)=(k/(2.d0*li))*(zmen+z(1))
							%                                  calcolo (42)
							ff=kk*lcil*mi*(ec-1.d0)               
							ff1vv=ev(0)*(b(0)*ctetsl*(-ctetil)-stetsl*z(0))
							ff1vh=0.d0
							ff1hv=0.d0
							ff1hh=b(0)*eta*hh(0)
							DO n=1,nmax
								a(n)=(k/(2*li))*(z(n-1)-z(n+1))      % (38)
								b(n)=(k/(2*li))*(z(n-1)+z(n+1))      % (38)
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
								%  Cambiamento di riferimento
								fv1=-ff1vv*tvs+ff1hv*ths       % (49)
								fv2=-ff1vh*tvs+ff1hh*ths
								fh1=-ff1vv*ths-ff1hv*tvs
								fh2=-ff1vh*ths-ff1hh*tvs
								dd=tt*tts    
								fvv=(-fv1*tvi+fv2*thi)/dd      % (48)
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