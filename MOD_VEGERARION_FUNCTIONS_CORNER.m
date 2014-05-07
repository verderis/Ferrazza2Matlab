%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUBROUTINE CORNER
% CALLED BY: MAIN
% Computes scattering coefficients in forward scattering, backscattering and in
% specular direction
%
%SUBROUTINE CORNER(avert,hvert,amoi,dsdw,sigmae,sigmad,sigmau)
function [sigmae,sigmad,sigmau] = MOD_VEGERARION_FUNCTIONS_CORNER(avert,hvert,amoi,dsdw)

%global nm;
%global nn;
global nij;
%global nm1;
global f;
global PERMITTIVITY;
%global chi;
%global om;
global 

% Dummey arguments declaration
%REAL, DIMENSION(:,:,:), INTENT(out)	:: sigmad,sigmau
sigmad = zeros(nij,2,2);
sigmau = zeros(nij,2,2);

%REAL, DIMENSION(:,:), INTENT(out)		:: sigmae
sigmae = zeros(nij,2);

%REAL, INTENT(IN)										:: amoi,avert,hvert,dsdw
% Local variables declaration
%REAL, DIMENSION(7,2,2)							:: ssigmad, ssigmau
ssigmad = zeros(7,2,2); ssigmau = zeros(7,2,2);
%REAL, DIMENSION(7,2)								:: ssigmae
ssigmae = zeros(7,2);

%COMPLEX(8), DIMENSION(0:NMAX0 + 2)	:: ju
ju = zeros(0:NMAX0 + 2,1);
%COMPLEX(8), DIMENSION(0:NMAX0 + 1)	:: z
z = zeros(0:NMAX0 + 1,1);
%COMPLEX(8), DIMENSION(0:NMAX0)			:: hni,hni1,ju1,r,ev,eh,hv,hh,a,b
hni = zeros(0:NMAX0,1);
hni1 = zeros(0:NMAX0,1);
ju1 = zeros(0:NMAX0,1);
r = zeros(0:NMAX0,1);
ev = zeros(0:NMAX0,1);
eh = zeros(0:NMAX0,1);
hv = zeros(0:NMAX0,1);
hh = zeros(0:NMAX0,1);
a = zeros(0:NMAX0,1);
b = zeros(0:NMAX0,1);
%COMPLEX(8)  :: fv1fs,fh1fs,fvfs,fhfs,li,u,uu,uunn,rr1,rr2,rr3,rr4,rr5,zmeno,ff,ff0
%COMPLEX(8)  :: ff1vv0,ff1vh0,ff1hv0,ff1hh0,svv,svh,shv,shh,ff1vv,ff1vh,ff1hv,ff1hh
%COMPLEX(8)  :: fv1,fv2,fh1,fh2,fvv,fvh,fhv,fhh,ec
%REAL(8), DIMENSION(0:NMAX0 + 2) :: js
js = zeros(0:NMAX0 + 2,1);
%REAL(8)     :: lcil,lambda,k,kk,ni,nini,mi,ls,teta,alfa,beta,gamma,salfa,calfa,sbeta,cbeta,sgamma,cgamma
%REAL(8)			:: conv,deltet,delfis,eta,roh2o,ak1,ak1k1,tetagraz,abscteta,steta,tetas,ak1steta,cteta,tvi,thi
%REAL(8)			:: ctetil,cctetil,stetil,tt,cfil,sfil,fil,ure,uim,stetas,ctetas,absctetas,fis,sfis,cfis
%REAL(8)			:: fisg,tvs,ths,ctetsl,cctetsl,stetsl,tts,cfisl,sfisl,fisl,arg2,arg,dd,cks,cka,ccv,cch
%REAL(8)			:: fvvfvv,fhvfhv,fvhfvh,fhhfhh,fvvj,fhhj,alfa1,dalfa,beta1,dbeta,gamma1,dgamma
%REAL, DIMENSION(2)	:: jcos,icos
%REAL				:: er,ej
%INTEGER			:: ialfa,ibeta,igamma,nalfa,nbeta,ngamma,nabg,j,imax,nmax,n,icode,k3
%LOGICAL			:: test
% Variables for library functions
%COMPLEX(8), DIMENSION(0:NMAX0 + 1)	:: JB,YB,dJB,dYB,SIG
JB = zeros(0:NMAX0,1);
YB = zeros(0:NMAX0,1);
dJB = zeros(0:NMAX0,1);
dYB = zeros(0:NMAX0,1);
SIG = zeros(0:NMAX0,1);

%COMPLEX(8)													:: ETAA,ZLMIN,cni

%INTEGER															:: KFN,MODE,JFAIL,JPR

% Initiation of input parameter for WCLBES function 
ETA = 0;
ZLMIN = 0;
KFN = 2;
MODE = 1;

jcos = [-1, -1];
icos = [1, -1];

alfa1 = ALPHA1STE;
nalfa = NALPHASTE;
dalfa = DALPHASTE;
beta1 = BETA1STE;
nbeta = NBETASTE;
dbeta = DBETASTE;
gamma1 = GAMMA1STE;
ngamma = NGAMMASTE;
dgamma = DGAMMASTE;
	
nabg=nalfa*nbeta*ngamma;
conv = 180.d0/pi;
deltet = pi/(2*nij);
delfis = 2.0*pi/nm;
lambda = 30.0/f;
k = 2.0*pi/lambda;
kk = k*k;
eta = 1.0;  
roh2o = 1.0; 
lcil = hvert/2.;

%CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)     % cylinder permittivity
[er, ej] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,PERMITTIVITY);

ec = er -1i*ej;
ff0 = kk*lcil*(ec - 1.0);
ak1 = (2.*pi)/lambda;           
ak1k1 = ak1*ak1;

for j = 1:nij %jloop: DO
	teta = deltet*(j - .5);
	tetagraz = pi/2. - teta;
	abscteta = cos(teta);
	steta = sin(teta);
	tetas = teta;
	ak1steta = ak1*steta;
	cteta = jcos(1)*abscteta;
    test = false;
    imax = 2;
	while (imax <=  7 && test  ==  false) %bessel: DO WHILE
        imax = imax + 1;
	  
    if (imax == 7)
        nmax = 100;
        disp('The Bessel functions cannot be computed with further accuracy, the order is 100 in CORNER')
    else
        nmax = 2^imax;
		end
		ssigmae(imax,:) = 0.0;
		ssigmad(imax,:,:) = 0.0;
		ssigmau(imax,:,:) = 0.0;

		for ialfa = 1:nalfa %alfaloop: DO 
			alfa = alfa1 + (ialfa - 1)*dalfa;
			salfa = sin(alfa);
			calfa = cos(alfa);

			for ibeta = 1:nbeta %betaloop: DO 
				beta = beta1 + (ibeta - 1)*dbeta;
				sbeta = sin(beta);
				cbeta = cos(beta);

				for igamma = 1:ngamma %gammaloop: DO 
					gamma = gamma1 + (igamma - 1)*dgamma;
					sgamma = sin(gamma);
					cgamma = cos(gamma);
					% evaluation of enp,hnp (15), (16).
					% tey depend only from incidence angle
					tvi = (salfa*sgamma + calfa*sbeta*cgamma)*cteta - cbeta*cgamma*steta;
					thi = salfa*sbeta*cgamma - calfa*sgamma;
					ctetil = (salfa*sgamma + calfa*sbeta*cgamma)*steta + cbeta*cgamma*cteta;
					cctetil = ctetil*ctetil;
					stetil = sqrt(1.0 - cctetil);
					tt = sqrt(tvi*tvi + thi*thi);
					if(beta == 0.0 && gamma == 0.0)
						cfil = calfa;
						sfil =  - salfa;
					else
						cfil = (cbeta*calfa*steta - sbeta*cteta)/tt;
						sfil = ((sbeta*sgamma*calfa - cgamma*salfa)*steta + cbeta*sgamma*cteta)/tt;
					end
					fil = acos(cfil);
					if (sfil < 0.0) 
                        fil = 2.0*pi - fil;
                    end
					li = k*sqrt(ec - cctetil);
					u = avert*li;
					ure = real(u);
					uim = imag(u);
					uu = u*u;
					ni = k*avert*stetil;
					nini = ni*ni;
					uunn = 1.0/nini - 1.0/uu;
					% First kind of Bessel functions
					CALL WBSJA(u, 0.d0, nmax, 3, ju)
					% First, Second and Third kinds Bessel functions

					cni = ni;
					CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
					if(JFAIL ~= 0) STOP "Error in WCLBES function in CORNER"
					% Hankel second kind functions computation [0,NMAX0]
					hni(:) = JB(0:NMAX0) - 1i * YB(0:NMAX0)
					% Derivates of Hankel functions 
					hni1(:) = dJB(0:NMAX0) - 1i * dYB(0:NMAX0)
					% Derivates of Bessel functions 
					ju1(0) = -ju(1);                    % first derivatives (9.1.28 of [9])
					for n = 1:nmax
						ju1(n) = ju(n - 1) - (n/u)*ju(n) % first derivatives (9.1.27 of [9])
                    end

					for n = 0:nmax
						rr1 = ju1(n)/(u*ju(n));
						rr2 = hni1(n)/(ni*hni(n));
						rr3 = pi*nini*hni(n)/2.0;
						rr4 = (rr2 - rr1)*(rr2 - ec*rr1);
						rr5 = n*n*cctetil*uunn*uunn;
						r(n) = rr3*(rr4 - rr5);
						ev(n) = 1i*stetil*(rr2 - rr1)/(r(n)*ju(n));
						hv(n) = (stetil/eta)*uunn*n*( -ctetil)/(r(n)*ju(n));
						hh(n) = 1i*(stetil/eta)*(rr2 - ec*rr1)/(r(n)*ju(n));
						eh(n) =  -stetil*uunn*n*( -ctetil)/(r(n)*ju(n));
                    end
					% icode = 1  : upward semispace scattering
					% icode = 2  : downward semispace scattering
					%icodeloop: DO icode = 2,1, - 1
                    for icode = 2:1:-1 %icodeloop: DO 
						%     
						stetas = sin(tetas);
						absctetas = cos(tetas);
						ctetas = icos(icode)*absctetas;

						for k3 = 1,2        %nm1,nm1 - 1      % calcoli per 0<fis<180 %fisloop: DO 
							fis = (k3 - 1)*pi;                 
							sfis = sin(fis);
							cfis = cos(fis);
							fisg = fis*conv;
							tvs = (sin(alfa - fisg)*sgamma + cos(alfa - fisg)*sbeta*cgamma)*ctetas - cbeta*cgamma*stetas;
							ths = sin(alfa - fisg)*sbeta*cgamma - cos(alfa - fisg)*sgamma;
							ctetsl = (sin(alfa - fisg)*sgamma + cos(alfa - fisg)*sbeta*cgamma)*stetas + cbeta*cgamma*ctetas;
							cctetsl = ctetsl*ctetsl;
							stetsl = sqrt(1.0 - cctetsl);
							tts = sqrt(tvs*tvs + ths*ths);

							if(beta == 0.d0 && gamma == 0.d0)
								cfisl = cos(alfa - fisg);
								sfisl =  - sin(alfa - fisg);
							else
								cfisl = (cbeta*stetas*cos(alfa - fisg) - sbeta*ctetas)/tts;
								sfisl = ((sbeta*sgamma*calfa - cgamma*sin(alfa - fisg))*stetas + cbeta*sgamma*ctetas)/tts;
							end
							fisl = acos(cfisl);
							if(sfisl < 0.0)
                                fisl = 2.0*pi - fisl;
                            end

							ls = k*stetsl;                               % (31)
							arg2 = ls*avert;
                            CALL DBSJA(arg2, 0.d0, nmax, 3, js)

							z(0) = (avert/(li*li - ls*ls))*(li*js(0)*ju(1) - ls*ju(0)*js(1));  % (a10)
							z(1) = (avert/(li*li - ls*ls))*(li*js(1)*ju(2) - ls*ju(1)*js(2));  % (a11)
							zmeno = z(1);
							if(ls == 0.0)
								z(2) = z(0) - 2*ju(1)/li*avert*0.5;
							else
								z(2) = z(0) - 2*ju(1)*js(1)/(ls*li);  % (a13)
							end
					
							for n = 3:nmax + 1
								if(ls == 0.0)
									z(n) = z(n - 2);
								else
									z(n) = z(n - 2) - 2*(n - 1)*ju(n - 1)*js(n - 1)/(ls*li);  % (a13)
								end
                            end
							a(0) = (k/(2.0*li))*(zmeno - z(1));  %  (38)
							b(0) = (k/(2.0*li))*(zmeno + z(1));
							%                                  calcolo (42)
							ff1vv0 = ev(0)*(b(0)*ctetsl*( - ctetil) - stetsl*z(0));    
							ff1vh0 = 0.0;
							ff1hv0 = 0.0;
							ff1hh0 = b(0)*eta*hh(0);
							for n = 1:nmax
								a(n) = (k/(2*li))*(z(n - 1) - z(n + 1));      % (38)
								b(n) = (k/(2*li))*(z(n - 1) + z(n + 1));      % (38)
								svv = 2.d0*((ev(n)*( - ctetil)*b(n) - 1i*eta*hv(n)*a(n))*ctetsl - stetsl*ev(n)*z(n)); 
								svh = ((eh(n)*( - ctetil)*b(n) - 1i*eta*hh(n)*a(n))*ctetsl - stetsl*eh(n)*z(n));
								shv = (eta*hv(n)*b(n) + 1i*ev(n)*( - ctetil)*a(n));
								shh = 2.d0*(eta*hh(n)*b(n) + 1i*eh(n)*( - ctetil)*a(n));
								ff1vv0 = ff1vv0 + svv*cos(n*(fisl - fil));
								ff1vh0 = ff1vh0 + svh*sin(n*(fisl - fil));
								ff1hv0 = ff1hv0 + shv*sin(n*(fisl - fil));
								ff1hh0 = ff1hh0 + shh*cos(n*(fisl - fil));
                            end % n
							arg = k*lcil*(ctetsl - ctetil);              % (28)

							if(arg == 0.0)
								mi = 1.0;
							else
								mi = sin(arg)/arg;
							end
							ff = ff0*mi;
							ff1vv = ff1vv0*ff;
							ff1vh = 2.0*1i*ff1vh0*ff;
							ff1hv = 2.d0*1i*ff1hv0*ff;
							ff1hh = ff1hh0*ff;
							%
							fv1 =  -ff1vv*tvs + ff1hv*ths;
							fv2 =  -ff1vh*tvs + ff1hh*ths;
							fh1 =  -ff1vv*ths - ff1hv*tvs;
							fh2 =  -ff1vh*ths - ff1hh*tvs;
							dd = tt*tts;
							if(beta == 0.d0 && gamma == 0.d0)
								fvv = ff1vv;
								fvh = ff1vh;
								fhv = ff1hv,
								fhh = ff1hh;
							else
								fvv = ( -fv1*tvi + fv2*thi)/dd;
								fvh = ( -fv2*tvi - fv1*thi)/dd;
								fhv = ( -fh1*tvi + fh2*thi)/dd;
								fhh = ( -fh1*thi - fh2*tvi)/dd;
							end
							%
							fvvfvv = fvv*fvv';
							fvhfvh = fvh*fvh';
							fhvfhv = fhv*fhv';
							fhhfhh = fhh*fhh';
							if (k3  ==  2) 
								if (icode  ==  2) 
									ssigmad(imax,1,1) = ssigmad(imax,1,1) + 4*pi*fvvfvv/nabg;
									ssigmad(imax,2,1) = ssigmad(imax,2,1) + 4*pi*fhvfhv/nabg;
									ssigmad(imax,1,2) = ssigmad(imax,1,2) + 4*pi*fvhfvh/nabg;
									ssigmad(imax,2,2) = ssigmad(imax,2,2) + 4*pi*fhhfhh/nabg;
								else
									cks = 4*pi/(nabg*abscteta);
									ssigmau(imax,1,1) = ssigmau(imax,1,1) + cks*fvvfvv;
									ssigmau(imax,2,1) = ssigmau(imax,2,1) + cks*fhvfhv;
									ssigmau(imax,1,2) = ssigmau(imax,1,2) + cks*fvhfvh;
									ssigmau(imax,2,2) = ssigmau(imax,2,2) + cks*fhhfhh;
                                end
                            end

							if(k3 == 1 && icode == 2)
								fv1fs = real(ff1vv)+1i*abs(imag(ff1vv));
								fh1fs = real(ff1hh)+1i*abs(imag(ff1hh));
								if(beta == 0. && gamma == 0.)
									fvfs = fv1fs;
									fhfs = fh1fs;
								else
									fvfs = (fv1fs*tvs*tvi + fh1fs*ths*thi)/dd;          % (49)
									fhfs = (fv1fs*ths*thi + fh1fs*tvs*tvi)/dd;          % (49)
								end
								%
								fvvj = imag(fvfs);
								fhhj = imag(fhfs);
								cka = (4*pi/k)/(abscteta*nabg);
								ssigmae(imax,1) = ssigmae(imax,1) + cka*abs(fvvj);
								ssigmae(imax,2) = ssigmae(imax,2) + cka*abs(fhhj)
                            end % if
                        end %END DO fisloop
                    end % DO icodeloop
                    end %ENDDO gammaloop
                end % DO betaloop
            end %ENDDO alfaloop
		%
		%
		if(imax > 3) 
			ccv = (ssigmae(imax,1) - ssigmae(imax - 1,1))/ssigmae(imax,1);
			cch = (ssigmae(imax,2) - ssigmae(imax - 1,2))/ssigmae(imax,2);
			if(ccv < .01 && cch < .01) test = true
            end %END if
        end %END DO bessel

	sigmae(j,:) = ssigmae(imax,:);
	sigmad(j,:,:) = ssigmad(imax,:,:);
	sigmau(j,:,:) = ssigmau(imax,:,:);

        end %END DO jloop

%RETURN

%END SUBROUTINE CORNER