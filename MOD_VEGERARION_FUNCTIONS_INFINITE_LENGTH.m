%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUBROUTINE INFINITE_LENGTH
% CALLED BY: MAIN
% Computes the scatter matrix and the extinction vector using the Infinite-length approximation for an ensemble of cylinders.
% Requires double precision.
%
%SUBROUTINE INFINITE_LENGTH(acyl,lcyl,amoi,dsdw,s,ke,select)

function [s,ke] = MOD_VEGERARION_FUNCTIONS_INFINITE_LENGTH(acyl,lcyl,amoi,dsdw,select)

global nm;
global nn;
global nij;
global nm1;
global f;
global PERMITTIVITY;
global chi;
global om;
global nmax0;
global ns;
global nf;

% Leaves, used with Raylegh-Gans and Pysical Optics approximation
 global ALPHA1DIS; 
 global NALPHADIS; 
 global DALPHADIS; 
 global BETA1DIS; 
 global NBETADIS; 
 global DBETADIS; 
 global GAMMA1DIS; 
 global NGAMMADIS; 
 global DGAMMADIS; 
% Petioles, needles and secondary branches used with Rayleigh-Gans, Infinte length
% and Pysical Optics approximation
 global ALPHA1PET;
 global NALPHAPET;
 global DALPHAPET; 
 global BETA1PET; 
 global NBETAPET;
 global DBETAPET;
 global GAMMA1PET; 
 global NGAMMAPET;
 global DGAMMAPET; 
% Ears used with Infinite length approximation
 global ALPHA1STE; 
 global NALPHASTE; 
 global DALPHASTE; 
 global BETA1STE; 
 global NBETASTE; 
 global DBETASTE; 
 global GAMMA1STE; 
 global NGAMMASTE; 
 global DGAMMASTE; 

% Dummy variables declaration
%real, DIMENSION(:,:,:,:,:,:), INTENT(OUT)   :: s
s = zeros(nij,nij,2,nm1,2,2); %warning: no estoy seguro de nm1!
%real, DIMENSION(:,:), INTENT(OUT)		        :: ke
ke = zeros(nn/2,2);

% real, INTENT(IN)						                :: acyl,lcyl,amoi,dsdw
% CHARACTER, INTENT(IN)					              :: select
% Lopcal variables declaration

%COMPLEX(8), DIMENSION(0:nmax0+1)	:: z
z = zeros(nmax0+1+1,1);
%COMPLEX(8), DIMENSION(0:nmax0)		:: hni,hni1,ju,ju1,r,ev,eh,hv,hh,a,b
hni = zeros(nmax0+1+1,1); hni1 = zeros(nmax0+1+1,1); ju = zeros(nmax0+1+1,1); ju1 = zeros(nmax0+1+1,1);
r = zeros(nmax0+1+1,1); ev = zeros(nmax0+1+1,1); hv = zeros(nmax0+1+1,1); hh = zeros(nmax0+1+1,1);
a = zeros(nmax0+1+1,1); b = zeros(nmax0+1+1,1); eh = zeros(nmax0+1+1,1);

%COMPLEX(8)							          :: li,u,uu,uunn,rr1,rr2,rr3,rr4,rr5,zmin,ff,ec,svv,svh,shv,shh
%COMPLEX(8)							          :: fv1,fv2,fh1,fh2,fvv,fvh,fhv,fhh,fvfs,fhfs,ff1vv,ff1vh,ff1hv,ff1hh,fv1fs,fh1fs

%COMPLEX, DIMENSION(NM)				    :: xx

%real(8), DIMENSION(NM,7,2,2)		  :: xff	
xff = zeros(nm,7, 2,2);

%real(8), DIMENSION(0:nmax0)				:: js
js = zeros(nmax0+1+1,1);

%real, DIMENSION(2)				        :: icos,jcos

%real(8)								:: lambda,k,kk,ni,nini,mi,ls,conv,delthet,delphis,eta,theta,thetaj,thetais,alpha,beta,gamma
%real(8)								:: salpha,calpha,sbeta,cbeta,sgamma,cgamma,stheta,ctheta,absctheta,tvi,thi,cthetil,ccthetil
%real(8)								:: sthetil,sphil,tt,cphil,phil,ure,uim,thetas,sthetas,cthetas,abscthetas,o0,cdirez
%real(8)								:: phis,sphis,cphis,ths,cthetsl,ccthetsl,sthetsl,tts,argmi,phisg,tvs,cphisl,sphisl
%real(8)								:: phisl,arg2,dd,fvvfvv,fhhfhh,fhvfhv,fvhfvh,cke,rap1,rap2

%real							  	:: lcyl2,er,ej,alpha1,dalpha,beta1,dbeta,gamma1,dgamma,aww,sww,ww

%INTEGER								:: nabg,i,j,k1,k2,k3,k3a,iflag,imax,ialpha,ibeta,igamma,icode,nmax,n,istoki,istoks
%INTEGER								:: nalpha,nbeta,ngamma,icosn

% Variables for library functions
%%COMPLEX(8), DIMENSION(0:nmax0 + 1)	:: JB,YB,dJB,dYB,SIG
JB = zeros(nmax0 + 2,1); YB = zeros(nmax0 + 2,1); dJB = zeros(nmax0 + 2,1); dYB = zeros(nmax0 + 2,1); SIG = zeros(nmax0 + 2,1);
jb = zeros(nmax0+1+1,1); yb = zeros(nmax0+1+1,1);  dyb = zeros(nmax0+1+1,1);
djb = zeros(nmax0+1+1,1); sig = zeros(nmax0+1+1,1);

%COMPLEX(8)													:: ETAA,ZLMIN,cni

%INTEGER															:: KFN,MODE,JFAIL,JPR
jcos = [-1, -1];
icos = [1, -1];
toDebug = 0;

% Initiation of input parameter for WCLBES function 
%ETAA = 0;
%ZLMIN = 0;
%KFN = 2;
%MODE = 1;

icosn = 1;
aww = 1.;

%SELECT CASE (select)
switch select
    %CASE ('S')
	case 's'
    % Stem or Trunks selection
    alpha1 = ALPHA1STE;
    nalpha = NALPHASTE;
    dalpha = DALPHASTE;
    beta1 = BETA1STE;
    nbeta = NBETASTE;
    dbeta = DBETASTE;
    gamma1 = GAMMA1STE;
    ngamma = NGAMMASTE;
    dgamma = DGAMMASTE;
%CASE('P','N','R')
    case 'P''N''R'
    % Petioles, needle, secondary branches
    alpha1 = ALPHA1PET;
    nalpha = NALPHAPET;
    dalpha = DALPHAPET;
    beta1 = BETA1PET;
    nbeta = NBETAPET;
    dbeta = DBETAPET;
    gamma1 = GAMMA1PET;
    ngamma = NGAMMAPET;
    dgamma = DGAMMAPET;
%CASE('E')
    case 'E'
    % Ears
    alpha1 = ALPHA1EAR;
    nalpha = NALPHAEAR;
    dalpha = DALPHAEAR;
    beta1 = BETA1EAR;
    nbeta = NBETAEAR;
    dbeta = DBETAEAR;
    gamma1 = GAMMA1EAR;
    ngamma = NGAMMAEAR;
    dgamma = DGAMMAEAR;
%CASE('B')
    case 'B'
    % Branches selection, in this case a weight function is used
    alpha1 = ALPHA1BRA;
    nalpha = NALPHABRA;
    dalpha = DALPHABRA;
    beta1 = BETA1BRA;
    nbeta = NBETABRA;
    dbeta = DBETABRA;
    gamma1 = GAMMA1BRA;
    ngamma = NGAMMABRA;
    dgamma = DGAMMABRA;
    icosn=2;
    sww=0.;
    for ibeta = 1:nbeta
        beta = beta1 + (ibeta - 1)*dbeta;
        sww = sww + cos(beta - 90)^2;
    end
    aww = nbeta/sww;
end

nabg = nalpha*nbeta*ngamma;
conv = 180/pi;
delthet = pi/(2*nij);	    	% interval amplitude in theta and theta_s
delphis = 2*pi/nm;	    	% interval amplitude in phi_s - phi
%
lambda = 30.d0/f;				    % wavelength, in cm
k = 2.0*pi/lambda;			    % wavenumber, in cm^ - 1
kk = k*k;
eta = 1.0;					        % free-space impedance
lcyl2 = lcyl/2.;
%
%CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)
   	                        % permittivity computation
[er, ej] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,PERMITTIVITY);
ec = er -1i*ej;

for j = 1:nij %jLoop: DO			  	% j = incidence angular interval index
	thetaj = (j - 1)*delthet;	% lower limit of angular interval in theta
    ke(j,:) = 0.;				    	% initialization to compute extinction coefficients   
	for icode = 1:2 %icodeLoop: DO 
    % icode = 1  : upper half - space scattering
    % icode = 2  : lower half - space scattering
        for i = 1:nij %iLoop: DO          % i = scattering angular interval index
            thetais = (i - 1)*delthet;		% lower limit of angular interval in theta_s
            xff(:,:,:,:) = 0;						      % scattering function initialization
            iflag = 0;						        % Controls order of Bessel functions series
            imax = 1;
            for ialpha = 1:nalpha %alphaLoop: DO 
                alpha = alpha1 + (ialpha - 1)*dalpha;
                salpha = sind(alpha);
                calpha = cosd(alpha);
                for ibeta = 1: nbeta %betaLoop: DO 
                    
                    beta = beta1 + (ibeta - 1)*dbeta;
                    sbeta = sind(beta);
                    cbeta = cosd(beta);
                    if (icosn == 2) 
                        ww = aww*cosd(beta - 90)^2;
                    else
                        ww = 1.;
                    end
                    for igamma = 1:ngamma %gammaLoop: DO 
                        
                        gamma = gamma1 + (igamma - 1)*dgamma;
                        sgamma = sind(gamma);
                        cgamma = cosd(gamma);
                        for k1 = 1:ns %k1Loop: DO      % integration within jth theta interval
                            % using Gauss - Legendre technique
                            
                            %iflag = 0; %de este no estoy seguro!
                            
                            theta = thetaj + .50*delthet*(1. + chi(k1));
                            stheta = sin(theta);
                            absctheta = cos(theta);
                            ctheta = jcos(icode)*absctheta;
                            % computation of enp,hnp (15) of [7].(dependent
                            % on incidence angle only)
                            tvi = (salpha*sgamma + calpha*sbeta*cgamma)*ctheta - cbeta*cgamma*stheta;                    % (7) of [7]
                            thi = salpha*sbeta*cgamma - calpha*sgamma;													% (7) of [7]
                            cthetil = (salpha*sgamma + calpha*sbeta*cgamma)*stheta + cbeta*cgamma*ctheta;                % (6) of [7]
                            ccthetil = cthetil*cthetil;
                            sthetil = 0.d0;
                            if(ccthetil <= 1.d0) 
                            	sthetil = sqrt(1.d0 - ccthetil);
                            end
                            tt = sqrt(tvi*tvi + thi*thi);
                            if (beta == 0.0 && gamma == 0.0)
                                cphil = calpha;
                                sphil =  - salpha;
                            else
                                cphil = (cbeta*calpha*stheta - sbeta*ctheta)/tt;    % (6) of [7]
                                sphil = ((sbeta*sgamma*calpha - cgamma*salpha)*stheta + cbeta*sgamma*ctheta)/tt;
                            end
                            phil = 0.d0;
                            if (cphil < 0.d0) 
                                phil = pi;
                            end
                            if (abs(cphil) <= 1.d0) 
                                phil = acos(cphil);
                            end
                            if(sphil < 0.d0) 
                                phil = 2.d0*pi - phil;
                            end
                            li = k*sqrt(ec - ccthetil);                      % (15) of [7]
                            u = acyl*li;                                       % (16) of [7]
                            ure = real(u);
                            uim = imag(u);
                            uu = u*u;
                            ni = k*acyl*sthetil;                               % (16) of [7]
                            nini = ni*ni;
                            uunn = 1.0/nini - 1.0/uu;
                            
                            % The number of Bessel functions is increased until  relative error is  <  1/1000 (see control later)
                            %% Error aca! el iflag se resetea recien
                            %% arriba, y no entra nunca mas al loop
                            condition = 1;
                            while condition 
                                %while iflag == 0
                                if (iflag == 0)
                                    imax = imax + 1;
                                    nmax = 2^imax;
                                end
                                if (imax == 7)
                                    nmax = 100;
                                    iflag = 1;
                                    disp('The Bessel functions cannot be computed with further accuracy, the order is 100 in INFINITE_LENGTH')
                                end
                                % First kind of Bessel functions
                                %%%%%%%%%%%%%%
                                %Z (COMPLEX) Argument z.
                                %A (REAL) Order a of the first Bessel function in the computed sequence.
                                %NL (INTEGER) Specifies the order a+NL of the last Bessel function in the computed sequence.
                                %ND (INTEGER) Requested number of correct significant decimal digits.
                                %CB (COMPLEX) One-dimensional array with dimension (0:d) where  tex2html_wrap_inline157. On exit, CB(n),  tex2html_wrap_inline159 , contains  tex2html_wrap_inline161 .
                                %CALL WBSJA(u, 0.d0, nmax, 3, ju)
                                %%%%%%%%%%%%%%                                  
                                % First, Second and Third kinds Bessel functions
                                %CALL WCLBES(cni,ETAA,ZLMIN,nmax,JB(0:nmax+1),YB(0:nmax+1),dJB(0:nmax+1),dYB(0:nmax+1),SIG,KFN,MODE,JFAIL,JPR)
                                %CALL WCLBES(Z,ETA,ZLMIN,NL,F,G,FP,GP,SIG,KFN,NODE,JFAIL,JPR)
                                % http://cmd.inp.nsk.su/old/cmd2/manuals/cernlib/shortwrups/node48.html
                                % Initiation of input parameter for WCLBES function 
                                %ETAA = 0;
                                %ZLMIN = 0;
                                %KFN = 2;
                                %MODE = 1;
                                %for ii = 1:nmax
                                for ii = 0:nmax
                                    cni = ni;
                                    ju(ii+1) = besselj(ii,u); %REVISAR!                                            
                                    JB(ii+1) = besselj(ii,cni); %first kind! 
                                    YB(ii+1) = bessely(ii,cni); %second kind!
                                    H(ii+1) = besselh(ii,cni); %third kind!
                                end
                                %if(JFAIL /= 0) STOP "Error in WCLBES function in INFINITE_LENGTH"
                                % Hankel second kind functions computation [-NH,NH]
                                %hni(0:) = JB(0:nmax0) - (0.d0, 1.d0) * YB(0:nmax0);

                                %http://mathworld.wolfram.com/HankelFunctionoftheSecondKind.html
                                hni = JB - 1i*YB;   
                                % Derivates of Hankel functions 
                                %syms z;
                                %syms nu;
                                %diff(besselj(nu,z)-1i*bessely(nu,z))
                                %hni1(:) = dJB(0:nmax0) - (0.d0, 1.d0) * dYB(0:nmax0);
                                for ii = 0:nmax                                       
                                	hni1(ii+1) = (ii*besselj(ii, cni))/cni + bessely(ii + 1, cni)*1i - besselj(ii + 1, cni) - (ii*bessely(ii, cni)*1i)/cni; % Derivates of Hankel functions
                                    % Derivates of Bessel functions
                                    %ju1(0) = -ju(1);                         % first derivatives (9.1.28 of [9])
                                    %
                                    %for n = 1,nmax
                                    %  ju1(n) = ju(n - 1) - (n/u)*ju(n);      % first derivatives (9.1.27 of [9])
                                    %end
                                    %diff(besselj(nu,z))
                                    % ju1(ii+1) = (ii*besselj(ii,
                                    % u))/cni - besselj(ii + 1, u);
                                    % %Fran: este no anda, vuelvo a lo de Paolo

                                end
                                ju1(1) = -ju(2); 
                                for n = 2:nmax+1
                                    ju1(n) = ju(n - 1) - ((n-1)/u)*ju(n);     %first derivatives (9.1.27 of [9])
                                end
                                for n = 1:nmax+1
                                    rr1 = ju1(n)/(u*ju(n));
                                    rr2 = hni1(n)/(ni*hni(n));
                                    rr3 = pi*nini*hni(n)/2.0;
                                    rr4 = (rr2 - rr1)*(rr2 - ec*rr1);
                                    rr5 = (n-1)*(n-1)*ccthetil*uunn*uunn;
                                    r(n) = rr3*(rr4 - rr5);
                                    ev(n) = 1i*sthetil*(rr2 - rr1)/(r(n)*ju(n));           % (15) of [7]
                                    hv(n) = (sthetil/eta)*uunn*(n-1)*( - cthetil)/(r(n)*ju(n));         % (15) of [7]
                                    hh(n) = 1i*(sthetil/eta)*(rr2 - ec*rr1)/(r(n)*ju(n));  % (23) of [7]
                                    eh(n) =  - sthetil*uunn*(n-1)*( - cthetil)/(r(n)*ju(n));            % (23) of [7]
                                end
                                for k2 = 1:ns %k2Loop: DO      % integration within ith theta_s interval using Gauss - Legendre technique
                                    thetas = thetais + .50*delthet*(1. + chi(k2));
                                    sthetas = sin(thetas);
                                    abscthetas = cos(thetas);
                                    cthetas = icos(icode)*abscthetas;
                                    o0 = om(k1)*om(k2)/4.0;  % to be used in the double integration in theta and theta_s
                                    % The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor,
                                    % to be used in formulas (3) and (4) of [7].
                                    cdirez = (delthet*stheta/abscthetas)*o0;
                                    % Computations for 0 < phi_s - phi < 180
                                    for k3 = 1:nm1 %k3Loop: DO 
                                        if (0)
                                            disp(strcat('j=',num2str(j)))
                                            disp(strcat('i=',num2str(i)))
                                            disp(strcat('ialpha=',num2str(ialpha)))
                                            disp(strcat('k1=',num2str(k1)))
                                            disp(strcat('k2=',num2str(k2)))
                                            disp(strcat('k3=',num2str(k3)))
                                            disp(strcat('iflag=',num2str(iflag)))
                                            disp(strcat('condition=',num2str(condition)))
                                            disp('-----------------------------')
                                        end
                                        
                                        phis = (k3 - 1)*delphis;    % phi_s - phi
                                        sphis = sin(phis);
                                        cphis = cos(phis);
                                        phisg = phis*conv;
                                        % (6) and (7) of [7] in the scattering direction
                                        tvs = (sind(alpha - phisg)*sgamma + cosd(alpha - phisg)*sbeta*cgamma)*cthetas - cbeta*cgamma*sthetas;
                                        ths = sind(alpha - phisg)*sbeta*cgamma - cosd(alpha - phisg)*sgamma;
                                        cthetsl = (sind(alpha - phisg)*sgamma + cosd(alpha - phisg)*sbeta*cgamma)*sthetas + cbeta*cgamma*cthetas;
                                        ccthetsl = cthetsl*cthetsl;
                                        sthetsl = 0.d0;
                                        if (ccthetsl <= 1.d0) 
                                            sthetsl = sqrt(1.d0 - ccthetsl);
                                        end
                                        tts = sqrt(tvs*tvs + ths*ths);
                                        argmi = k*lcyl2*(cthetsl - cthetil);              % (28 of [7])
                                        if (argmi == 0.0)
                                            mi = 1.d0;
                                        else
                                            mi = sin(argmi)/argmi;
                                        end
                                        if (beta == 0.d0 && gamma == 0.d0) 
                                            cphisl = cosd(alpha - phisg);
                                            sphisl =  - sind(alpha - phisg);
                                        else
                                            cphisl = (cbeta*sthetas*cosd(alpha - phisg) - sbeta*cthetas)/tts;
                                            sphisl = ((sbeta*sgamma*cosd(alpha - phisg) - cgamma*sin(alpha - phisg))*sthetas + cbeta*sgamma*cthetas)/tts;
                                        end
                                        phisl = 0.0;
                                        if (cphisl < 0.d0) 
                                            phisl = pi;
                                        end
                                        if (abs(cphisl) <= 1.d0) 
                                            phisl = acos(cphisl);
                                        end
                                        if (sphisl  <  0.d0) 
                                            phisl = 2.d0*pi - phisl;
                                        end
                                        ls = k*sthetsl;                        % (31 of [7])
                                        arg2 = ls*acyl;
                                        %CALL DBSJA(arg2, 0.d0, nmax, 3, js)
                                        %CALL WBSJA(u, 0.d0, nmax, 3, ju)
                                        for ii = 0:nmax0
                                            js(ii+1) = besselj(ii,arg2);
                                        end
                                        %z(0) = (acyl/(li*li - ls*ls))*(li*js(0)*ju(1) - ls*ju(0)*js(1));  % (a10 of [7])
                                        z(1) = (acyl/(li*li - ls*ls))*(li*js(1)*ju(2) - ls*ju(1)*js(2));  % (a10 of [7])
                                        %z(1) = (acyl/(li*li - ls*ls))*(li*js(1)*ju(2) - ls*ju(1)*js(2));  % (a11 of [7])
                                        z(2) = (acyl/(li*li - ls*ls))*(li*js(2)*ju(3) - ls*ju(2)*js(3));  % (a11 of [7])
                                        zmin = z(2);
                                        if (ls == 0.0)
                                            z(3) = z(1) - 2*ju(2)/li*acyl*0.5;
                                        else
                                            z(3) = z(1) - 2*ju(2)*js(2)/(ls*li);  % (a13 of [7])
                                        end
                                        for n = 4:nmax0 + 2
                                            if (ls == 0.0)
                                                z(n) = z(n - 2);
                                            else
                                                z(n) = z(n - 2) - 2*(n-1-1)*ju(n - 1)*js(n - 1)/(ls*li);  % (a13 of [7])
                                            end
                                        end
                                        a(1) = (k/(2.0*li))*(zmin - z(2));    %  (38 of [7])
                                        b(1) = (k/(2.0*li))*(zmin + z(2));
                                        % Computations in 42 of [7]
                                        ff = kk*lcyl2*mi*(ec - 1.0);    % to be used in (42) of [7]
                                        ff1vv = ev(1)*(b(1)*cthetsl*( - cthetil) - sthetsl*z(1));
                                        ff1vh = 0.d0;
                                        ff1hv = 0.d0;
                                        ff1hh = b(1)*eta*hh(1);
                                        %for n = 1:nmax0+1
                                        for n = 2:nmax0+1
                                            a(n) = (k/(2*li))*(z(n-1) - z(n+1));      % (38 of [7])
                                            b(n) = (k/(2*li))*(z(n-1) + z(n+1));      % (38 of [7])
                                            svv = 2.0*((ev(n)*(-cthetil)*b(n) - 1i*eta*hv(n)*a(n))*cthetsl - sthetsl*ev(n)*z(n));
                                            svh = ((eh(n)*( - cthetil)*b(n) - 1i*eta*hh(n)*a(n))*cthetsl - sthetsl*eh(n)*z(n));
                                            shv = (eta*hv(n)*b(n) + 1i*ev(n)*( - cthetil)*a(n));
                                            shh = 2.0*(eta*hh(n)*b(n) + 1i*eh(n)*( - cthetil)*a(n));
                                            ff1vv = ff1vv + svv*cos((n-1)*(phisl - phil));
                                            ff1vh = ff1vh + svh*sin((n-1)*(phisl - phil));
                                            ff1hv = ff1hv + shv*sin((n-1)*(phisl - phil));
                                            ff1hh = ff1hh + shh*cos((n-1)*(phisl - phil));
                                        end
                                        % (42) of [7]
                                        ff1vv = ff1vv*ff;
                                        ff1vh = 2.0*1i*ff1vh*ff;
                                        ff1hv = 2.0*1i*ff1hv*ff;
                                        ff1hh = ff1hh*ff;
                                        % From local frame to absolute reference frame
                                        if (beta==0.0 && gamma==0.0)
                                           fvv = ff1vv;
                                           fvh = ff1vh;
                                           fhv = ff1hv;
                                           fhh = ff1hh;
                                        else
                                           fv1 =  - ff1vv*tvs + ff1hv*ths;       % (49 of [7])
                                           fv2 =  - ff1vh*tvs + ff1hh*ths;
                                           fh1 =  - ff1vv*ths - ff1hv*tvs;
                                           fh2 =  - ff1vh*ths - ff1hh*tvs;
                                           dd = tt*tts;
                                           fvv = ( - fv1*tvi + fv2*thi)/dd;      % (48 of [7])
                                           fvh = ( - fv2*tvi - fv1*thi)/dd;
                                           fhv = ( - fh1*tvi + fh2*thi)/dd;
                                           fhh = ( - fh1*thi - fh2*tvi)/dd;
                                        end
                                        % Computations of scattering functions and averaging over Eulerian angles and angular intervals of theta and theta_s
                                        fvvfvv = fvv*fvv';
                                        fvhfvh = fvh*fvh';
                                        fhvfhv = fhv*fhv';
                                        fhhfhh = fhh*fhh';
                                        xff(k3,imax,1,1) = xff(k3,imax,1,1) + ww*cdirez*fvvfvv/nabg;
                                        xff(k3,imax,1,2) = xff(k3,imax,1,2) + ww*cdirez*fvhfvh/nabg;
                                        xff(k3,imax,2,1) = xff(k3,imax,2,1) + ww*cdirez*fhvfhv/nabg;
                                        xff(k3,imax,2,2) = xff(k3,imax,2,2) + ww*cdirez*fhhfhh/nabg;
                                        %Extinction cross sections (51 of [7])
                                        if (icode == 2 && i == j && k2 == k1 && k3 == 1) % forward direction
                                            fv1fs = real(ff1vv)+1i*abs(imag(ff1vv));
                                            fh1fs = real(ff1hh)+1i*abs(imag(ff1hh));
                                            if (beta == 0.0 && gamma == 0.0)
                                                fvfs = fv1fs;
                                                fhfs = fh1fs;
                                            else
                                                fvfs = (fv1fs*tvs*tvi + fh1fs*ths*thi)/dd;          % (49 of [7])
                                                fhfs = (fv1fs*ths*thi + fh1fs*tvs*tvi)/dd;          % (49 of [7])
                                            end
                                            % ke values RETURNed to the main program are divided by absctheta
                                            % as requested in formula (8) of [1]
                                            cke = om(k1)/(2*absctheta*nabg);
                                            ke(j,1) = ke(j,1) + ww*cke*abs(4.d0*pi/k*imag(fvfs));
                                            ke(j,2) = ke(j,2) + ww*cke*abs(4.d0*pi/k*imag(fhfs));
                                        end
                                        if (iflag == 0)
                                            rap1 = (xff(k3,imax,1,1) - xff(k3,imax - 1,1,1))/xff(k3,imax,1,1);
                                            rap2 = (xff(k3,imax,2,2) - xff(k3,imax - 1,2,2))/xff(k3,imax,2,2);
                                            if (abs(rap1) < 1.d-3 && abs(rap2) < 1.d-3) 
                                                iflag = 1;
                                                imax = imax - 1;
                                                nmax = 2^imax;
                                                condition = 0;
                                            else
                                                if (icode == 2 && i == j && k2 == k1 && k3 == 1) 
                                                    ke(j,:) = 0.d0;
                                                end
                                                break
                                            end
                                        end
                                    end %k3Loop
                                    if (iflag == 0)
                                        break;
                                    end
                                end %k2Loop 
                                
                                % de este no estoy super seguro
%                                 if (iflag == 0)
%                                     break;
%                                 end
                                % fin de este no estoy super seguro
                            
                                if (k3==nm1)
                                    condition = 0;
                                    %break
                                    if (k2==ns)
                                        condition = 0;
                                        if (k1==ns) 
                                            condition = 0;
                                            if (igamma==ngamma) 
                                                condition = 0;
                                                if (ibeta==nbeta) 
                                                    condition = 0;
                                                    if (ialpha==nalpha)
                                                        condition = 0;
                                                        if (j==nij)
                                                            condition = 0;
                                                            if (i==nij)
                                                                condition = 0;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                
                            end %while

                        end %END DO k1Loop

                    end %END DO gammaLoop

                end %END DO betaLoop

            end %END DO alphaLoop
                %  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
            %     DO k3 = nm/2 + 2,nm
            %       k3a = 2*nm1 - k3
            %       xff(k3,imax,:,:) = xff(k3a,imax,:,:)
            %     END DO
            %     %
            %     %      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
            %     %
            %     DO istoki = 1,2
            %       DO istoks = 1,2
            %         xx(:) = CMPLX(xff(:,imax,istoks,istoki),0.)
            %         CALL CFSTFT(NF,xx)
            %         % normalization
            %         xx(1) = xx(1) / NM
            %         xx(2:) = xx(2:)*2/NM
            %         s(i,j,icode,1:nm1,istoks,istoki) = real(xx(1:nm1))	   
            %       END DO
            %     END DO

            for k3 = nm/2 + 2:nm
                k3a = 2*nm1 - k3;
                xff(k3,imax,:,:) = xff(k3a,imax,:,:);
            end
            %      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
            for istoki = 1:2
                for istoks = 1:2
                    xx = xff(:,imax,istoks,istoki);
                    %CALL CFSTFT(NF,xx)
                    xx = fft(xx);
                    % normalization
                    xx(1) = xx(1)/nm;
                    %xx(2:) = xx(2:)*2/NM
                    xx(2:end) = xx(2:end)*2/nm;
                    s(i,j,icode,1:nm1,istoks,istoki) = real(xx(1:nm1));	
                end
            end    
        end %DO iLoop
	end % DO icodeLoop
end % DO jLoop

