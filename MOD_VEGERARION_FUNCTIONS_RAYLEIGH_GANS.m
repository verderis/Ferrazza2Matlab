%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUBROUTINE RAYLEIGH_GANS
% CALLED BY: MAIN
% Computes the scatter matrix and the extinction vector uSINg the
% Rayleigh - Gans approximation for:
% 1) An ensemble of circular discs (a > l);
% 2) An ensemble of needles (l > a).
%

function [s,ke] = MOD_VEGERARION_FUNCTIONS_RAYLEIGH_GANS(a,b,l,amoi,dsdw)

global nm;
global nn;
global nij;
global nm1;
global f;
global ns;
global PERMITTIVITY;
global chi;
global om;

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

%SUBROUTINE RAYLEIGH_GANS(a,b,l,amoi,dsdw,s,ke)
%Dummy varaiables declaration
%REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
%s (i,j,icode,1:NM1,istoks,istoki)
s = zeros(nij,nij,2,nm1,2,2); %warning: no estoy seguro de nm1!
%REAL, DIMENSION(:,:), INTENT(OUT)			    :: ke
ke = zeros(nn/2,2);
% REAL, INTENT(IN)							            :: a,b,l,amoi,dsdw
% Local variabled declaration
%COMPLEX, DIMENSION(3,3)		:: ua,u2,uaa,u1,ucom
ua = zeros(3,3); u2 = zeros(3,3); uaa = zeros(3,3); u1 = zeros(3,3); ucom = zeros(3,3);
%COMPLEX, DIMENSION(3,1)		:: cv,ch,av,ah
cv = zeros(3,1); ch = zeros(3,1); av = zeros(3,1); ah = zeros(3,1);
%COMPLEX, DIMENSION(NM)		:: xx
xx = zeros(nm);

%COMPLEX						        :: cscatc,ecom,fvv,fvh,fhv,fhh

%REAL, DIMENSION(NM,2,2)		:: xff
xff = zeros(nm,2,2);
%REAL, DIMENSION(NN/2,2)		:: ks,ka
ks = zeros(nn/2,2);ka = zeros(nn/2,2);
%REAL, DIMENSION(3,3)		  :: u
u = zeros(3,3);
%REAL, DIMENSION(3,1)		  :: v,q,q2,h
v = zeros(3,1); q = zeros(3,1); q2 = zeros(3,1); h = zeros(3,1);
%REAL, DIMENSION(2,2)		  :: xdirez
xdirez = zeros(2,2);
%REAL, DIMENSION(3)			  :: hs,vs
hs = zeros(3); vs = zeros(3);
%REAL, DIMENSION(2)			  :: jcos,icos,kstoro
%jcos = zeros(2); icos  = zeros(2); 
kstoro = zeros(2,1);

%REAL						:: lambda,k,kk,j1aqt,lqt,kav,kah,kavs,kahs,kdirez1,kdirez2,cscatr,cscatj,ccabs,thetaj,thetais,theta
%REAL						:: stheta,ctheta,absctheta,thetas,sthetas,cthetas,abscthetas,o0,cdirez,phis,sphis,cphis
%REAL						:: alpha1,beta1,gamma1,dalpha,dbeta,dgamma,delthet,delphis,c,er,ej,alpha,salpha,calpha,beta,sbeta,cbeta
%REAL						:: gamma,sgamma,cgamma,qt,aqt,vscat,fvvfvv,fhhfhh,fhvfhv,fvhfvh,rl,rm,rn,BESJ1

%INTEGER					:: nalpha,nbeta,ngamma,nabg,i,j,icode,k1,k2,k3,k3a,ialpha,ibeta,igamma,istoki,istoks

%DATA ua/1.,0.,0.,0.,1.,0.,0.,0.,1./
ua = [1.,0.,0.;0.,1.,0.;0.,0.,1.];
%DATA jCOS/ -1, -1/
jCOS = [-1, -1];
%DATA iCOS/1, -1/
iCOS = [1, -1];

if (a > l)
  % Discs selection 
  alpha1 = ALPHA1DIS;
  nalpha = NALPHADIS;
  dalpha = DALPHADIS;
  beta1 = BETA1DIS;
  nbeta = NBETADIS;
  dbeta = DBETADIS;
  gamma1 = GAMMA1DIS;
  ngamma = NGAMMADIS;
  dgamma = DGAMMADIS;
else
  % Petioles or ears selection
  alpha1 = ALPHA1PET;
  nalpha = NALPHAPET;
  dalpha = DALPHAPET;
  beta1 = BETA1PET;
  nbeta = NBETAPET;
  dbeta = DBETAPET;
  gamma1 = GAMMA1PET;
  ngamma = NGAMMAPET;
  dgamma = DGAMMAPET;
end

nabg = nalpha*nbeta*ngamma;
delthet = pi/(2*nij);	% interval amplitude in theta and theta_s
delphis = 2.*pi/nm;		% interval amplitude in phi_s - phi

lambda = 30./f;				% wavelength, in cm
k = 2.*pi/lambda;			% wavenumber, in cm^ - 1
kk = k*k;
c = 3.*l/4.;					  % ellipsoid small semi - axis (discs), in cm ellipsoid large semi - axis (needles), in cm

%CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY)  % Permittivity computation
[er, ej] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,PERMITTIVITY);
ecom = er + 1i*ej;

%CALL OSBORN(a,b,c,rl,rm,rn)         % Demagnetization factors 2.23 of [3]
[rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_OSBORN(a,b,c);
ua(1,1) = 1./(1. + rl*(ecom - 1));   % Diagonal elements in (24) of [4] and (2) of [6]
ua(2,2) = 1./(1. + rm*(ecom - 1));      
ua(3,3) = 1./(1. + rn*(ecom - 1));
%
cscatr = kk*a*b*l*(er - 1);          % factors in (28) of [4], (5) of [6]
cscatj = kk*a*b*l*ej;
cscatc = cscatr + 1i*cscatj;
ccabs = k*pi*a*b*l*ej;			          % fattore coeff. assorbimento Ref.1 appendice 

for j = 1:nij  %jLoop: DO               % j = incidence angular interval index

  thetaj = (j - 1)*delthet;          % lower limit of angular interval in theta

  ks(j,:) = 0. ; ka(j,:) = 0.;

  for icode = 1:2 %icodeLoop: DO
    % icode = 1  : upper half - space scattering
    % icode = 2  : lower half - space scattering

    for i = 1:nij %iLoop: DO              % i = scattering angular interval index
      
      thetais = (i - 1)*delthet;     % lower limit of angular interval in theta_s
      %
      kstoro(:) = 0.;
      xff(:,:,:) = 0.;

      for k1 = 1:ns %k1Loop: DO               % integration within jth theta interval using Gauss - Legendre technique

        theta = thetaj + .5*delthet*(1. + chi(k1));
        stheta = sin(theta);
        absctheta = cos(theta);
        ctheta = jCOS(icode)*absctheta;

        for k2 = 1:ns %k2Loop: DO 		        % integration within ith theta_s interval using Gauss - Legendre technique

          thetas = thetais + .5*delthet*(1. + chi(k2));
          sthetas = sin(thetas);
          abscthetas = cos(thetas);
          cthetas = iCOS(icode)*abscthetas;
          o0 = om(k1)*om(k2)/4.;         % to be used in the double integration in theta and theta_s
          %
          % The "s" functions returned to the main program include the
          % (delthet*stheta/abscthetas) factor, to be used in formulas (3)
          % and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0;

          % Computations for 0 < phi_s - phi < 180
          for k3 = 1:nm1 %k3Loop: DO 

            phis = (k3 - 1)*delphis;   % phi_s - phi
            sphis = sin(phis);
            cphis = cos(phis);
            %
            %  Computations at page 193 of [4], 141 of [6]. Incidence and scattering polarization vectors.
            v(1,1) = ctheta;           % V pol, incidence
            v(2,1) = 0;
            v(3,1) =  - stheta;
            vs(1) = cthetas*cphis;     % V pol, scattering
            vs(2) = cthetas*sphis;
            vs(3) =  - sthetas;
            h(1,1) = 0;                % H pol, incidence
            h(2,1) = 1.;
            h(3,1) = 0;
            hs(1) =  - sphis;          % H pol, scattering
            hs(2) = cphis;
            hs(3) = 0;
            q(1,1) = sthetas*cphis - stheta;       % wavenumber difference vector (page 192 of [4])
            q(2,1) = sthetas*sphis;
            q(3,1) = cthetas - ctheta;
            %
            xdirez(:,:) = 0.;
            if (k2 == 1 && k3 == 1 && icode == 1 && i == 1)
              kavs = 0.;
              kahs = 0.;
            end
            %                  (averaging over scatterer orientation)
            for ialpha = 1:nalpha %alphaLoop: DO 
              alpha = alpha1 + (ialpha - 1)*dalpha;
              salpha = sin(alpha);
              calpha = cos(alpha);

              for ibeta = 1:nbeta %betaLoop: DO 
                beta = beta1 + (ibeta - 1)*dbeta;
                sbeta = sin(beta);
                cbeta = cos(beta);

                for igamma = 1:ngamma %gammaLoop: DO 
                  gamma = gamma1 + (igamma - 1)*dgamma;
                  sgamma = sin(gamma);
                  cgamma = cos(gamma);
                  % Reference system transformation matrix, (22) of [4]
                  u(1,1) = cbeta*calpha;
                  u(1,2) = salpha*cbeta;
                  u(1,3) = sbeta;
                  u(2,1) =  - calpha*sbeta*sgamma - salpha*cgamma;
                  u(2,2) =  - salpha*sbeta*sgamma + calpha*cgamma;
                  u(2,3) = cbeta*sgamma;
                  u(3,1) =  - calpha*sbeta*cgamma + salpha*sgamma;
                  u(3,2) =  - salpha*sbeta*cgamma - calpha*sgamma;
                  u(3,3) = cbeta*cgamma;
                  % (24) of [4]
                  ucom(:,:) = u(:,:);
                  %u1(:,:) = TRANSPOSE(u(:,:))            % u^ - 1
                  u1 = u';            % u^ - 1
                  %u2(:,:) = MATMUL(u1(:,:),ua(:,:))      % u^ - 1[ - ]
                  u2 = u1*ua;      % u^ - 1[ - ]
                  %uaa(:,:) = MATMUL(u2(:,:),ucom(:,:))   % u^ - 1[ - ]u
                  uaa = u2*ucom;   % u^ - 1[ - ]u
                  %q2(:,:) = MATMUL(u(:,:),q(:,:))
                  q2 = u*q;

                  if (a > c) 
                    % Disc CASE (form factor in (28) of [4])
                    qt = k*sqrt((q2(1,1)*q2(1,1)*a*a) + (q2(2,1)*q2(2,1)*b*b));
                    aqt = qt;
                    if (qt == 0.) 
                      %vscat = .5/2.
                      vscat = .5/2;
                    else					  
                      %j1aqt = BESJ1(aqt)
                      j1aqt = besselj(1,aqt);
                      vscat = j1aqt/(2.*qt);
                    end
                    %
                    % Needle CASE (form factor in (5) of [6])
                  else
                      lqt = k*l*q2(3,1)/2.;
                      if (q2(3,1) == 0.)
                        vscat = 1./4.;
                      else
                        vscat = sin(lqt)/(4.*lqt);
                      end
                  end
                  %
                  cv(:,:) = v(:,:);
                  ch(:,:) = h(:,:);
                  % Transformation to principal reference frame ((28) of [4], (5) of [6])
                  av = uaa*cv; 
                  ah = uaa*ch;
                  % Elements of scattering amplitude matrix
                  fvv = cscatc*vscat*(av(1,1)*vs(1) + av(2,1)*vs(2) + av(3,1)*vs(3));
                  fhv = cscatc*vscat*(av(1,1)*hs(1) + av(2,1)*hs(2) + av(3,1)*hs(3));
                  fhh = cscatc*vscat*(ah(1,1)*hs(1) + ah(2,1)*hs(2) + ah(3,1)*hs(3));
                  fvh = cscatc*vscat*(ah(1,1)*vs(1) + ah(2,1)*vs(2) + ah(3,1)*vs(3));

                  if (k2 == 1 && k3 == 1 && icode ==  1 && i == 1)
                    kav = av(1,1)*av(1,1)' + av(2,1)*av(2,1)' + av(3,1)*av(3,1)';
                    kah = ah(1,1)*ah(1,1)' + ah(2,1)*ah(2,1)' + ah(3,1)*ah(3,1)';
                    kavs = kavs + kav;
                    kahs = kahs + kah;
                  end

                  fvvfvv = fvv*fvv';
                  fvhfvh = fvh*fvh';
                  fhvfhv = fhv*fhv';
                  fhhfhh = fhh*fhh';

                  xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nabg;
                  xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nabg;
                  xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nabg;
                  xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nabg;

                  xdirez(1,1) = xdirez(1,1) + fvvfvv;
                  xdirez(1,2) = xdirez(1,2) + fvhfvh;
                  xdirez(2,1) = xdirez(2,1) + fhvfhv;
                  xdirez(2,2) = xdirez(2,2) + fhhfhh;

                end

              end

            end

            kdirez1 = o0*sthetas*(xdirez(1,1) + xdirez(2,1));
            kdirez2 = o0*sthetas*(xdirez(1,2) + xdirez(2,2));
            %     vedi appendice
            if (k3 == 1 || k3 == nm1) 
                kdirez1 = kdirez1/2.;
                kdirez2 = kdirez2/2.;
            end
            kstoro(1) = kstoro(1) + kdirez1/(nabg*absctheta);
            kstoro(2) = kstoro(2) + kdirez2/(nabg*absctheta);

          end

        end

        if (icode == 1 && i == 1) 
          ka(j,1) = ka(j,1) + kavs*ccabs*om(k1)/(2*nabg*absctheta);   
          ka(j,2) = ka(j,2) + kahs*ccabs*om(k1)/(2*nabg*absctheta);
        end

      end
      kstoro = kstoro*2.;
      ks(j,:) = ks(j,:) + kstoro'*delthet*delphis;  % somma contributi dei tori di scattering 
      %
      %  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      for k3 = nm/2 + 2:nm
        k3a = 2*nm1 - k3;
        xff(k3,:,:) = xff(k3a,:,:);
      end % k3
      %
      %      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
      %	      
      
      for istoki = 1:2
        for istoks = 1:2
          xx = xff(:,istoks,istoki);
          %CALL CFSTFT(NF,xx)
          xx = fft(xx);
          % normalization
          xx(1) = xx(1)/nm;
          %xx(2:) = xx(2:)*2/NM
          xx(2:end) = xx(2:end)*2/nm;
          s(i,j,icode,1:nm1,istoks,istoki) = real(xx(1:nm1));	
        end
      end    

    end % iLoop

  end % icodeLoop

end % jLoop

% Extinction cross senction
ke = ks + ka;