%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUBROUTINE PHYSICAL_OPTICS
% CALLED BY: MAIN
% Computes the scatter matrix and the extinction vector using the
% Physical - Optics approximation for an ensemble of circular discs.
%
%SUBROUTINE PHYSICAL_OPTICS(a,b,l,amoi,dsdw,s,ke)

function [s,ke] = MOD_VEGERARION_FUNCTIONS_PHYSICAL_OPTICS(a,b,l,amoi,dsdw)

global nm;
global nn;
global nij;
global nm1;
global f;
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

% Dummy variables declaration
%REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
s = zeros(nij,nij,2,nm1,2,2); %warning: no estoy seguro de nm1!
%REAL, DIMENSION(:,:), INTENT(OUT)			    :: ke
ke = zeros(nn/2,2);
%REAL, INTENT(IN)							            :: a,b,l,amoi,dsdw
% Local declaration variables
%COMPLEX, DIMENSION(NM)	:: xx
%xx = zeros(nm);
%COMPLEX, DIMENSION(3,3)	:: v2p,v2m,h2,h2ehp,v2evp,h2ehm,v2evm,alev,epiu,emen,epsin,emsin
v2p = zeros(3,3); v2m = zeros(3,3); h2 = zeros(3,3); h2ehp= zeros(3,3); v2evp= zeros(3,3); h2ehm= zeros(3,3); 
v2evm = zeros(3,3); alev= zeros(3,3); epiu= zeros(3,3); emen= zeros(3,3); epsin= zeros(3,3); emsin= zeros(3,3);
%COMPLEX, DIMENSION(3,1) :: uvepsp,uvepsm,fv,fh
uvepsp= zeros(3,1); uvepsm = zeros(3,1); fv = zeros(3,1); fh= zeros(3,1); 
%COMPLEX, DIMENSION(2,2) :: 
eq= zeros(2,2);
%COMPLEX, DIMENSION(2)   :: r,t,sig,sinc,tet
r = zeros(2,1); t = zeros(2,1);sig = zeros(2,1);sinc = zeros(2,1);tet = zeros(2,1);
%COMPLEX					        :: ubetz,betmen,betpiu,cscatc,alfexp,psi,denom,ecom,radeps,fvv,fhv,fvh,fhh
%REAL, DIMENSION(NM,2,2) :: xff
xff = zeros(nm,2,2);
%REAL, DIMENSION(3,1)    :: h,v,av,ah
h = zeros(3,1); v = zeros(3,1); av = zeros(3,1); ah = zeros(3,1);
%REAL, DIMENSION(1,3)	  :: ht,vt
ht = zeros(1,3); vt = zeros(1,3);
%REAL, DIMENSION(3)		  :: vs,hs,q,o,ri,,hvn,rx,ry
vs = zeros(3,1);hs = zeros(3,1); q = zeros(3,1); o = zeros(3,1); ri = zeros(3,1); rn= zeros(3,1); hvn = zeros(3,1); rx = zeros(3,1); ry= zeros(3,1);
%REAL, DIMENSION(2)		  :: icos,jcos,q2
q2 = zeros (2,1);

%REAL					:: delthet,delphis,er,ej,cscatr,cscatj,BESJ1
%REAL					:: lambda,k,kk,j1qt,gamma,beta,alpha,theta,thetas,thetaj,thetais,stheta,ctheta,absctheta,sthetas,abscthetas
%REAL					:: cthetas,o0,cdirez,phis,sphis,cphis,salpha,calpha,sbeta,cbeta,sgamma,cgamma,cthetau,tvi,thi,tt,sthetau
%REAL					:: on,qt,aqt,vscat,cke,fvvfvv,fvhfvh,fhvfhv,fhhfhh,alpha1,beta1,gamma1,dalpha,dbeta,dgamma

%INTEGER				:: ialpha,ibeta,igamma,nalpha,nbeta,ngamma,i,j,k1,k2,k3,icode,istoki,istoks,kj,nabg,k3a

%DATA jcos/ -1, -1/
jcos = [-1, -1];
%DATA icos/1, -1/
icos = [1, -1];

if(a > l) 
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

delthet = pi/(2*nij);     % interval amplitude in theta and theta_s
delphis = 2.*pi/nm;       % interval amplitude in phi_s - phi
%
lambda = 30./f;           % wavelength, in cm
k = 2.*pi/lambda;         % wavenumber, in cm^ - 1
kk = k*k;
% 
%CALL VEGEDIEL(dsdw,amoi,er,ej,PERMITTIVITY) 
[er, ej] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,PERMITTIVITY);
                               % permittivity computation
%ecom = CMPLX(er,ej)
ecom = er + 1i*ej;
%
cscatr = kk*a*b*l*(er - 1);     % (22) of [5]
cscatj = kk*a*b*l*ej;
%cscatc = CMPLX(cscatr,cscatj)
cscatc = cscatr + 1i*cscatj;

%
sig(1) = 1./ecom;               % (20c) of [5]
sig(2) = 1;
radeps = sqrt(ecom);

for j = 1:nij %jLoop: DO           % j = incidence angular interval index
    tic;
  thetaj = (j - 1)*delthet;     % lower limit of angular interval in theta
  ke(j,:) = 0.;                 % initialization to compute extinction coefficients
  %
  for icode = 1:2 %icodeLoop: DO
    % icode = 1  : upper half - space scattering
    % icode = 2  : lower half - space scattering
    %

    for i = 1:nij %iLoop: DO           % i = scattering angular interval index
      thetais = (i - 1)*delthet;  % lower limit of angular interval in theta_s
      xff(1:nm1,:,:) = 0.;

      for k1 = 1:2 %k1Loop: DO         % integration within jth theta interval using Gauss - Legendre technique
        theta = thetaj + .5*delthet*(1. + chi(k1));
        stheta = sin(theta);
        absctheta = cos(theta);
        ctheta = jcos(icode)*absctheta;
        % Incidence direction vector
        ri(1) = stheta;
        ri(2) = 0.;
        ri(3) = ctheta;
        % Incidence vertical polarization vector
        v(1,1)  = ctheta;
        v(2,1) = 0.0;
        v(3,1) = -stheta;
        % Incidence horizontal polarization vector
        h(1,1) = 0.0;
        h(2,1)  = 1.0;
        h(3,1)  = 0.0;

        for k2 = 1:2 %k2Loop:      % integration within ith theta_s interval using Gauss - Legendre technique
          thetas = thetais + .5*delthet*(1. + chi(k2));
          sthetas = sin(thetas);
          abscthetas = cos(thetas);
          cthetas = icos(icode)*abscthetas;
          o0 = om(k1)*om(k2)/4.;  % to be used in the double integration in theta and theta_s
          % The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor, 
          % to be used in formulas (3) and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0;

          % Computations for 0 < phi_s - phi < 180
          for k3 = 1:nm1 %k3Loop: DO 
            phis = (k3 - 1)*delphis;    % phi_s - phi
            sphis = sin(phis);
            cphis = cos(phis);
            %
            % Scattering direction vector
            o(1) = sthetas*cphis;
            o(2) = sthetas*sphis;
            o(3) = cthetas;
            % Scattering vertical polarization vector
            vs(1) = cthetas*cphis;
            vs(2) = cthetas*sphis;
            vs(3) =  -sthetas;
            % Scattering horizontal polarization vector
            hs(1) =  -sphis;
            hs(2) = cphis;
            hs(3) = 0.;

            for ialpha = 1:nalpha %alphaLoop: DO 
              alpha = alpha1+(ialpha - 1)*dalpha;
              salpha = sin(alpha);
              calpha = cos(alpha);

                for ibeta = 1:nbeta %betaLoop: DO 
                  beta = beta1+(ibeta - 1)*dbeta;
                  sbeta = sin(beta);
                  cbeta = cos(beta);

                  for igamma = 1:ngamma %gammaLoop: DO 
                    gamma = gamma1+(igamma - 1)*dgamma;
                    sgamma = sin(gamma);
                    cgamma = cos(gamma);
                    % Unit vector perpendicular to disc
                    rn(1) = calpha*sbeta*cgamma + salpha*sgamma;
                    rn(2) = salpha*sbeta*cgamma - calpha*sgamma;
                    rn(3) = cbeta*cgamma;
                    % Direction of x axis in local reference frame
                    rx(1) = cbeta*calpha;
                    rx(2) = cbeta*salpha;
                    rx(3) = -sbeta;
                    % Direction of y axis in local reference frame
                    ry(1) = calpha*sbeta*sgamma - salpha*cgamma;
                    ry(2) = calpha*cgamma + salpha*sbeta*sgamma;
                    ry(3) = cbeta*sgamma;
                    % Cosine of incidence direction in local reference frame
                    cthetau = ri(1)*rn(1) + ri(3)*rn(3);
                    if(cthetau > 0.) 
                      cthetau = - cthetau;
                      rn(:) =  - rn(:);
                      rx(:) =  - rx(:);
                      ry(:) =  - ry(:);
                    end
                    sthetau = sqrt(1. - cthetau*cthetau);
                    ubetz = sqrt(ecom - sthetau*sthetau); % page 1257 of [5]
                    alfexp = 0. + 1i*k*l/2.;
                    psi = 4*alfexp*ubetz;
                    betmen = cthetau - ubetz;
                    betpiu = cthetau + ubetz;

                    tvi = (salpha*sgamma + calpha*sbeta*cgamma)*ctheta - cbeta*cgamma*stheta;
                    thi = salpha*sbeta*cgamma - calpha*sgamma;    
                    tt = sqrt(tvi*tvi + thi*thi);
                    % Vertical and horizontal polarization vectors in local reference frame
                    if(beta == 0. && gamma == 0.)
                      av = v(:,:);
                      ah = h(:,:);
                    else
                      av(:,:) = -(tvi*v(:,:) + thi*h(:,:))/tt;
                      ah(:,:) = (thi*v(:,:) - tvi*h(:,:))/tt;
                    end
                    ht(:,:) = ah(:,:)';
                    vt(:,:) = av(:,:)';
                    %CALL PRODV(ah,rn,hvn)									        % 18a of [5]
                    hvn = cross(ah,rn);
                    uvepsp(:,1) = (ubetz - cthetau)*hvn(:);
                    uvepsp(:,1) = (av(:,1) + uvepsp(:,1))/radeps;	% 18b of [5]
                    uvepsm(:,1) =  - (ubetz + cthetau)*hvn(:);
                    uvepsm(:,1) = (av(:,1) + uvepsm(:,1))/radeps;	% 18b of [5]

                    for kj = 1:2
                      denom = cthetau - sig(kj)*ubetz;
                      r(kj) = (cthetau + sig(kj)*ubetz)/denom;			% 20a of [5]
                      t(kj) = 2.*sqrt(sig(kj))*cthetau/denom;			% 20b of [5]
                    end
                    for kj = 1:2
                      denom = 1. - r(kj)*r(kj)*exp(psi);
                      eq(kj,1) =  - t(kj)*r(kj)*exp(psi)*exp(alfexp*betmen)/denom;   % 19a of [5]
                      eq(kj,2) = t(kj)*exp(alfexp*betpiu)/denom;					            % 19b of [5]
                    end
                    %h2(:,:) = MATMUL(ah(:,:),ht(:,:))
                    h2 = ah*ht;
                    %v2p(:,:) = MATMUL(uvepsp(:,:),vt(:,:))
                    v2p = uvepsp*vt;
                    %v2m(:,:) = MATMUL(uvepsm(:,:),vt(:,:))
                    v2m = uvepsm*vt;
                    h2ehp(:,:) = h2(:,:)*eq(2,1);
                    v2evp(:,:) = v2p(:,:)*eq(1,1);
                    h2ehm(:,:) = h2(:,:)*eq(2,2);
                    v2evm(:,:) = v2m(:,:)*eq(1,2);
                    epiu(:,:) = h2ehp(:,:) + v2evp(:,:);								    % 23a of [5]
                    emen(:,:) = h2ehm(:,:) + v2evm(:,:);
                    %on = DOT_PRODUCT(o(:),rn(:))
                    on = o'*rn;
                    tet(1) = k*l*(on + ubetz)/2.;						% 23c of [5]
                    tet(2) =  - k*l*(on - ubetz)/2.;
                    %sinc(:) = CSIN(tet(:))/tet(:);
                    sinc = sin(tet)./tet;
                    epsin(:,:) = epiu(:,:)*sinc(2);
                    emsin(:,:) = emen(:,:)*sinc(1);
                    % Scattering amplitude matrix in local reference frame
                    alev(:,:) = epsin(:,:) + emsin(:,:);
                    % Wavenumber difference vector
                    q(:) = o(:) - ri(:);
                    % Transformation to local reference frame
                    %q2(1) = DOT_PRODUCT(q(:),rx(:))
                    q2(1) = q'*rx;
                    q2(2) = q'*ry;
                    qt = k*sqrt(a*a*q2(1)*q2(1) + b*b*q2(2)*q2(2));
                    % Form factor of circular disc ((23d) + (23e) + (23f) of [5])
                    aqt = qt;
                    if(qt == 0.) 
                      vscat = .5/2;
                    else
                      %j1qt = BESJ1(aqt)
                      j1qt = besselj(1,aqt);
                      vscat = j1qt / (2. * aqt);
                    end
                    % Transformation to principal reference frame  ((28)of [4])
                    %fv(:,:) = MATMUL(alev(:,:),v(:,:))
                    fv(:,:) = alev*v;
                    fh(:,:) = alev*h;
                    % Elements of the scattering amplitude matrix
                    fvv = vs(1)*fv(1,1) + vs(2)*fv(2,1) + vs(3)*fv(3,1);
                    fhv = hs(1)*fv(1,1) + hs(2)*fv(2,1) + hs(3)*fv(3,1);
                    fvh = vs(1)*fh(1,1) + vs(2)*fh(2,1) + vs(3)*fh(3,1);
                    fhh = hs(1)*fh(1,1) + hs(2)*fh(2,1) + hs(3)*fh(3,1);
                    fvv = fvv*vscat*cscatc;
                    fhv = fhv*vscat*cscatc;
                    fvh = fvh*vscat*cscatc;
                    fhh = fhh*vscat*cscatc;
                    % Extinction cross - section computation using the forward scattering theorem
                    if(icode == 2 && i == j && k2 == k1 && k3 == 1) 
                      % ke values RETURNed to the main program are divided by absctheta as requested in formula (8) of [1]
                      cke = om(k1)/(2*absctheta*nabg);
                      ke(j,1) = ke(j,1) + cke*abs(4*pi/k*imag(fvv));
                      ke(j,2) = ke(j,2) + cke*abs(4*pi/k*imag(fhh));
                    end
                    %
                    %fvvfvv = fvv*CONJG(fvv)
                    fvvfvv = fvv*fvv';
                    fvhfvh = fvh*fvh';
                    fhvfhv = fhv*fhv';
                    fhhfhh = fhh*fhh';
                    %
                    xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nabg;
                    xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nabg;
                    xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nabg;
                    xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nabg;

                  end % DO gammaLoop
                end % DO betaLoop
            end % DO alphaLoop
          end % DO k3Loop
        end % DO k2Loop
        toc
      end % DO k1Loop
      %
      %  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      %
      for k3 = nm/2 + 2:nm
        k3a = 2*nm1 - k3;
        xff(k3,:,:) = xff(k3a,:,:);
      end
      %
      %      Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
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

    end %DO iLoop

end % DO icodeLoop
end % DO jLoop




% %*********************************************************
% SUBROUTINE PRODV(a,b,c)
% % called by PHYSICAL_OPTICS
% % vector product;
% REAL, INTENT(IN)  :: a(3,1),b(3)
% REAL, INTENT(OUT) :: c(3)
% 
% c(1) = a(2,1)*b(3) - a(3,1)*b(2)
% c(2) = a(3,1)*b(1) - a(1,1)*b(3)
% c(3) = a(1,1)*b(2) - a(2,1)*b(1)
% 
% RETURN 
% END SUBROUTINE PRODV
% 
% END SUBROUTINE PHYSICAL_OPTICS
%