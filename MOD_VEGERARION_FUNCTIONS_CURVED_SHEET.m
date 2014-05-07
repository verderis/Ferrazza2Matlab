%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUBROUTINE CurvedSheet
% CALLED BY: MAIN
% Computes the scatter matrix and the extinction vector using the
% curved sheet approximation for a whole leaf.
%
%SUBROUTINE CURVED_SHEET(ac,bc,lc,amoi,dsdw,s,ke)
function [s,ke] = MOD_VEGERARION_FUNCTIONS_CURVED_SHEET(ac,bc,lc,amoi,dsdw)

% Dummy variables declaration
%REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT)	:: s
s = zeros(nij,nij,2,nm1,2,2); %warning: no estoy seguro de nm1!

%REAL, DIMENSION(:,:),INTENT(OUT)			    :: ke
ke = zeros(nn/2,2);

%REAL,INTENT(IN)														:: ac,bc,lc,amoi,dsdw
% Local variables declaration
%COMPLEX, DIMENSION(3,2)		:: SS,Sl,Sm
SS = zeros(3,2); Sl = zeros(3,2); Sm = zeros(3,2);

%COMPLEX, DIMENSION(NM)		:: xx
%COMPLEX										:: fvv,fvh,fhv,fhh,molt,RR

%REAL, DIMENSION(NM,2,2)		:: xff
xff = zeros(nm,2,2);

%REAL, DIMENSION(3,2)			:: es
es = zeros(3,2);

%REAL		:: delta_theta,delta_phi,thetamin,phi_l,theta,thetas,phis,theta_f,phi_f,stheta
%REAL		:: thetaj,thetais,delthet,delphis,abscthetas,cdirez,ctheta,sthetas,cthetas,sphis,cphis
%REAL		:: alpharad,den_r,den_i,thetamax,o0,cke,fvvfvv,fhhfhh,fvhfvh,fhvfhv,r,k,z,soglia,er,ei,lambda,a,b,l

%INTEGER	:: nphi,dphi,l2,m2_1,m2,ck,test,i,j,k1,k2,k3,istoki,istoks,k3a,icode,i_phi,i_elev,i_azim

%***************************************************************************************
% ANGULAR CURVATURE OF DIELECTRIC SHEET
%***************************************************************************************
%INTEGER, PARAMETER  :: ANGLE_CURVATURE = 90			% Curvature angle of dielectric curved sheet

ANGLE_CURVATURE = 90;

%fvv = (0.,0.)
%fvh = (0.,0.)
%fhv = (0.,0.)
%fhh = (0.,0.)

% conversion from centimeter to meter of sheet dimension
a = ac/100.;                
b = bc/100.;
l = lc/100.;

delthet = pi/(2*nij);      % interval amplitude in theta and theta_s
delphis = 2.*pi/nm;        % interval amplitude in phi_s - phi

z = sqrt(14400.*pi*pi);    % Vacuum impedance  377 ohm	
alpharad = ANGLE_CURVATURE*pi/180.; 
r = b/alpharad;            % sphere radius 
lambda = 0.3/f;		        
k = 2.*pi/lambda;		      

% Integration step
soglia = a/5.;	

% Permittivity constant routine
%CALL VEGEDIEL(dsdw,amoi,er,ei,PERMITTIVITY)
[er, ei] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,PERMITTIVITY);

% Definition of constant
%molt = (0.,1.)*((k*r*r)/(2.*pi));    % simplification in k in from SS to Fpq
molt = 1i*((k*r*r)/(2.*pi));    % simplification in k in from SS to Fpq
den_r = k*l*(er - 1);                % (1)  pag.655
den_i = k*l*ei;						          % (1)  pag.655 
RR = 1i*z/(den_r+1i*den_i); % (1)  pag.655
%
% Determination of the integration point with trapezium method, al least two point must be selected
%
%l2 = CEILING(b/soglia)
l2 = ceil(b/soglia);
%m2_1 = CEILING(a/soglia)
m2_1 = ceil(a/soglia);

delta_theta = alpharad/l2;        % Angular amplitude along the meridian
%thetamin = ASIN(a/(2*pi*r))      % Theta value where the sheet begin
thetamin = asin(a/(2*pi*r));      % Theta value where the sheet begin
thetamax = thetamin + alpharad;   % Finishing theta

dphi = 10;                        % Angular step amplitude for the rotation of the leaf
nphi = 360/dphi;		               % Number	of leaf position in Azimuth

for j = 1:nij %jLoop: DO			         % j = incidence angular interval index

  thetaj = (j - 1)*delthet;		   % lower limit of angular interval in theta   
  ke(j,:) = 0.;				           % initialization to compute extinction coefficients

  for icode = 1:2
    % icode = 1  : upper half - space scattering
    % icode = 2  : lower half - space scattering

    for i = 1:nij %iLoop: DO           % i = scattering angular interval index
      thetais = (i - 1)*delthet;  % lower limit of angular interval in theta_s
      % Initialization to compute scatter matrix
      xff(1:nm1,:,:) = 0.;	   

      for k1 = 1:2 %k1Loop: DO 	       % integration within jth theta interval using Gauss - Legendre technique
        theta = thetaj + .5*delthet*(1. + chi(k1));
        stheta = sin(theta);
        ctheta = cos(theta);
        for k2 = 1:2 %k2Loop: DO 	    % integration within ith theta_s interval using Gauss - Legendre technique

          thetas = thetais + .5*delthet*(1. + chi(k2));
          if (icode==2) 
              thetas = pi - thetas; % changing of reference systems [1] e [4]
          end
          sthetas = sin(thetas);
          cthetas = cos(thetas);
          abscthetas = abs(cthetas);
          o0 = om(k1)*om(k2)/4.; % to be used in the double integration in theta and theta_s
          %
          % The "s" functions RETURNed to the main program include the (delthet*stheta/abscthetas) factor, 
          % to be used in formulas (3) and (4) of [1].
          cdirez = (delthet*stheta/abscthetas)*o0;
          %
          % Computations for 0 < phi_s - phi < 180
          for k3 = 1:nm1 %k3Loop: DO 

            phis = (k3 - 1)*delphis;     % phi_s - phi	
            phis = phis + pi;			      % changing of reference systems  [1] e [4]
            sphis = sin(phis);
            cphis = cos(phis);
            %
            % Initialization of scattering polarization vector
            %
            es(1,1) = -cthetas*cphis;
            es(2,1) = -cthetas*sphis;
            es(3,1) = sthetas;
            es(1,2) = -sphis;
            es(2,2) = cphis;
            es(3,2) = 0.;
            %	
            % I suppose 36 position in azimuth plane
					  %
            for i_phi = 1:nphi %iphiLoop: DO 

              phi_l = (i_phi - 1)*dphi*pi/180.;	    % Centre azimuth of the leaf
              theta_f = thetamin - delta_theta;      % Minimum value in elevation

              % S Initialization  
              %SS(:,:) = (0.,0.)
              SS(:,:) = 0;
              %
              % Integration respect to theta
              %        
              for i_elev = 0:l2 %elevation: DO 

                theta_f = theta_f + delta_theta;    

                m2 = ceil(sin(thetamax)/sin(theta_f));
                if(m2 < m2_1) 
                    m2 = m2_1;
                end

                % Determination of angular amplitude of title along the parallel
                delta_phi = a/(m2*r*sin(theta_f));	
                phi_f = phi_l - a/(2.*r*sin(theta_f)) - delta_phi; 

                % Sl (Azimuthal integral) and Sm (point value) initialization
                Sl(:,:) = 0;
                Sm(:,:) = 0;	

                for i_azim = 0:m2 %azimuth: DO    

                  % Integration respect to phi
                  phi_f = phi_f + delta_phi;
                  %
                  % The following procedure has the task to test if the part of considered sheet is in shadow
									% or it is illuminated by the field
                  %
                  %	Test = 1  This is a singularity of the function (Q=0 at denominator), I take the contribution of the next
									%						tile as valid contribution
									%
                  %	Test = 2  The tile is not considered, it is in the shadow part
                  %
                  CALL VERIFICA(theta,theta_f,phi_f,test)

                  if(test == 1) 
                    Sl(:,:) = Sl(:,:) + ck*Sm(:,:);
                    CYCLE
                  elseif(test == 2) 
                    CYCLE
                  END if
                  %   Estimation of S matrix
                  %	
                  CALL TRAPEZIO(k,r,RR,z,theta_f,phi_f,stheta,ctheta,sthetas,cthetas,sphis,cphis,Sm)
                  %	
                  ck = 2;
                  if((i_azim == 0.)||(i_azim == m2))  
                      ck = 1;
                  end
                  Sl(:,:) = Sl(:,:) + ck*Sm(:,:);
                  end % azimuth  

                % Here is concluded the integration for row, or with fixed theta and phi variable
                Sl = .5*delta_phi*Sl;
                % Integration in theta
                ck = 2;
                if((i_elev==0)||(i_elev==l2))  
                    ck = 1;
                end

                SS(:,:) = SS(:,:) + ck*Sl(:,:);
                end % DO elevation

              % Constant moltiplication factor (18) di [1]
              SS(:,:) = .5*delta_theta*SS(:,:);
              % The amplitude scattering function are in [cm]

              fvv = 100.*molt*(SS(1,1)*es(1,1) + SS(2,1)*es(2,1) + SS(3,1)*es(3,1));
              fhv = 100.*molt*(SS(1,1)*es(1,2) + SS(2,1)*es(2,2) + SS(3,1)*es(3,2));
              fvh = 100.*molt*(SS(1,2)*es(1,1) + SS(2,2)*es(2,1) + SS(3,2)*es(3,1));
              fhh = 100.*molt*(SS(1,2)*es(1,2) + SS(2,2)*es(2,2) + SS(3,2)*es(3,2));

              % Extinction cross - section computation using the forward scattering theorem
              if(icode==2 && i==j && k2==k1 && k3==1) 
                % ke values returned to the main program are divided by ctheta as requested in formula (8) of [1]
                cke = om(k1)/(2*ctheta*nphi);
                ke(j,1) = ke(j,1) + cke*abs(4*pi/k*AIMAG(fvv))*100;   % k in cm^ - 1
                ke(j,2) = ke(j,2) + cke*abs(4*pi/k*AIMAG(fhh))*100;   % k in cm^ - 1
              end
              %
              fvvfvv = fvv*fvv';
              fvhfvh = fvh*fvh';
              fhvfhv = fhv*fhv';
              fhhfhh = fhh*fhh';
              %
              xff(k3,1,1) = xff(k3,1,1) + cdirez*fvvfvv/nphi;
              xff(k3,1,2) = xff(k3,1,2) + cdirez*fvhfvh/nphi;
              xff(k3,2,1) = xff(k3,2,1) + cdirez*fhvfhv/nphi;
              xff(k3,2,2) = xff(k3,2,2) + cdirez*fhhfhh/nphi;

              end %END DO iphiLoop

            end %END DO k3Loop

          end %END DO k2Loop

        end %END DO  k1Loop
      %
      %  180 < phi_s - phi < 360 values (azimuthal symmetry is assumed)
      %
      for k3 = nm/2 + 2:nm
        k3a = 2*nm1 - k3;
        xff(k3,:,:) = xff(k3a,:,:);
      end   % k3
      %
      % Fourier transform and "s" matrix computation ("s" values are in radians^ - 1)
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

      end %END DO iLoop

    end %END DO icodLoop

  end % END DO jLoop

end %END SUBROUTINE CURVED_SHEET