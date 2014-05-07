%This subroutine contains a subroutine to calculate litter permittivity
%
%SUBROUTINE LITTER(vmoist,lmoist,biomass,etc)
function etc = MOD_SOIL_FUNCTIONS_LITTER(vmoist,lmoist,biomass)

global nn;
global f;
global dsdwf;
global PERMITTIVITY;
global NIJ;

% REAL, INTENT(IN)      :: vmoist,lmoist
% REAL, INTENT(IN)      :: biomass
% COMPLEX, INTENT(OUT)  :: etc

BIOMAX = 201;  % bio = (0 - 2) [g/cm^2] => 20 [Kg/m^2]
NUM_SMC = 201;  % smc_eq = (0 - 0.4)
NUM_K_MAX = 201;       % k = (0 - 200)
%
% The maximum value that the exponential value can assume is limited to 200,
% higher values do not give further improvements to the function trends 
%
DELTA_BIO = .01;
DELTA_K = 1.;
%
% The values of delta k is set to 1, because finer values do not give better
% sensitivity to the function trend
%
DELTA_SMC = .002;
LITTER_SLOPE = 8.5;    % Slope of line which link the litter biomass [g]
                                               % to litter thickness [cm] 
%REAL, DIMENSION(NN)            :: RS,RL
RS = zeros(nn); RL = zeros(nn);

%REAL, DIMENSION(NN,BIOMAX)     :: r0   
r0 = zeros(nn,BIOMAX);
%COMPLEX, DIMENSION(NN,BIOMAX)  :: RT0
RT0 = zeros(nn,BIOMAX);

% Variables used for the interpolation
%REAL, DIMENSION(NN,BIOMAX)     :: fexp
fexp = zeros(nn,BIOMAX);
%REAL, DIMENSION(NUM_K_MAX)     :: k,RMS_K
k = zeros (NUM_K_MAX); RMS_K = zeros (NUM_K_MAX);
%REAL, DIMENSION(NN)            :: a,b
a = zeros (nn); b = zeros(nn);
%REAL, DIMENSION(BIOMAX)        :: deltax
deltax = zeros(BIOMAX);
%INTEGER, DIMENSION(NN,1)       :: minimum
minimum = zeros(nn,1);
%INTEGER                        :: h

% Variables used for the equivalent permittivity
%COMPLEX, DIMENSION(NUM_SMC)    :: ec
ec = zeros(NUM_SMC);
%REAL, DIMENSION(NUM_SMC)       :: RMS_SMC
RMS_SMC = zeros(NUM_SMC);
%REAL, DIMENSION(NN)            :: r, Req
r = zeros(nn); Req = zeros(nn);

% Generic variables
% COMPLEX :: elc,k1z,k2z,R1,R2,esp,egc,egceq
% REAL    :: elr,elj,beta21,app,dsdwl,Bio_tot,deltet,theta,smceq,Vv
% REAL    :: deltah,lambda,sin2,alfa2,ka,L2,costheta2,beta2,omega,epsilon0,k02,k0z,mu0
% INTEGER :: j,ja,istoki,ibio,i,ismceq

deltet = pi/nn;
lambda = 30/f;	
dsdwl = dsdwf;

% Soil permittivity
%CALL DOBSOIL(vmoist,egc)
egc = MOD_PERMITTIVITY_DOBSOIL(vmoist);

% Vegetation permittivity 
%CALL VEGEDIEL(dsdwl,lmoist,elr,elj,PERMITTIVITY)
[elr,elj] = MOD_PERMITTIVITY_VEGEDIEL(dsdwl,lmoist,PERMITTIVITY);
elc = (elr -1i*elj);

% Litter permittivity
Vv = (1 - lmoist*(1 - dsdwl))/ (dsdwl * LITTER_SLOPE);

elc = (1 + Vv*(sqrt(elc) - 1))^2;

% Estimation of litter reflectivity (coherent approach)
for ibio = 1:BIOMAX %Biomassa_tot:  
      Bio_tot = (ibio - 1)*DELTA_BIO;    % [g/cm^2]  
      deltah = Bio_tot * LITTER_SLOPE;   % [cm]

      for istoki = 1:2  
            for j = 1:NIJ
                  ja = NIJ*(istoki - 1) + j;
                  theta = (j - .5)*deltet; 
                  sin2 = (sin(theta))^2;
                  costheta2 = sqrt(1 - (1/elc)*sin2);
                  alfa2 = (2*pi/lambda)*abs(imag(sqrt(elc)));   
                  ka = 2*alfa2;                                
                  L2 = exp((ka*deltah)/costheta2);              
                  beta2 = (2*pi/lambda)*real(sqrt(elc));	 
                  beta21 = beta2/costheta2;
                  app = 2*beta21*deltah;
                  esp = (0. -1i*app);
                  omega = 2*pi*f*1.e9;
                  epsilon0 = 8.854e-12;
                  mu0 = 4*pi*1.e-7;
                  k02 = omega^2*epsilon0*mu0;
                  k0z = sqrt(k02 - k02*sin2);
                  k1z = sqrt(k02*elc - k02*sin2);
                  k2z = sqrt(k02*egc - k02*sin2);
                  %
                  % R1 air - litter reflectivity
                  % R2 litter - soil reflectivity
                  %
                  if (istoki == 1)
                        RS(ja) = (abs((egc*cos(theta) - sqrt(egc - sin2))/(egc*cos(theta) + sqrt(egc - sin2))))^2;
                        RL(ja) = (abs((elc*cos(theta) - sqrt(elc - sin2))/(elc*cos(theta) + sqrt(elc - sin2))))^2;
                        R1 = (elc*k0z - k1z)/(elc*k0z + k1z);
                        R2 = (egc*k1z - elc*k2z)/(egc*k1z + elc*k2z);
                  else
                        RS(ja) = (abs((cos(theta) - sqrt(egc - sin2))/(cos(theta) + sqrt(egc - sin2))))^2;
                        RL(ja) = (abs((cos(theta) - sqrt(elc - sin2))/(cos(theta) + sqrt(elc - sin2))))^2;
                        R1 = (k0z - k1z)/(k0z + k1z);
                        R2 = (k1z - k2z)/(k1z + k2z);
                  end
                  RT0(ja,ibio) = (R1 + (R2/L2)*exp(esp))/(1 + ((R1*R2*exp(esp))/L2));
                  r0(ja,ibio) = RT0(ja,ibio)*conj(RT0(ja,ibio));
            end
      end
end
%
% An exponential function has been used to reproduct the coherent reflectivity function
% f(x) = A e^(-kx) + B
% f(0) = A + B   is the reflectivity coefficient of soil
% f(inf) = B     is the reflectivity coefficient of litter
%
% Interpolation of coherent reflectivity function
for istoki = 1:2   
    for j = 1:NIJ
        ja = NIJ*(istoki - 1) + j;
        b(ja) = RL(ja);             
        a(ja) = r0(ja,1) - b(ja);
        h = 1;
        for i = 1:NUM_K_MAX
              k(h) = (i - 1)*DELTA_K;
              for ibio = 1:BIOMAX
                    Bio_tot = (ibio - 1)*DELTA_BIO;                         
                    fexp(ja,ibio) = a(ja)*exp(-k(h)*Bio_tot) + b(ja); 
                    deltax(ibio) = (r0(ja,ibio) - fexp(ja,ibio))^2;
              end
              %RMS_K(h) = sqrt(SUM(deltax(:))/BIOMAX);
              RMS_K(h) = sqrt(sum(deltax(:))/BIOMAX);
              h = h + 1;
        end
        %minimum(ja,:) = MINLOC(RMS_K(:))
        % minimum(ja,:) = min(RMS_K(:)); %version posta
        minimum(ja) = min(RMS_K(:))+1;
    end
end
%
% Estimation of equivalent permittivity
i = 1;
for ismceq = 1:NUM_SMC
      smceq = 0.002 + (ismceq - 1)*DELTA_SMC;
      %CALL DOBSOIL(smceq,egceq)
      egceq = MOD_PERMITTIVITY_DOBSOIL(smceq);      
      ec(i) = egceq;
      for istoki = 1:2   
            for j = 1:NIJ 
                  ja = NIJ*(istoki - 1) + j;
                  theta = (j - .5)*deltet;    
                  sin2 = sin(theta)*sin(theta);
                  % Reflectivity value for the input biomass value
                  r(ja) = a(ja)*exp(-k(minimum(ja,1))*biomass) + b(ja);
                  % Reflectivity of bare soil with equivalent SMC
                  if (istoki == 1) 
                        Req(ja) = ((abs((ec(i)*cos(theta) - sqrt(ec(i) - sin2))/(ec(i)*cos(theta) + sqrt(ec(i) - sin2))))^2); 
                  else
                        Req(ja) = ((abs((cos(theta) - sqrt(ec(i) - sin2))/(cos(theta) + sqrt(ec(i) - sin2))))^2); 
                  end
            end
      end   
      RMS_SMC(i) = sqrt(sum((Req(:) - r(:)).^2)/nn);
      i = i + 1;
end
%minimum(1,:) = MINLOC(RMS_SMC(:))
%%minimum(1,:) = min(RMS_SMC(:)); %original
minimum(1,:) = min(RMS_SMC(:))+1;
etc = ec(minimum(1,1));
