% SUBROUTINE SOILDIEL
% CALLED BY: IEM, CORNER
%
% This SUBROUTINE computes soil permittivity uSINg
% the semiemPIrical model described in [12], pp. 2101 - 2104.
% Inputs: volumetric moisture, frequency.
% Outputs: real (egr) and imaginary (egj) parts of soil permittivity.

function [egr,egj] = MOD_PERMITTIVITY_SOILDIEL(vmoi)

% Dummy variables
% REAL, INTENT(IN)  :: vmoi
% REAL, INTENT(OUT) :: egr,egj

% Local variables
% COMPLEX           :: g2,g3,g4,egc
% REAL			        :: lambda,alpha,beta,ROB,ross,epsss,g1,f0,fhz,g2b

 lambda = 30./F;        % wavelength, in cm
 alpha = 0.65;          % pag. 2103 of [12]
 beta = 1.1;            % pag. 2103 of [12] (tyPIcal value)
 ROB = 1.1;             % tyPIcal bulk density(g/cm^3)
 ross = 2.65;           % soil porosity (g/cm^3)
 epsss = 4.7;           % solid soil permittivity
 g1 = 1 + (ROB/ross)*(epsss^alpha - 1);
 f0 = F/18.64;
 fhz = F*1.E+09;

 if(F <= 4)
       g2 = 4.9 + 74.1/(CMPLX(1.,f0));
       g2b = .107/(2.*PI*fhz*8.854E-12)*(ross - ROB)/(ross*vmoi);
       g2 = (g2 - CMPLX(0.,g2b))^alpha;
 else
        g2 = (4.9 + 74.1/(CMPLX(1.,f0)))^alpha;
 end

 g3 = g1 + (vmoi)^beta*(g2 - 1);
 g4 = (CLOG(g3))/alpha;
 egc = CEXP(g4);
 egr = REAL(egc);
 egj = IMAG(egc);
