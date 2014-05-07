% SUBROUTINE DOBSOIL
% CALLED BY: IEM, CORNER
%
% This SUBROUTINE computes soil permittivity using
% the semiempirical model described in Dobson , 85, en IEEE, vol.23, no.1
% Inputs 
%					TSOL	: ground temperature in Kelvin
%					xmv	  : volumetric soil moisture
%					SAND	: SAND percentage 
%					CLAY	: caly percentage
%					ROB	  : dry matter density
%					epsol : dielectric constant
%
% Outputs 
%				  eg    : complex soil permittivity.
%
% Modello 
%SUBROUTINE DOBSOIL(xmv,eg)
function eg = MOD_PERMITTIVITY_DOBSOIL(xmv)

global f;
global Tsol;
global Rob;
global Sand;
global Clay;

% Dummy argument declaration
% 	REAL, INTENT(IN)	  :: xmv
% 	COMPLEX,INTENT(OUT)	:: eg
%
% Local variables declaration
%
% 	COMPLEX   ::  epsol,epfw,cx
% 	REAL	    ::  permit,alp,ros,epwi,salsol,ts,epx,epy,epsi,epss,xa,xb,f0
% 	REAL	    ::  epw0,toPI,seff,bet1,bet2,eps,x,y,fr

fr = f*1.e9;
% Air permitivitty
permit = 8.854E-12;
% Alfa (in refractive model, ulaby appendix E.62)
alp = 0.65;
% Density of the solid soil material
ros = 2.66; 
% e in the high freq limit (indep of salinity, Stogryn, 71, ulaby appendix E.112)
epwi = 4.9;
salsol = 0.65;
ts = Tsol - 273.157;
if (ts <= -0.5) 
  % SOIL   =   FROZEN SOIL	(Ulaby et al., 1986, p2101) 	
    epx = 5.;
    epy =  - 0.5;
    epsol = epx+1i*epy;
elseif (xmv < .02 && Sand >= .9) 
  % DRY SAND (ref Matzler, 1998)
    epsi = 2.53;
    epss = 2.79;
    xa = .002;
    f0 = .27e9;
    epsol = epsi + (epss - epsi)/(1 - 1i*fr/f0) + (0+1i*xa);
else
    %	SALINE WATER
    % ref: saline water (Ulaby/M/F p2024)
    %
    % E.24 ulaby 2024
    xa = 1 + 1.613E-5*ts*salsol - 3.656E-3*salsol + 3.21E-5*salsol*salsol - 4.232E-7*(salsol^3);
    %	E.23 dependence of e saline water
    epw0 = xa*(87.134 - 1.949E-1*ts - 1.276E-2*ts*ts + 2.491E-4*(ts^3)); 
    % E.26: correction to the free water relaxation time due to salinity
    xb = 1 + 2.282E-5*ts*salsol - 7.638E-4*salsol - 7.76E-6*salsol*salsol + 1.105E-8*(salsol^3);
    % E.25, E.17(2*PI*relaxation time of free water)
    toPI = xb*(1.1109E-10 - 3.824E-12*ts + 6.938E-14*ts*ts - 5.096E-16*(ts^3));
    % Dobson E.32 corrected, effective conductivity, function of soil texture
    seff =  - 1.645 + 1.939*Rob - 2.256*Sand + 1.594*Clay;
    % Dobson E.29
    cx = (epw0 - epwi)/(1.+1i*toPI*fr);
    % modified Debye equation, Dobson E.29		 
    epfw = epwi + cx - (0+1i*seff*(ros - Rob))/(2*pi*fr*permit*ros*xmv);
    % Dobson E.30 (related to real(e))
    bet1 = 1.275 - 0.519*Sand - 0.152*Clay;
    % Dobson E.31 (related to imag(e))
    bet2 = 1.338 - 0.603*Sand - 0.166*Clay;
    % Dobson E.22 soil permit
    eps = (1.01 + 0.44*ros)^2 - 0.062;
    % Ulaby E.111, Dobson E.28
    x = 1 + Rob*(eps^alp - 1)/ros + (xmv^bet1)*(real(epfw)^alp) - xmv;
    y = (xmv^bet2)*(abs(imag(epfw))^alp);
    epx = x^(1/alp);
    epy = y^(1/alp);
    % frozen soil
    epsol = (epx-1i*epy);
end
eg = epsol;

