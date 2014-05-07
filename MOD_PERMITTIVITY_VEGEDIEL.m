% SUBROUTINE VEGEDIEL
% CALLED BY: RAYGANS, PHO, CYL, CORNER
%
% Computes vegetation permittivity uSINg:
%   -   empirical fitting of experimental data in [10] for volumetric moisture <= 50%;
%   -   The model of [11] for volumetric moisture > 50%.
%
% Input  : gravimetric moisture (amoi)
% output : real (er) and imaginary (ej) parts of permittivity.
%
function [er, ej] = MOD_PERMITTIVITY_VEGEDIEL(dsdw,amoi,select)

global f;
% Local variables
%REAL	  :: eps0,fhz,t,sw,eoo,eo,atsw,eos,btsw,tau2PI,tau2PIs,dd,ffi,ss25,ss,erw,ejw,amd,zz,vw,er0
 
 eps0 = 8.854E-12; % vacuum absolute permittivity
 fhz = f*1.E+09;   % frequency in Hz
 t = 20.;          % tyPIcal vegetation temperature (in C)
 sw = 10.;         % tyPIcal salinity (parts per thousand)
 %
 % Liquid water permittivity by [12]
 %
 eoo = 4.9;																                                          % e16
 eo = 87.134 - .1949*t - .01276*t*t + .0002491*t*t*t;							                  % e23
 atsw = 1 + t*sw*1.613E-05 - sw*3.656E-03 + sw*sw*3.21E-05 - sw*sw*sw*4.232E-07;     % e24
 eos = eo*atsw;															                                        % e22
 btsw = 1 + t*sw*2.282E-05 - sw*7.638E-04 - sw*sw*7.76E-06 + sw*sw*sw*1.105E-08;	    % e26
 tau2PI = 1.1109E-10 - t*3.824E-12 + t*t*6.938E-14 - t*t*t*5.096E-16;			          % e17
 tau2PIs = tau2PI*btsw;													                                    % e25
 dd = 25 - t;
 ffi = dd*(2.033E-02 + dd*1.266E-04 + dd*dd*2.464E-06 - sw*(1.849E-05 - dd*2.551E-07 + dd*dd*2.551E-08)); % e28b
 ss25 = sw*(0.18252 - sw*1.4619E-03 + sw*sw*2.093E-05 - sw*sw*sw*1.282E-07);		      % e28a
 ss = ss25*exp( -ffi);														                                    % e27
 erw = eoo + (eos - eoo)/(1 + (fhz*tau2PIs)^2);									                % liquid water pemittivity (real)
 ejw = fhz*tau2PIs*(eos - eoo)/(1 + (fhz*tau2PIs)^2) + ss/(2*pi*fhz*eps0);		  % liquid water permittivity (imaginary)
 %
 % The model of [11] is used for higher moistures.
 %
 if (select == 1) 
 % Matzler
       amd = 1. - amoi;
       zz = 0.522*(1 - 1.32*amd);
       er = zz*erw + 0.51 + 3.84*amd;
       ej = zz*ejw;
     %
 % For lower moistures an emPIrical formula is used, obtained by
 % fitting data published in [10] and impoSINg the liquid water
 % permittivity to be obtained when the moisture is 100%.
 %
 else
 % Ulaby
       vw = amoi*dsdw/(1 - amoi*(1 - dsdw));  % volumetric moisture
       er0 = 1.75;
       er = er0 + (erw - er0)*(1 - exp(-0.55*vw/(1. - vw)));
       ej = ejw*(1. - exp(-0.4*vw/(1. - vw))); 
 end
