%% pruebas
clc;
%clear all
MOD_GLOBAL_PARAMETERS

%etc = MOD_SOIL_FUNCTIONS_LITTER(0.2,0.3,1);                            OK!
%[Sgg,rr] = MOD_SOIL_FUNCTIONS_IEM(0.01,0.1,0.2,0.2,1);                 OK!
%[s,ke] = MOD_VEGERARION_FUNCTIONS_RAYLEIGH_GANS(1,1,0.02,0.3,0.4);  	OK!
%[s,ke] = MOD_VEGERARION_FUNCTIONS_PHYSICAL_OPTICS(1,1,0.02,0.3,0.4)    OK!

% CALL INFINITE_LENGTH(1.0,100.,0.5,0.3,s,ke,'S')
tic
[s,ke] = MOD_VEGERARION_FUNCTIONS_INFINITE_LENGTH(1.0,100,0.5,0.3,'s');
toc
a = 1;


