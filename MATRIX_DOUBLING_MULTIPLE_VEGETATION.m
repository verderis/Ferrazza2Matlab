% *************************************************************************************
% MOD_MATRIX_DOUBLING contains the algorithms that compute the matrix doubling among
% vegetation sublayers and between vegetation and soil. Both the procedure are called
% throught MATRIX_DOUBLING interface
%
% MATRIX_DOUBLING_VEGETATION
% MATRIX_DOUBLING_SOIL
% % *************************************************************************************
% MODULE MOD_MATRIX_DOUBLING
% 
% USE MOD_GLOBAL_PARAMETERS
% 
% IMPLICIT NONE
% 
% PRIVATE
% 
% INTERFACE MATRIX_DOUBLING
%   MODULE PROCEDURE MATRIX_DOUBLING_MULTIPLE_VEGETATION ,&
% 									 MATRIX_DOUBLING_DOUBLE_VEGETATION ,&
%                    MATRIX_DOUBLING_SOIL
% END INTERFACE 
% 
% PUBLIC MATRIX_DOUBLING 
% 
% CONTAINS

% SUBROUTINE MATRIX_DOUBLING_MULTIPLE_VEGETATION(sv,tv,N_iter)

function [sv, tv] = MATRIX_DOUBLING_MULTIPLE_VEGETATION(sv,tv,N_iter)

global nn;
	% Dummy variables declaration
	%REAL, DIMENSION(:,:), INTENT(INOUT) :: sv,tv
	%INTEGER, INTENT(IN)									:: N_iter
	% Local variables declaration
%REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2,sv1,tv1

% ami = zeros(nn,nn);am1 = zeros(nn,nn);
% am2 = zeros(nn,nn);sv1 = zeros(nn,nn);tv1 = zeros(nn,nn);

%REAL, DIMENSION(NN)			:: work
% work = zeros (nn,1);
%INTEGER	:: info,i

%id(:,:) = 0.
%FORALL(i = 1:NN) id(i,i) = 1.
id = eye(nn,nn);
for i = 1:N_iter
    am1 = id - (sv*sv);
    ami = am1;
    %CALL RINV(NN,ami,NN,work,info)
    ami = inv(ami);
    %IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
    
    am1 = tv*sv;
    am2 = am1*ami;
    sv1 = am2*tv;
    am1 = tv*ami;
    tv1 = am1*tv;
    sv = sv1 + sv;
    tv = tv1;
end

