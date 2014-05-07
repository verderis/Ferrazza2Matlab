
function [s_up,t_up] = MATRIX_DOUBLING_DOUBLE_VEGETATION(s_up,t_up,s_down,t_down)
%SUBROUTINE MATRIX_DOUBLING_DOUBLE_VEGETATION(s_up,t_up,s_down,t_down)

global nn;

% Dummy variables declaration
%REAL, DIMENSION(:,:), INTENT(INOUT) :: s_up,t_up
%REAL, DIMENSION(:,:), INTENT(IN)	  :: s_down,t_down
% Local variables declaration
%REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2,sv1,tv1
%REAL, DIMENSION(NN)			:: work
%INTEGER	:: info,i
% ami = zeros(nn,nn);am1 = zeros(nn,nn);
% am2 = zeros(nn,nn);sv1 = zeros(nn,nn);tv1 = zeros(nn,nn);

%REAL, DIMENSION(NN)			:: work
work = zeros (nn,1);


%id(:,:) = 0.
%FORALL(i = 1:NN) id(i,i) = 1.
id = eye(nn,nn);
am1 = id - s_up*s_down;
ami = am1;
%CALL RINV(NN,ami,NN,work,info)
ami = inv(ami);
%IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
am1 = t_up*s_down;
am2 = am1*ami;
sv1 = am2*t_up;
am1 = t_down*ami;
tv1 = am1*t_up;
s_up = sv1 + s_up;
t_up = tv1;