%SUBROUTINE MATRIX_DOUBLING_SOIL(sv_star,tv_star,tv,sg0,svg1)
function svg1 = MATRIX_DOUBLING_SOIL(sv_star,tv_star,tv,sg0)

% Dummy variables declaration
%REAL, DIMENSION(:,:), INTENT(INOUT) :: svg1
%REAL, DIMENSION(:,:), INTENT(IN)	  :: sv_star,tv_star,tv,sg0

% Local variables declaration
%REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2
%REAL, DIMENSION(NN)			:: work
%INTEGER	:: info,i

%    id(:,:) = 0.
%    FORALL(i = 1:NN) id(i,i) = 1.
id = eye(nn,nn);
am1 = id - sv_star*sg0;
ami = am1;
%CALL RINV(NN,ami,NN,work,info)
ami = inv(ami);
%IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
am1 = sg0*ami;
am2 = tv_star*am1;
svg1 = am2*tv;            % Interaction between vegetation and soil