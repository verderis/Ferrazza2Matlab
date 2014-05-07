%SUBROUTINE DISC(a,c,rl,rm,rn)
	
function [rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_DISC(a,c)
% CALLed by OSBORN

% Dummy variables declaration
%REAL, INTENT(IN)  :: a,c
%REAL, INTENT(OUT) :: rm,rl,rn

% Local variables declaration
%REAL			        :: aa,cc,diff,rad,atg,ac3

aa = a*a;
cc = c*c;
diff = aa - cc;
rad = sqrt(diff);
atg = atan(sqrt(cc/diff));
ac3 = (diff)^1.5;
rl = ( -c*rad/aa + pi/2. - atg)/ac3;
rl = aa*c/2.*rl;
rm = rl;
rn = 2.*(rad/c - pi/2. + atg)/ac3;
rn = aa*c/2.*rn;