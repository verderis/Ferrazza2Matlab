	% CALLed by OSBORN
    %SUBROUTINE NEEDLE(a,c,rl,rm,rn)
    
function [rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_NEEDLE(a,c)
    
% Dummy variables declaration
%REAL, INTENT(IN)  ::a,c
%REAL, INTENT(OUT) :: rl,rm,rn

% Local variables declaration
%REAL			        :: aa,cc,ac1,ac2,ac3,aa0,aa1

aa = a*a;
cc = c*c;
ac1 = sqrt(cc - aa);
%ac2 = ALOG((c - ac1)/(c + ac1));
ac2 = log((c - ac1)/(c + ac1));
ac3 = 1/(cc - aa);
aa0 = ac3*(c/aa + (0.5/ac1)*ac2);
aa1 = -ac3*((1/ac1)*ac2 + 2./a);
rl = aa*c/2.*aa0;
rm = rl;
rn = aa*c/2.*aa1;