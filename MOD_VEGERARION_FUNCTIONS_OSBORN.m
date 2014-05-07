%******************************************************************
%SUBROUTINE OSBORN(a,b,c,rl,rm,rn)
function [rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_OSBORN(a,b,c)
% CALLed by RAYLEIGH_GANS
% 2.23 - 2.25 of [3] (demagnetization factors)

% Dumyy variables declaration
%REAL, INTENT(IN)    :: a,b,c
%REAL, INTENT(OUT)   :: rl,rm,rn

% Local variables declaration
%REAL				        :: rapba,rapca,e,e2,e1,ellk,elle

if ((a == b) && (a > c)) 
    %CALL DISC(a,c,rl,rm,rn)
    [rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_DISC(a,c);
    return
elseif ((a == b) && (c > a)) 
    %CALL NEEDLE(a,c,rl,rm,rn)
    [rl,rm,rn] = MOD_VEGERARION_FUNCTIONS_NEEDLE(a,c);
    return
else
    rapba = b/a;
    rapca = c/a;
    e = sqrt(1. - (rapba*rapba));
    e2 = e*e;
    e1 = 1. - e2;
    %CALL INTEGR(e2,1,ellk)
    ellk = MOD_VEGERARION_FUNCTIONS_INTEGR(e2,1);
    %CALL INTEGR(e2,-1,elle)
    elle = MOD_VEGERARION_FUNCTIONS_INTEGR(e2,-1);
    rl = rapca*sqrt(e1)*(ellk - elle)/(e2);
    rm = rapca*(elle-e1*ellk)/(e2*sqrt(e1));
    rn = 1. - rapca*elle/sqrt(e1);
end
	%******************************************************************