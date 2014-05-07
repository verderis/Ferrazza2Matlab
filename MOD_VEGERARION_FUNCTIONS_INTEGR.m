%SUBROUTINE INTEGR(rm,iflag,sum)
function sum = MOD_VEGERARION_FUNCTIONS_INTEGR(rm, iflag)	
% CALLed by OSBORN
	%
	%  Elliptic integrals computation [9]
	%  iflag = 1 : first kind
	%  iflag =  - 1: second kind

	% Dummy variables declaration
%	REAL, INTENT(IN)    :: rm
%	REAL, INTENT(OUT)   :: sum
%	INTEGER, INTENT(IN) :: iflag

	% Local variables declaration
%	REAL				  :: denom,coeff,param
%	INTEGER				:: n,i,par,disp

sum = pi/2.;
denom = 1.;
n = 10;
for i = 1:n
	disp = 2*i - 1;
	par = 2*i;
	if (iflag == - 1) 
        denom = disp;
    end
		coeff = (ffatt(disp)/ffatt(par))*(ffatt(disp)/ffatt(par));
		param = rm^i;
		sum = sum + pi/2.*(iflag*coeff*param/denom);    % eq. 17.3.11 and 17.3.12 of [9]
end
