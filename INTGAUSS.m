%
% ************ Integral of Gauss function ************
%
%REAL FUNCTION INTGAUSS(inf,sup,mean,std)

function [INTGAUSS] = INTGAUSS(inf,sup,mean,std)

% Dummy variables declaration
%REAL	:: inf,sup,mean,std
% Local variables declaration
%REAL	:: b,passo,ck
%INTEGER :: i,nint

nint = 10;

INTGAUSS = 0.;
passo = (sup - inf)/nint;
nint = nint + 1;
b = inf;

for i = 1:nint
    if (i  ==  1 || i  ==  nint)
		ck = 1.;
    else
		ck = 2.;
    end

	INTGAUSS = INTGAUSS + ck * exp(( - (b - mean)^2)/(2*std*std));
	b = b + passo;
end

INTGAUSS = (1/(std*sqrt(2*pi))) * INTGAUSS * passo/2.;
