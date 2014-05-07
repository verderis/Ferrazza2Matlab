% ************ Integral of Karam function ************
%
%REAL FUNCTION INTKARAM(inf,sup,A,n,b0,bm)
function [INTKARAM] = INTKARAM(inf,sup,A,n,b0,bm)

% % Dummy variables declaration
% REAL	:: inf,sup,A,n,b0,bm
% % Local variables declaration
% REAL	:: b,passo,ck,arg
% INTEGER :: i,nint

nint = 10;
INTKARAM = 0.;
passo = (sup - inf)/nint;
nint = nint + 1;
b = inf;

for i = 1:nint
    if (i  ==  1 || i  ==  nint)
    	ck = 1.;
    else
 		ck = 2.;
    end

    arg = abs((b - bm)/(b0 - bm));
		% Essendo n un numero irrazionale il radicando (arg) deve essere positivo, inoltre per motivi di
		% precisione di calcolo quando l'argomento tende ad 1. il prodotto arg*PI/2 pu� eccedere PI/2. portando
		% il coseno ad un valore negativo e non regolare per un indice irrazionale.
	if (arg < 1.) 
        INTKARAM = INTKARAM + ck*A*cos((PI/2.)*arg)^n;	   
    end
    b = b + passo;
end
INTKARAM = INTKARAM*passo/2;

