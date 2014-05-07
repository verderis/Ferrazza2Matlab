%REAL FUNCTION ERF(teta,c1,ms)
function ERF(teta,c1,ms)

% Dummy arguments declaration
% REAL, INTENT(IN)	:: teta, c1, ms
% Local variable delaration
% REAL		x,xx,x1,x2,z,d,sum,y
% INTEGER i
if (teta == 0) 
	erf = 1.;
	return
end
x1 = (cos(teta)/sin(teta));
x2 = sqrt(2.)*ms;
z = x1/x2;
d = z/100;
sum = 0;
for i = 1:100    
	x = (2*i - 1)*(z/200);
	xx = x*x;
	if (xx > 180.)
		y = 0.;
    else
		y = EXP(-xx);
    end
	sum = sum + y;
end
erf = c1*d*sum;

END FUNCTION ERF
%
