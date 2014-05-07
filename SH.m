%REAL FUNCTION SH(teta,PI,ms,erfc)
function SH(teta,PI,ms,erfc)
% Dummy arguments declaration
%REAL, INTENT(IN)	:: teta, PI, ms, erfc 
% Local variable delaration
%REAL a0,a1,aa1,a2,aa3,a3,f0,f1,ff,mm2,mm

if (teta  ==  0.)
	sh = 1;
	return
end
a0 = sqrt(2/PI);
a1 = cos(teta)/sin(teta);
aa1 = a1*a1;
a2 = ms/a1;
mm = ms*ms;
mm2 = mm*2;
aa3 = aa1/mm2;
if (aa3 > 180.) 
	a3 = 0;
else
	a3 = EXP(-(aa1/mm2));
end
f0 = 0.5*(a0*a2*a3 - erfc);
f1 = 1 + f0;
ff = 1/f1;
sh = (1 - (0.5*erfc))*ff;