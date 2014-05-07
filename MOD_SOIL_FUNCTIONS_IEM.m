%SUBROUTINE IEM(sigmaz,lcorr,vmoist,Sgg,rr,litter_gmoist,biomass)
function [Sgg,rr] = MOD_SOIL_FUNCTIONS_IEM(sigmaz,lcorr,vmoist,litter_gmoist,biomass)

global nm;
global nm1;
global nij;
global f;
global chi;
global om;
global CORRELATION

% Dummy arguments declaration 
%real, INTENT(IN)														  :: sigmaz,vmoist,lcorr
%real, DIMENSION(:,:,:,:,:),INTENT(OUT)			  :: Sgg
Sgg = zeros(2,2,nm1,2,2);

%real, OPTIONAL, DIMENSION(:,:),INTENT(OUT)  	:: RR
rr = zeros(nij,2);

%real, OPTIONAL, INTENT(IN)                    :: litter_gmoist, biomass
% Local variable delaration
%COMPLEX, DIMENSION(nm)	:: GX

%COMPLEX					:: eta2,ff,rhh0,rvv0,egc,fvv,fhh,f1vv,f1hh,fhv,fvh,f1hv,f1vh,f1hhs,f1vvs,f1hvs,f1vhs
%COMPLEX					:: tn,tnm,th,thm,tp,tm,r,sq,sqs,ihhn,ihvn,ivhn,ivvn,rpa,ror,ihh2,ivv2,ihh1,ivv1
%COMPLEX					:: r21,r22,rt,rii,ctet2,ivh2,ihv2,ihh3,ivv3,ivh3,ihv3,ivh1,ihv1

%real(8)					:: Sommahh,Sommavv,Sommahv,Sommavh

%real, DIMENSION(nm,2,2)	:: Gtoro
Gtoro = zeros(nm,2,2);
%real, DIMENSION(2,2)		:: Gdirez
Gdirez = zeros(2,2);
%real, DIMENSION(64)			:: FFss
FFss = zeros(64);

%real		:: kx,ky,kz,ksx,ksy,ksz,kdx,kdy,wn,lambda,k,kk,llcorr,fi,NM1,conv,delfis,sqrt2,ak1k1,eta1
%real		:: alfag,betag,ppg,qqg1,cgg,tetaj,tetais,teta,steta
%real		:: cteta,abscteta,steta2,cteta2,ssteta,ccteta,tetagraz,ak1steta,qqg,tet2,c1,ak1,stet2,r12,beckspiz1
%real		:: beckspiz2,beckspiz,tteta,tteta2,r11,tetas,stetas,ctetas,sstetas,cctetas,absctetas,o0,cdirez,cdirezg
%real		:: fis,sfis,cfis,s,ss,cs,css,sf,sf1,csf,c2,c1s,c2s,sigma1,a,b,sigmahv,sigmavh,sigmavv,sigmahh,WN1,WN2
%real		:: deltet

%INTEGER	:: nsup,k1,k2,k3,n,i,j,istoki,istoks

%LOGICAL :: RR_present,Litter_present

%if(PRESENT(RR)) THEN
%if nargin>1
	RR_present = true;
%else
%	RR_present = false;
%end

Litter_present = false;
if(nargin>3) 
    if(biomass > 0.) 
      Litter_present = true;
    end
end

fi = 0;
%NM1 = nm/2 + 1; %ojo aca! la global es pisada por NM1
conv = 180/pi;
deltet = pi./(2*nij);
delfis = 2*pi./nm;
sqrt2 = sqrt(2.);
c1 = 2./sqrt(pi);
nsup = 32;
lambda = 30/f;
k = 2*pi/lambda;
kk = k^2;
ak1 = (2*pi)/lambda ;
ak1k1 = ak1*ak1;
eta1 = 1.;

if(Litter_present)
  %CALL LITTER(vmoist,litter_gmoist,biomass,egc) 
  egc = MOD_SOIL_FUNCTIONS_LITTER(vmoist,litter_gmoist,biomass);
else
  %CALL DOBSOIL(vmoist,egc) 
  egc = MOD_PERMITTIVITY_DOBSOIL(vmoist);
end
eta2 = sqrt(1/egc);
alfag = k*imag(sqrt(egc));
betag = k*real(sqrt(egc));
ppg = 2*alfag*betag;
qqg1 = betag*betag - alfag*alfag;
cgg = abs((eta1 - eta2)/(eta1 + eta2));
cgg = cgg*cgg;

llcorr = lcorr*lcorr;
for j = 1:nij           % Incidence loop, toroi: 
	tetaj = (j - 1)*deltet;
	if (RR_present) 
        rr(j,:) = 0.;
	end

	for i = 1:nij         % Scattering loop, toros: 
		tetais = (i - 1)*deltet;
		Gtoro(:,:,:) = 0.;

		for k1 = 1:2     % integrale (gauss) in teta, tetaloop: 
			%teta = .5*(deltet*Chi(k1) + 2*tetaj + deltet);
            teta = .5*(deltet*chi(k1) + 2*tetaj + deltet);
			steta = sin(teta);
			cteta = cos(teta);
			abscteta = cos(teta);
			steta2 = steta*steta;
			cteta2 = cteta*cteta;
			ssteta = steta2;
			ccteta = cteta2;
			if(i == 1) 
				tetagraz = pi/2. - teta;
				ak1steta = ak1*steta;
				qqg = qqg1 - ak1steta*ak1steta;
				tet2 = atan(ak1steta*sqrt2/sqrt(sqrt(ppg*ppg + qqg*qqg) + qqg));
				ctet2 = cos(tet2);
				ctet2 = sqrt(1 - (steta*steta)/egc);
				stet2 = sin(tet2);
				r11 = eta1*abscteta;
				r12 = eta1*ctet2;
				r21 = eta2*abscteta;
				r22 = eta2*ctet2;
				rt = (r21 - r12)/(r21 + r12);
				rii = (r11 - r22)/(r11 + r22);
				beckspiz1 = 4*pi*sigmaz*sin(tetagraz)/lambda;
				beckspiz2 =  beckspiz1*beckspiz1;
				beckspiz = exp( - beckspiz2);
				if(RR_present) 
					rr(j,1) = rr(j,1) + (om(k1)/2)*beckspiz*rii*conj(rii);
					rr(j,2) = rr(j,2) + (om(k1)/2)*beckspiz*rt*conj(rt);
				end
				tteta = steta/abscteta;
				tteta2 = tteta*tteta;
			end
			ff = sqrt(egc - ssteta);

			for k2 = 1:2       % integrale (gauss) in tetas, tetasloop: 
				tetas = .5*(deltet*chi(k2) + 2*tetais + deltet);
				stetas = sin(tetas);
				ctetas = cos(tetas);
				sstetas = stetas*stetas;
				cctetas = ctetas*ctetas;
				absctetas = cos(tetas);
				o0 = om(k1)*om(k2)/4.;
				cdirez = (deltet*steta/absctetas)*o0;
				cdirezg = cdirez/(4.*pi);

				for k3 = 1:nm1      % calcoli per 0<fis<180, fisloop: 
					fis = (k3 - 1)*delfis;
					sfis = sin(fis);
					cfis = cos(fis);
					% ^^^^^^^^^^^^^^^^^^^^^^^*
					% Soil bistatic scattering coefficient
					% ^^^^^^^^^^^^^^^^^^^^^^^*
					kx = k*sin(teta)*cos(fi);
					ky = k*sin(teta)*sin(fi);
					kz = k*cos(teta);
					ksx = k*sin(tetas)*cos(fis);
					ksy = k*sin(tetas)*sin(fis);
					ksz = k*cos(tetas);
					kdx = ksx - kx;
					kdy = ksy - ky;
					% ^^^^^^ pag 206,211 fung*^^^^*
					s = sin(teta);
					ss = sin(tetas);
					cs = cos(teta);
					css = cos(tetas);
					sf = sin(fis - fi);
					sf1 = sin(fi - fis);
					csf = cos(fis - fi);
					sq = sqrt(egc - (s^2));
					c1 = (csf - s*ss)/(sq*css);
					c2 = s*(ss - s*csf)/css;
					sqs = sqrt(egc - (ss^2));
					c1s = (csf - s*ss)/(sqs*cs);
					c2s = ss*(s - ss*csf)/cs;
					% ^^^^^ coeff. di fresnel ^^^^*
					rhh0 = (cs - sq)/(cs + sq);
					rvv0 = (egc*cs - sq)/(egc*cs + sq);
					r = (rvv0 - rhh0)/2;
					tn = 1 + rvv0;
					tnm = 1 - rvv0;
					th = 1 + rhh0;
					thm = 1 - rhh0;
					tp = 1 + r;
					tm = 1 - r;
					% ^^^^ kirchhoff field coefficent ^^^^
					rpa = rvv0;
					ror = rhh0;
					fvv = ((2*rpa)/(cs + css))*(s*ss - (1 + cs*css)*csf);
					fhh =  - ((2*ror)/(cs + css))*(s*ss - (1 + cs*css)*csf);
					fhv = 2*r*sf;
					fvh = 2*r*sf1;
					% ^^^^ complementary field coefficent ^^^^
					f1hhs = (css*thm - sqs*th)*(th*csf + thm*c1s) - (th^2 - (th*thm*css)/sqs)*c2s;
					f1hh = (cs*thm - sq*th)*(th*csf + thm*c1) - (thm^2 - (th*thm*cs)/sq)*c2;
					f1vv =  - (cs*tnm - sq*tn/egc)*(tn*csf + tnm*egc*c1) + (tnm^2 - (tn*tnm*cs)/sq)*c2;
					f1vvs =  - (css*tnm - sqs*tn/egc)*(tn*csf + tnm*egc*c1s) + (tn^2 - (tn*tnm*css)/sqs)*c2s;
					f1vhs =  - (css*tm - sqs*tp/egc)*((tp/cs) + tm*egc/sqs)*sf - (tp^2 - (tp*tm*css)/sqs)*(ss^2)*sf;
					f1hvs =  - (css*tp - sqs*tm)*(tm/cs + tp/sqs)*sf - (tm^2 - (tp*tm*css)/sqs)*(ss^2)*sf;
					f1vh = (cs*tp - sq*tm)*((tm/css) + tp/sq)*sf + ((tp^2) - (tp*tm*cs)/sq)*(s^2)*sf;
					f1hv = (cs*tm - sq*tp/egc)*((tp/css) + (tm*egc)/sq)*sf + (tm^2 - tp*tm*cs/sq)*(s^2)*sf;
					% ^^^* scattering coefficient computation ^^^^

					sigma1 = 0.5*(k^2)*exp( - (sigmaz^2)*(kz*kz + ksz*ksz));

					Sommahh = 0.;
					Sommahv = 0.;
					Sommavh = 0.;
					Sommavv = 0.;

					a = (ksz + kz);
					b = exp( - (sigmaz*sigmaz)*kz*ksz);
					ihh1 = fhh*b;
					ihv1 = fhv*b;
					ivh1 = fvh*b;
					ivv1 = fvv*b;
					ihh2 = f1hh;
					ihv2 = f1hv;
					ivh2 = f1vh;
					ivv2 = f1vv;
					ihh3 = f1hhs;
					ihv3 = f1hvs;
					ivh3 = f1vhs;
					ivv3 = f1vvs;
					for n = 1:nsup %nloop: 
						ihh1 = ihh1*a;
						ihv1 = ihv1*a;
						ivh1 = ivh1*a;
						ivv1 = ivv1*a;
						ihh2 = ihh2*ksz;
						ihv2 = ihv2*ksz;
						ivh2 = ivh2*ksz;
						ivv2 = ivv2*ksz;
						ihh3 = ihh3*kz;
						ihv3 = ihv3*kz;
						ivh3 = ivh3*kz;
						ivv3 = ivv3*kz;

						ihhn = ihh1 + 0.5*(ihh2 + ihh3);
						ihvn = ihv1 + 0.5*(ihv2 + ihv3);
						ivhn = ivh1 + 0.5*(ivh2 + ivh3);
						ivvn = ivv1 + 0.5*(ivv2 + ivv3);

						% Roughness spectrum
						if(CORRELATION ==  1)
							% 1) Exponential
							wn1 = llcorr/(n*n);
							wn2 = (1. + (llcorr*(kdx*kdx + kdy*kdy))/(n*n))^( - 1.5);
							wn = wn1*wn2;
						else
							% 2) Gaussian
							wn = (llcorr/(2*n))*exp( - llcorr*((kdx*kdx + kdy*kdy)/(4*n)));
						end

						if(n  ==  1)
							FFss(n) = sigmaz^2;
						else
							FFss(n) = FFss(n - 1)*sigmaz^2/n;
						end
						%
						Sommahh = Sommahh + wn*abs(ihhn)*FFss(n)*abs(ihhn);
						Sommahv = Sommahv + wn*abs(ihvn)*FFss(n)*abs(ihvn);
						Sommavh = Sommavh + wn*abs(ivhn)*FFss(n)*abs(ivhn);
						Sommavv = Sommavv + wn*abs(ivvn)*FFss(n)*abs(ivvn);

                    end

					sigmahh = (sigma1*Sommahh);
					sigmahv = (sigma1*Sommahv);
					sigmavh = (sigma1*Sommavh);
					sigmavv = (sigma1*Sommavv);

					Gdirez(1,1) = sigmavv;
					Gdirez(1,2) = sigmavh;
					Gdirez(2,1) = sigmahv;
					Gdirez(2,2) = sigmahh;
					%
					%  ^^^^^^^^^^^^
					%  * Gdirez(1,1) = sigma0vv *
					%  * Gdirez(1,2) = sigma0vh *
					%  * Gdirez(2,1) = sigma0hv *
					%  * Gdirez(2,2) = sigma0hh *
					%  ^^^^^^^^^^^^
					%
                    
					% Gtoro(k3,:,:) = Gtoro(k3,:,:) + cdirezg*Gdirez(:,:); %original
                    Gtoro(k3,:,:) = Gtoro(k3,:,:) + reshape(cdirezg*Gdirez(:,:), 1, 2, 2);
                    
					% It is enough to compute the averaged value in fis, either incident radiation and scattered one
					% are uniformly distribuited in the space
					%
                end
            end
        end
		%
		% contributions for 180<fis<360
		%Gtoro(nm/2 + 2:nm,:,:) = Gtoro(nm/2:2: - 1,:,:);
        Gtoro(nm/2 + 2:nm,:,:) = Gtoro(nm/2:-1:2,:,:); %se da vuelta porque en fortran es start, stop, step y en matlab start, step, stop
		%
        for istoki = 1:2
            for istoks = 1:2
				%GX(:) = CMPLX(Gtoro(:,istoks,istoki),0.);
                gx = Gtoro(:,istoks,istoki);
				%CALL CFSTFT(NF,gx)
                Gx = fft(gx);
				% normalization
				Gx(1) = Gx(1) / nm;
				%Gx(2:) = Gx(2:)*2/nm
                Gx(2:end) = Gx(2:end)*2/nm;
				Sgg(i,j,1:nm1,istoks,istoki) = real(Gx(1:nm1));
            end     % istoks
        end     % istoki
    end
end
% RETURN
% END SUBROUTINE IEM
%
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*