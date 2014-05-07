
%
%******************************************************************************
%
% *********** INPUT
%
% lai	    : LAI values																									[cm^2/cm^2] 
% species : (1)  Hardood,  (2)  Softwood
% name    : Kind of selected forest species
% dbh_min : Minimum value for tree diameter used for the distribution   [cm]
% dbh_man : Maximum value for tree diameter used for the distribution   [cm]
% dbh_mean: Mean value for tree diameter used for the distribution      [cm]
% std     : Gaussian Standard deviation
% thick   : Needles length or leaves thickness								          [cm]
% radius  : Needles or leaves radius																		[cm]
%
% *********** OUTPUT
%
% an0nn       : Leaves surface density					                        [cm^2/cm^2] 
% Trees_group	: Trees number for at a specific dbh									    [#/cm^2] 
% an0br				: Primary branches densities                              [cm^2/cm^2]
% an0br2			: Secondary branches densities                            [cm^2/cm^2]
% htr					: Trees height at a specific dbh						              [cm]
% dbh					: Distribution diameters trees														[cm]
% an0tr				: Trees densities									                        [#/cm^2]
%
%
%SUBROUTINE FOREST_GROWTH(lai,name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius,)

function [an0nn,Trees_group,htr,dbh,an0tr,Ncyl,an0br,an0br2,Dry_Litter_biomass] = FOREST_GROWTH(lai,name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius)


% Dummy argument declaration
%
% INPUT
%REAL, INTENT(IN)								:: lai,dbh_min,dbh_max,dbh_mean,std,thick,radius
%INTEGER, INTENT(IN)							:: name,species
%
% OUTPUT
%REAL, DIMENSION(:,:), POINTER	  :: an0br
%REAL, DIMENSION(:), POINTER		  :: htr,dbh,an0tr
%INTEGER, DIMENSION(:), POINTER	:: Ncyl

%REAL, INTENT(OUT)								:: an0nn,an0br2,Dry_Litter_biomass
%INTEGER, INTENT(OUT)						:: Trees_group

% Local variables declaration
%REAL, DIMENSION(:), ALLOCATABLE	:: Bio_tot,Bio_tr,Bio_f,Bio_b,Vol_tot,Vol_tr,Vol_f,Vol_b,dmax
%REAL, DIMENSION(:), ALLOCATABLE	:: Bio_F_tot,Bio_F_tr,Bio_F_f,Bio_F_b

%REAL, DIMENSION(5,2,2) :: biomassa
%REAL, DIMENSION(4,2,2) :: component
%REAL, DIMENSION(2,3)   :: A,n,bm,b0

%REAL		:: area,dbh_step,rhoDB,rhoDT,amoiB,amoiT,inf,sup,rho_b,rho_t,step,db,vol
%REAL		:: Trees_Number,rhoDF,amoiF,rho_f,secondary_branches,lbran,Foliage_biomass

%INTEGER	:: AllocErr,up_boundary

% Constants values for dry and fresh matter
dbh_step = 5.;
rhoDB = 0.4;
rhoDT = 0.4;
rhoDF = 0.3;
amoiB = 0.6;
amoiT = 0.5;
amoiF = 0.6;
%5., 0.4, 0.4, 0.3, 0.6, 0.5, 0.6/ 

% 4 Hardwood species
 % 5 Softwood species
 
%REAL, DIMENSION(5,2,2) :: biomassa
%REAL, DIMENSION(4,2,2) :: component
biomassa = [- 2.2094, - 1.9123,  - 2.48,   - 2.0127, 0.;...   
2.3867, 2.3651,  2.4835, 2.4342, 0.]';
biomassa(:,:,2) = [- 2.0336, - 2.2304,  - 2.5384, - 2.5356, - 2.0773;...
 2.2592, 2.4435,  2.4814, 2.4349, 2.3323]';

% Hardwood components coefficients 
% Softwood components coefficients 

component = [ - 4.0813, - 1.6911,  - 2.0129,  - 0.3065;...         
5.8816, 0.816,   - 1.6805,  - 5.424]';
component(:,:,2) = [ - 2.9584, - 1.5619,  - 2.098,   - 0.3737;...         
							4.4766, 0.6614,  - 1.1432,  - 1.8055]';

% Constants for Karam function, which determines the branches orientation
%REAL, DIMENSION(2,3)   :: A,n,bm,b0
%DATA A,n,bm,b0/2.42,2.6,2.12,2.3,1.98,2.01,&
%			        &3.66,5.84,3.29,3.99,2.67,3.1,&
%					 	  &.47,.39,.58,.58,.57,.6,&
%			        &1.,1.,0.,0.,0.,0./  

A = [2.42,2.6,2.12;2.3,1.98,2.01];
n = [3.66,5.84,3.29,3.99,2.67,3.1];
bm = [47,.39,.58,.58,.57,.6];
b0 = [1.,1.,0.,0.,0.,0.];


% Complete leaf area
area = 2*pi*radius*radius*thick*(1/radius + 1/thick);  

% projected leaf area
if (species  ==  1) 
	% Hardwood
  area = area/2;
else
  % Softwood
  area = area/2.5;
end

% Leaves density
an0nn = lai/area;

% Number of differents consideres diameters values
%Trees_group = NINT((dbh_max - dbh_min)/dbh_step)    
Trees_group = round((dbh_max - dbh_min)/dbh_step);

% Memory allocation for used pointers
%ALLOCATE(an0tr(Trees_group),Bio_tot(Trees_group),Bio_tr(Trees_group),Bio_f(Trees_group),Bio_b(Trees_group),Bio_F_tot(Trees_group),&
%				&Bio_F_tr(Trees_group),Bio_F_f(Trees_group),Bio_F_b(Trees_group),Vol_tot(Trees_group),Vol_tr(Trees_group),Vol_b(Trees_group),&
%				&Vol_f(Trees_group),htr(Trees_group),dbh(Trees_group),dmax(Trees_group),Ncyl(Trees_group),STAT = AllocErr)

%if(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria per an0tr'

% Diameter trunks step
dbh_step = (dbh_max - dbh_min)/Trees_group;

% Integration probability function at each interval

an0tr = zeros(Trees_group,1);
dbh = zeros(Trees_group,1);
Bio_tot = zeros(Trees_group,1);
Bio_tr = zeros(Trees_group,1);
Bio_f = zeros(Trees_group,1);

for i = 1:Trees_group  
	inf = dbh_min + (i - 1) * dbh_step;
	sup = inf + dbh_step;
    % [INTGAUSS] = INTGAUSS(inf,sup,mean,std)
	an0tr(i) = INTGAUSS(inf,sup,dbh_mean,std);
	dbh(i) = inf + dbh_step / 2;
	%
	% All the weight are expressed in grams
	%
	% Total dry biomass
	Bio_tot(i) = 1000 * exp(biomassa(name,1,species) + biomassa(name,2,species) * log(dbh(i)));
  % Trunks dry biomass 
	Bio_tr(i) = Bio_tot(i) * (exp(component(3,1,species) + component(3,2,species) / dbh(i)) + exp(component(4,1,species) + component(4,2,species) / dbh(i)));
  % Foliage dry biomass
	Bio_f(i) = Bio_tot(i) * exp(component(1,1,species) + component(1,2,species) / dbh(i));
end  

% Branches dry biomass
Bio_b = Bio_tot - Bio_tr - Bio_f;

% Number of tree at each interval
an0tr = an0tr + (1 - sum(an0tr))/Trees_group;

% Fresh biomass
Bio_F_tr = Bio_tr ./ (1 - amoiT);
Bio_F_b = Bio_b ./ (1 - amoiB);
Bio_F_f = Bio_f ./ (1 - amoiF);
Bio_F_tot = Bio_F_tr + Bio_F_b + Bio_F_f;

% Efficacious density of fresh matter
rho_b = rhoDB / ((1 - amoiB) + rhoDB*amoiB);
rho_t = rhoDT / ((1 - amoiT) + rhoDT*amoiT); 	
rho_f = rhoDF / ((1 - amoiF) + rhoDF*amoiF); 	

% Volumes of each dielectric component
Vol_tr = Bio_F_tr ./ rho_t;
Vol_b = Bio_F_b ./ rho_b;
Vol_f = Bio_F_f ./ rho_f;
Vol_tot = Vol_tr + Vol_b + Vol_f;

% Trunks height
htr = 4 * Vol_tr ./ (pi * dbh.^2);

% Foliage biomass link with LAI values
if (species  ==  1)
	Foliage_biomass = lai / 1.4889;                        % [t/Ha]
    Dry_Litter_biomass = .01 *(1.2463 * Foliage_biomass);  % [g/cm^2]
else
    Foliage_biomass = lai / 0.4314;                        % [t/Ha]
    Dry_Litter_biomass = .01 *(0.2498 * Foliage_biomass);  % [g/cm^2]
end

% Total number trees
Trees_Number = 1.e-2 * Foliage_biomass / (sum(Bio_f .* an0tr)); % [#/cm^2]

% Trees number for each selected dbh value
an0tr = an0tr * Trees_Number;

% Branches maximum value
if (species  ==  1)
  dmax = dbh/4;
else
	dmax = 10^( -0.468 + 0.803 * log10(1.0 + 1.2*dbh));
end

% Print of parameter values
%OPEN(192,file = "Biomasse.txt",FORM = 'FORMatted',status = 'unknown')
%save

%if(.true.) THEN
aux = [];

for i = 1:Trees_group
    aux = [aux;lai,dbh(i),htr(i)/1e2,round(an0tr(i)*1e8),Bio_tot(i)*...
    an0tr(i) *1e2,Bio_tr(i)* an0tr(i) *1e2,Bio_f(i)* an0tr(i) *...
    1e2,Bio_b(i)* an0tr(i) *1e2,Bio_F_tot(i)* an0tr(i)*... 
    1e2,Bio_F_tr(i)* an0tr(i) *1e2,Bio_F_f(i)* an0tr(i)*...
    1e2,Bio_F_b(i)* an0tr(i) *1e2,Vol_tot(i)* an0tr(i) *1e2,Vol_tr(i)...
    *an0tr(i) *1e2,Vol_f(i)* an0tr(i) *1e2,Vol_b(i)* an0tr(i) *1e2];
end

  %WRITE(192,*)
  %WRITE(192,*)

%11 FORMAT(F4.1,1X,2(F6.2,1X),I3,12(F9.3))    
  
%  if (lai  ==  10) CLOSE(192)

%ENDIF

step = 0.2;

% "secondary_branches * 100" is the percentage component of secondary branches presence
secondary_branches = .3;

% Secondary branches density
an0br2 = secondary_branches * sum(Vol_b .* an0tr) / (pi * 0.3 * 0.3 * 20);

Ncyl = round(dmax/step);
%up_boundary = MAXVAL(Ncyl(:))
up_boundary = max(Ncyl);

%ALLOCATE(an0br(Trees_group,up_boundary),STAT = AllocErr)
%if(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria'
an0br(:,:) = 0.; 
%
% Primary branches densities for each dbh value
%
for i = 1:Trees_group
	inf = 0.;
	sup = step;
	db = (sup + inf)/2.;					% dimameter value
    if(lb/db > 50)
        lbran = db*50.;
    else
        lbran = lb;
    end
	vol = pi*db*db*lbran/4.;   % volume cylinder
	for j = 1:Ncyl(i)	    	
		% function [INTKARAM] = INTKARAM(inf,sup,A,n,b0,bm)
        %an0br(i,j) = INTKARAM(inf,sup,A(1,3)/dmax(i),n(1,3),b0(1,3)*...
        %    dmax(i),bm(1,3)*dmax(i))* (1 - secondary_branches) * Vol_b(i)...
        %    * an0tr(i) / vol;	% [#/cm^2]      
        
        
        an0br(i,j) = INTKARAM(inf,sup,A(1,3)/dmax(i),n(1,3),b0(1,3)...
            *dmax(i),bm(1,3)*dmax(i)) *(1 - secondary_branches) * Vol_b(i)...
            * an0tr(i) / vol;
        
        inf = sup;
        sup = sup + step;
        db = db + step;
    if (lb/db > 50)
        lbran = db*50.;
    else
        lbran = lb;
    end
    vol = pi*db*db*lbran/4.;
    end
end

%DEALLOCATE(Bio_tot,Bio_tr,Bio_f,Bio_b,Vol_tot,Vol_tr,Vol_b,Vol_f,dmax,Bio_F_tot,Bio_F_tr,Bio_F_f,Bio_F_b)

%RETURN