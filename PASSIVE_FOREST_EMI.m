% PROGRAM Emiforest
% 
MOD_GLOBAL_PARAMETERS

global nn;
global DSDWW;
global CORRELATION_LENGTH;
global SOIL_ROUGHNESS;
global VOLMOIST_TRUNKS;
global nij;
% USE MOD_GLOBAL_PARAMETERS
% USE MOD_MATRIX_DOUBLING
% USE MOD_SOIL_FUNCTIONS
% USE MOD_VEGETATION_FUNCTIONS
% 
% IMPLICIT NONE

% Variables to describe forest architecture
% REAL, DIMENSION(:,:,:), ALLOCATABLE	:: Sigmaatr
% REAL, DIMENSION(:,:), POINTER				:: an0br
% REAL, DIMENSION(:), POINTER					:: htr,dbh,an0tr
% 
% INTEGER, DIMENSION(:), POINTER			:: Ncyl
% 
% REAL, DIMENSION(NN/2,NN/2,NM1,2,2)	:: Sgg
% 
% REAL, DIMENSION(NN/2,NN/2,2,2,2)		:: Sb,Sb0,Sb20,Snn,Snn0
% 
% REAL, DIMENSION(NN/2,2)							:: Emi,Emig,Emig1,Emitr,Emiv,keb,keb0,keb20,Kenn,Kenn0,RR
Emi = zeros (nn/2,2);
Emig = zeros (nn/2,2);
Emig1 = zeros (nn/2,2);
Emitr = zeros (nn/2,2);
Emiv = zeros (nn/2,2);
keb = zeros (nn/2,2);
keb0 = zeros (nn/2,2);
keb20 = zeros (nn/2,2);
Kenn = zeros (nn/2,2);
Kenn0 = zeros (nn/2,2);
RR = zeros (nn/2,2);

% REAL, DIMENSION(NN,NN)							:: Sv,Tv,Sg0,Stot,Id,Svg1
Sv = zeros (nn,nn);
Tv = zeros (nn,nn);
Sg0 = zeros (nn,nn);
Stot = zeros (nn,nn);
Id = zeros (nn,nn);
Svg1 = zeros (nn,nn);
% 
% % Variables which contain transmissivity values
% REAL, DIMENSION(NN)									:: Abstr,EN,EEN,EB,EEB,EC,EEC,Cot(NN/2)
Abstr = zeros(nn,1);
EN = zeros(nn,1);
EEN = zeros(nn,1);
EB = zeros(nn,1);
EEB = zeros(nn,1);
EC = zeros(nn,1);
EEC = zeros(nn,1);
% 
% % VGeneric variables
% REAL					  :: akznn,akzb,akatr,amvcil,vmoist,deltet,theta,an0br2,dsdw_trunks,litter_gmoist,Fresh_Litter_biomass
% REAL					  :: an0nn,fm,lcorr,sigmaz,cc,tau,dbh_min,dbh_max,dbh_mean,std,thick,radius,Dry_Litter_biomass
% 
% INTEGER					:: i,j,ja,ia,ivmoist,itau,ilai,istoki,istoks
% INTEGER					:: name,species,Ntr,itr,LAI_start,allocerr,ibiomass
            
% Path for scattering and extinction leaves matrix



PathF = 'C:\Fran\Varios\PassiveModel\MATRICI_FORESTE\MATRICI_DI_SCATTERING\FOGLIE_C\';

%CHARACTER(LEN = *), PARAMETER	:: PathF = "C:\Documents and Settings\andrea\Documenti\Universita\Dottorato\Codice\&
%                                          Matrici Foreste\Matrici di scattering\Foglie L\'


% Path for scattering and extinction branches matrix

PathB = 'C:\Fran\Varios\PassiveModel\MATRICI_FORESTE\MATRICI_DI_SCATTERING\&LEN_25\FUNZIONE_PESO_C\';

%CHARACTER(LEN = *), PARAMETER	:: PathB = "C:\Documents and Settings\andrea\Documenti\Universita\Dottorato\Codice\&
%                                          Matrici Foreste\Matrici di scattering\Len 25\Funzione Peso L\'

Chacil = ['002','004','006','008','010','012','014','016','018','020','022','024',...
'026','028','030','032','034','036','038','040','042','044','046','048',...
'050','052','054','056','058','060','062','064','066','068','070','072',...
'074','076','078','080','082','084','086','088','090','092','094','096',...
'098','100','102','104','106','108','110','112','114','116','118','120',...
'122','124','126','128','130','132','134','136','138','140','142','144',...
'146','148','150','152','154','156','158','160','162','164','166','168',...
'170','172','174','176','178','180','182','184','186','188','190','192',...
'194','196','198','200','202','204','206','208','210','212','214','216',...
'218','220','222','224','226','228','230','232','234','236','238','240',...
'242','244','246','248','250','252','254','256','258','260','262','264',...
'266','268','270','272','274','276','278','280','282','284','286','288',...
'290','292','294','296','298','300','302','304','306','308','310','312',...
'314','316','318','320','322','324','326','328','330','332','334','336',...
'338','340','342','344','346','348','350','352','354','356','358','360'];

%CHARACTER(LEN = 21)		:: leave_S_file_name, leave_E_file_name

%CHARACTER(LEN = 13)		:: NomeS = "S_P_Branches_",&
%                         NomeE = "E_B_Branches_"

% Output file with Emissivities and transmissivities values for both polarization
%OPEN(2 ,FILE = 'EV_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
%OPEN(22,FILE = 'EH_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
%OPEN(3 ,FILE = 'TV_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
%OPEN(33,FILE = 'TH_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')

deltet = pi./nn;

% Initialization from values defined in MOD_FOREST_PARAMETERS
dsdw_trunks = DSDWW;
lcorr = CORRELATION_LENGTH;
sigmaz = SOIL_ROUGHNESS;
amvcil = VOLMOIST_TRUNKS;

Cot = zeros(nij,1);
for i = 1:nij
  theta = (i - .5)*deltet;
  Cot(i) = cos(theta)*sin(theta);
end

species = 1;
name = 3;

if (species  ==  1)
	% *********************************************** HARDWOOD
	%	Possible name values:
	%
	%	Aspen/alder/cottonwood/willow 
	%	Soft maple/birch 
	%	Mixed hardwood 
	%	Hard maple/oak/hickory/beech 
	switch name
        case 3
		% Generic leaf dimension
            thick = .02;
            radius = 2.;            									
            leave_S_file_name = 'S_P_LEAVES_2_02_RG_06';
            leave_E_file_name = 'E_B_LEAVES_2_02_RG_06';
        otherwise
            disp('No leaf dimension has been selected %')
    end % SELECT
else
	% *********************************************** SOFTWOOD
	%	Possible name values:
	%
	% 1 - Cedar/larch 
	%	2 - Douglas - fir 
	%	3 - True fir/hemlock 
	%	4 - Pine 
	%	5 - Spruce
	switch name
        case 2
            % Douglas - Fir
            thick = 3.;
            radius = .1;              
            leave_S_file_name = 'S_P_NEEDLE_D_IL';
            leave_E_file_name = 'E_B_NEEDLE_D_IL';
        case 4
            % Pine needle dimension
            thick = 18.;
            radius = .1;
            leave_S_file_name = 'S_P_NEEDLE_P_IL';
            leave_E_file_name = 'E_B_NEEDLE_P_IL';
        otherwise
            disp('No needle dimension has been selected')
    end % SELECT
end %ENDIF

% Range distribution for trunks dimensions
dbh_min = 5.;
dbh_max = 70.;  
dbh_mean = 33.;

% Standard deviation for a Gaussian pdf
std = 18.2;

% % Leaves matrix reading
% OPEN(13,file = PathF//leave_S_file_name,FORM = 'unformatted',status = 'old')
% OPEN(14,file = PathF//leave_E_file_name,FORM = 'unformatted',status = 'old')
% READ(13) Snn0(:,:,:,:,:)
% READ(14) Kenn0(:,:)
% CLOSE(13)
% CLOSE(14)
% 
% % Secondary Branches matrix reading
% OPEN(13,file = PathB//"S_P_SECONDARI",FORM = 'unformatted',status = 'old')
% OPEN(14,file = PathB//"E_B_SECONDARI",FORM = 'unformatted',status = 'old')
% READ(13) Sb20(:,:,:,:,:)
% READ(14) Keb20(:,:)
% CLOSE(13)
% CLOSE(14)

LAI_start = 1;
for ilai = LAI_start:10 %Nlai: DO

	% Growth algorithms
	%CALL FOREST_GROWTH(REAL(ilai),name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius,an0nn,Ntr,htr,dbh,an0tr,Ncyl,an0br,an0br2,Dry_Litter_biomass)
    [an0nn,Ntr,htr,dbh,an0tr,Ncyl,an0br,an0br2,Dry_Litter_biomass] = FOREST_GROWTH(ilai,name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius);
    
	% Leaves contribution
	Snn(:,:,:,:,:) = Snn0(:,:,:,:,:) * an0nn / NSTRATI;
	Kenn(:,:) = Kenn0(:,:) * an0nn / NSTRATI;

% 	IF(.NOT.ALLOCATED(Sigmaatr)) THEN
% 	ALLOCATE(Sigmaatr(NN/2,2,Ntr),STAT = AllocErr)
% 	  IF(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria per Sigmaatr'
% 	ENDIF
    Sigmaatr = zeros(nn/2,2,Ntr);

	Sb(:,:,:,:,:) = Sb20(:,:,:,:,:) * an0br2 / NSTRATI;
	Keb(:,:) = Keb20(:,:) * an0br2 / NSTRATI;

	% Cycle wich take into account all the trunks with different diameter and height values
	for itr = 1:Ntr
        % Trunks absorption effect
        %CALL TRUNKS_ABS(dbh(itr)/2,htr(itr),amvcil,dsdw_trunks,Sigmaatr(:,:,itr))
        [Sigmaatr(:,:,itr),Kes] = MOD_VEGERARION_FUNCTIONS_TRUNKS_ABS(dbh(itr)/2,htr(itr),amvcil,dsdw_trunks);
        % Branches contribution 
        for i = 1:Ncyl(itr)		
              %OPEN(11,file = PathB//NomeS//Chacil(i),FORM = 'unFORMatted',status = 'old')
              %READ(11) Sb0(:,:,:,:,:)
              %CLOSE(11)
              Sb(:,:,:,:,:) = Sb(:,:,:,:,:) + Sb0(:,:,:,:,:)*an0br(itr,i)/NSTRATI;

              %OPEN(12,file = PathB//NomeE//Chacil(i),FORM = 'unFORMatted',status = 'old')
              %READ(12) Keb0(:,:)
              %CLOSE(12)
              Keb(:,:) = Keb(:,:) + Keb0(:,:)*an0br(itr,i)/NSTRATI;
        end % DO    
    end % DO

	% Soil moisture cycle
	for ivmoist = 1:7 %moisture: DO 

          % Volumetric soil moisture
          vmoist = .05 + (ivmoist - 1)*.05;
          litter_gmoist = 0;
        %original ferrazzoli litter moisture estimation

        % Litter gravimetric moisture, correlated with volumetric soil moisure
        % based on Les Landes ground measurements
        %IF(vmoist >= 0.1 .AND. vmoist <= .35) THEN
          %litter_gmoist = 3.1561. * vmoist - 0.17224               

        %ELSEIF (vmoist < 0.1 ) THEN
              %litter_gmoist = vmoist
            %ELSE
              %litter_gmoist = .85
        %END IF

        %end original ferrazzoli litter moisture estimation

        for ibiomass = 1:1   % [Kg/m^2] biomass: DO 

              % Soil scattering matrix evaluation  
              %Fresh_Litter_biomass = Dry_Litter_biomass / (1 - litter_gmoist) !original ferrazzoli
              Fresh_Litter_biomass=0.;

              %CALL IEM(sigmaz,lcorr,vmoist,Sgg,RR,litter_gmoist,Fresh_Litter_biomass)
              [Sgg,rr] = MOD_SOIL_FUNCTIONS_IEM(sigmaz,lcorr,vmoist,litter_gmoist,Fresh_Litter_biomass);
              tau = 0.;
              fm = 2.*pi;
              for istoki = 1:2
                    for j = 1:nij
                        ja = nij*(istoki - 1) + j;
                        akznn = Kenn(j,istoki);
                        akzb = Keb(j,istoki);
                       %
                       % VEGETATION TRANSMISSIVITIES
                       %
                       % Needles transmissivity
                         if (akznn >= 1.) 
                               EN(ja) = exp(-akznn);
                         else
                               EN(ja) = 1. - akznn;
                         end % IF
                         EEN(ja) = EN(ja)^NSTRATI;
                           %
                           % Branches transmissivity
                         if (akzb >= 1.) 
                               EB(ja) = exp(-akzb);
                         else
                               EB(ja) = 1. -akzb;
                         end %IF
                         EEB(ja) = EB(ja)^NSTRATI;
                           % 
                           % Crown transmissivity
                         if (akznn + akzb >= 1.) 
                            EC(ja) = exp(-akznn - akzb);
                         else
                            EC(ja) = 1. - akznn - akzb;
                         end % IF
                         EEC(ja) = EC(ja)^NSTRATI;
                           % 
                           % Trunks transmissivity
                         theta=(j - .5)*deltet;
                         akatr = sum(sigmaatr(j,istoki,:)*an0tr(:));
                         Abstr(ja) = exp(-akatr)*exp(-tau/cos(theta));
			  
                         for istoks = 1:2
                            for i = 1:nij
                                ia = nij*(istoks - 1) + i;
                                Sv(ia,ja) = (Sb(i,j,1,istoks,istoki) + Snn(i,j,1,istoks,istoki))*fm;
                                Tv(ia,ja) = (Sb(i,j,2,istoks,istoki) + Snn(i,j,2,istoks,istoki))*fm;
                                Sg0(ia,ja) = Sgg(i,j,1,istoks,istoki)*fm;
                                if (ia == ja) 
                                    Tv(ia,ja) = Tv(ia,ja) + EC(ja);
                                    Sg0(ia,ja) = Sg0(ia,ja) + RR(j,istoki);
                                end %IF
                            end %DO      
                         end %DO       
                    end %DO      
              end %DO   

		  % Matrix doubling amond crown sublayers 
		  %CALL MATRIX_DOUBLING(Sv,Tv,NSTRATI2)
          [sv, tv] = MATRIX_DOUBLING_MULTIPLE_VEGETATION(Sv,Tv,NSTRATI2);

		  % Attenuated soil scattering
		  for ja = 1:NN
            Sg0(:,ja) = Sg0(:,ja)*Abstr(ja)*Abstr(:);
          end %DO   

		  % Matrix doubling between ground and forest
		  %CALL MATRIX_DOUBLING(Sv,Tv,Tv,Sg0,Svg1)
          svg1 = MATRIX_DOUBLING_SOIL(Sv,Tv,Tv,Sg0);

		  % Total vegetation scattering matrix
		  Stot(:,:) = Sv(:,:) + Svg1(:,:);

		  % Emissivity estimation
		  for istoki = 1:2
            for j = 1:nij
               Emi(j,istoki) = 1.;
			   Emiv(j,istoki) = 1.;
			   Emig(j,istoki) = 1.;
			   ja = nij*(istoki - 1) + j;
			   for istoks = 1:2
                    for i = 1:nij
                        ia = nij*(istoks - 1) + i;
                        cc = Cot(i)/Cot(j);
					   % Total emissivity
                        Emi(j,istoki) = Emi(j,istoki) - cc*Stot(ia,ja);
					   % Crown emissivity
                        Emiv(j,istoki) = Emiv(j,istoki) - cc*(Sv(ia,ja) + Tv(ia,ja));
					   % Bare soil emissivity
					   Emig(j,istoki) = Emig(j,istoki) - cc*Sg0(ia,ja);    
                    end %DO    
               end %DO  
            end %DO   
          end %DO 

		  % Soil and trunks emissivities, attenuated by the crown
		  for istoks = 1:2
			  for i = 1:nij
				  ia = nij*(istoks - 1) + i;
				  Emig1(i,istoks) = 0.;
				  Emitr(i,istoks) = 0.;
				  for istoki = 1:2
					  for j = 1:nij
						  ja = nij*(istoki - 1) + j;
						  % Attenuated soil emissivity
						  Emig1(i,istoks) = Emig1(i,istoks) + emig(j,istoki)*tv(ia,ja);
						  % Attenuated trunks emissivity
						  Emitr(i,istoks) = Emitr(i,istoks) + (1 - abstr(ja))*tv(ia,ja);
                      end %DO   
                  end %DO     
              end %DO   
          end %DO     

		  %Original ferrazzoli
		  %WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emi(1:11,1)
		  %WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emiv(1:11,1)
		  %WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emitr(1:11,1)
		  %WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emig1(1:11,1)

		  %WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emi(1:11,2)
		  %WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emiv(1:11,2)
		  %WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emitr(1:11,2)
		  %WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emig1(1:11,2)
		  %Original ferrazzoli


		  % j = 6 => 27.5 deg incidence angle
		  disp(ilai,vmoist,Emi(6,1),Emiv(6,1),Emitr(6,1),Emig1(6,1))

		  %WRITE(2,101) ilai,vmoist,Emiv(6,1)
		  %WRITE(2,101) ilai,vmoist,Emitr(6,1)
		  %WRITE(2,101) ilai,vmoist,Emig1(6,1)

		  disp(ilai,vmoist,Emi(6,2),Emiv(6,2),Emitr(6,2),Emig1(6,2))
		  %WRITE(22,101) ilai,Emiv(6,2)
		  %WRITE(22,101) ilai,Emitr(6,2)
		  %WRITE(22,101) ilai,Emig1(6,2)

		  %WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEN(1:11)
		  %WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEB(1:11)
		  %WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Abstr(1:11)		
		  %WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEC(1:11)*Abstr(1:11)

		  %WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEN(19:29)
		  %WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEB(19:29)
		  %WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Abstr(19:29)
		  %WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEC(19:29)*Abstr(19:29)

  %101 FORMAT(I2, 2(1x,f4.2), 1x, f5.3, 1x,11(1x, f8.4))
  %102 FORMAT(I2, 5(1x,f5.3))

            end %DO biomass
        end %DO moisture 
    end %DO Nlai
%CLOSE(2)
%CLOSE(22)

%
