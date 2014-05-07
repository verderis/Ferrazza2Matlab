PROGRAM Emiforest

USE MOD_GLOBAL_PARAMETERS
USE MOD_MATRIX_DOUBLING
USE MOD_SOIL_FUNCTIONS
USE MOD_VEGETATION_FUNCTIONS

IMPLICIT NONE

! Variables to describe forest architecture
REAL, DIMENSION(:,:,:), ALLOCATABLE	:: Sigmaatr
REAL, DIMENSION(:,:), POINTER				:: an0br
REAL, DIMENSION(:), POINTER					:: htr,dbh,an0tr

INTEGER, DIMENSION(:), POINTER			:: Ncyl

REAL, DIMENSION(NN/2,NN/2,NM1,2,2)	:: Sgg

REAL, DIMENSION(NN/2,NN/2,2,2,2)		:: Sb,Sb0,Sb20,Snn,Snn0

REAL, DIMENSION(NN/2,2)							:: Emi,Emig,Emig1,Emitr,Emiv,keb,keb0,keb20,Kenn,Kenn0,RR
REAL, DIMENSION(NN,NN)							:: Sv,Tv,Sg0,Stot,Id,Svg1

! Variables which contain transmissivity values
REAL, DIMENSION(NN)									:: Abstr,EN,EEN,EB,EEB,EC,EEC,Cot(NN/2)

! VGeneric variables
REAL					  :: akznn,akzb,akatr,amvcil,vmoist,deltet,theta,an0br2,dsdw_trunks,litter_gmoist,Fresh_Litter_biomass
REAL					  :: an0nn,fm,lcorr,sigmaz,cc,tau,dbh_min,dbh_max,dbh_mean,std,thick,radius,Dry_Litter_biomass

INTEGER					:: i,j,ja,ia,ivmoist,itau,ilai,istoki,istoks
INTEGER					:: name,species,Ntr,itr,LAI_start,allocerr,ibiomass
            
! Path for scattering and extinction leaves matrix



CHARACTER(LEN = *), PARAMETER	:: PathF = "C:\Fran\Varios\PassiveModel\MATRICI_FORESTE\MATRICI_DI_SCATTERING\FOGLIE_C\"

!CHARACTER(LEN = *), PARAMETER	:: PathF = "C:\Documents and Settings\andrea\Documenti\Universita\Dottorato\Codice\&
!                                          Matrici Foreste\Matrici di scattering\Foglie L\"


! Path for scattering and extinction branches matrix

CHARACTER(LEN = *), PARAMETER	:: PathB = "C:\Fran\Varios\PassiveModel\MATRICI_FORESTE\MATRICI_DI_SCATTERING\&
LEN_25\FUNZIONE_PESO_C\"

!CHARACTER(LEN = *), PARAMETER	:: PathB = "C:\Documents and Settings\andrea\Documenti\Universita\Dottorato\Codice\&
!                                          Matrici Foreste\Matrici di scattering\Len 25\Funzione Peso L\"

CHARACTER(LEN = 3), DIMENSION(180):: Chacil = (/'002','004','006','008','010','012','014','016','018','020','022','024',&
																							 &'026','028','030','032','034','036','038','040','042','044','046','048',&
																							 &'050','052','054','056','058','060','062','064','066','068','070','072',&
																							 &'074','076','078','080','082','084','086','088','090','092','094','096',&
																							 &'098','100','102','104','106','108','110','112','114','116','118','120',&
																							 &'122','124','126','128','130','132','134','136','138','140','142','144',&
																							 &'146','148','150','152','154','156','158','160','162','164','166','168',&
																							 &'170','172','174','176','178','180','182','184','186','188','190','192',&
																							 &'194','196','198','200','202','204','206','208','210','212','214','216',&
																							 &'218','220','222','224','226','228','230','232','234','236','238','240',&
																							 &'242','244','246','248','250','252','254','256','258','260','262','264',&
																							 &'266','268','270','272','274','276','278','280','282','284','286','288',&
																							 &'290','292','294','296','298','300','302','304','306','308','310','312',&
																							 &'314','316','318','320','322','324','326','328','330','332','334','336',&
																					 &'338','340','342','344','346','348','350','352','354','356','358','360'/)
 
CHARACTER(LEN = 21)		:: leave_S_file_name, leave_E_file_name

CHARACTER(LEN = 13)		:: NomeS = "S_P_Branches_",&
                         NomeE = "E_B_Branches_"

! Output file with Emissivities and transmissivities values for both polarization
OPEN(2 ,FILE = 'EV_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
OPEN(22,FILE = 'EH_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
OPEN(3 ,FILE = 'TV_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')
OPEN(33,FILE = 'TH_15_HW_L_NT.txt',FORM = 'FORMatted',STATUS = 'unknown')

deltet = PI/NN

! Initialization from values defined in MOD_FOREST_PARAMETERS
dsdw_trunks = DSDWW
lcorr = CORRELATION_LENGTH
sigmaz = SOIL_ROUGHNESS
amvcil = VOLMOIST_TRUNKS

DO i = 1,NIJ
  theta = (i - .5)*deltet
  Cot(i) = COS(theta)*SIN(theta)
END DO

species = 1
name = 3

IF(species  ==  1) THEN
	! *********************************************** HARDWOOD
	!	Possible name values:
	!
	!	Aspen/alder/cottonwood/willow 
	!	Soft maple/birch 
	!	Mixed hardwood 
	!	Hard maple/oak/hickory/beech 
	SELECT CASE (name)
	CASE (3)
		! Generic leaf dimension
		thick = .02
		radius = 2.              									
		leave_S_file_name = "S_P_LEAVES_2_02_RG_06"
		leave_E_file_name = "E_B_LEAVES_2_02_RG_06"
  CASE DEFAULT
		STOP "No leaf dimension has been selected !"
	END SELECT
ELSE
	! *********************************************** SOFTWOOD
	!	Possible name values:
	!
	! 1 - Cedar/larch 
	!	2 - Douglas - fir 
	!	3 - True fir/hemlock 
	!	4 - Pine 
	!	5 - Spruce
	SELECT CASE (name)
	CASE (2)
		! Douglas - Fir
		thick = 3.
		radius = .1              
		leave_S_file_name = "S_P_NEEDLE_D_IL"
		leave_E_file_name = "E_B_NEEDLE_D_IL"
  CASE (4)
		! Pine needle dimension
		thick = 18.
		radius = .1
		leave_S_file_name = "S_P_NEEDLE_P_IL"
		leave_E_file_name = "E_B_NEEDLE_P_IL"
  CASE DEFAULT
		STOP "No needle dimension has been selected !"
	END SELECT
ENDIF

! Range distribution for trunks dimensions
dbh_min = 5.
dbh_max = 70.  
dbh_mean = 33.

! Standard deviation for a Gaussian pdf
std = 18.2

! Leaves matrix reading
OPEN(13,file = PathF//leave_S_file_name,FORM = 'unformatted',status = 'old')
OPEN(14,file = PathF//leave_E_file_name,FORM = 'unformatted',status = 'old')
READ(13) Snn0(:,:,:,:,:)
READ(14) Kenn0(:,:)
CLOSE(13)
CLOSE(14)

! Secondary Branches matrix reading
OPEN(13,file = PathB//"S_P_SECONDARI",FORM = 'unformatted',status = 'old')
OPEN(14,file = PathB//"E_B_SECONDARI",FORM = 'unformatted',status = 'old')
READ(13) Sb20(:,:,:,:,:)
READ(14) Keb20(:,:)
CLOSE(13)
CLOSE(14)


LAI_start = 1
Nlai: DO ilai = LAI_start,10

	! Growth algorithms
	CALL FOREST_GROWTH(REAL(ilai),name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius,an0nn,Ntr,htr,dbh,an0tr,Ncyl,an0br,an0br2,Dry_Litter_biomass)

	! Leaves contribution
	Snn(:,:,:,:,:) = Snn0(:,:,:,:,:) * an0nn / NSTRATI
	Kenn(:,:) = Kenn0(:,:) * an0nn / NSTRATI   

	IF(.NOT.ALLOCATED(Sigmaatr)) THEN
	ALLOCATE(Sigmaatr(NN/2,2,Ntr),STAT = AllocErr)
	  IF(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria per Sigmaatr'
	ENDIF

	Sb(:,:,:,:,:) = Sb20(:,:,:,:,:) * an0br2 / NSTRATI
	Keb(:,:) = Keb20(:,:) * an0br2 / NSTRATI

	! Cycle wich take into account all the trunks with different diameter and height values
	DO itr = 1,Ntr
	  ! Trunks absorption effect
    CALL TRUNKS_ABS(dbh(itr)/2,htr(itr),amvcil,dsdw_trunks,Sigmaatr(:,:,itr))
	  ! Branches contribution 
	  DO i = 1,Ncyl(itr)		

		  OPEN(11,file = PathB//NomeS//Chacil(i),FORM = 'unFORMatted',status = 'old')
		  READ(11) Sb0(:,:,:,:,:)
		  CLOSE(11)
		  Sb(:,:,:,:,:) = Sb(:,:,:,:,:) + Sb0(:,:,:,:,:)*an0br(itr,i)/NSTRATI

		  OPEN(12,file = PathB//NomeE//Chacil(i),FORM = 'unFORMatted',status = 'old')
		  READ(12) Keb0(:,:)
		  CLOSE(12)
		  Keb(:,:) = Keb(:,:) + Keb0(:,:)*an0br(itr,i)/NSTRATI
	  END DO    
	END DO

	! Soil moisture cycle
	moisture: DO ivmoist = 1,7

	  ! Volumetric soil moisture
	  vmoist = .05 + (ivmoist - 1)*.05
	  litter_gmoist = 0
    

	!original ferrazzoli litter moisture estimation
	
	! Litter gravimetric moisture, correlated with volumetric soil moisure
    ! based on Les Landes ground measurements
    !IF(vmoist >= 0.1 .AND. vmoist <= .35) THEN
      !litter_gmoist = 3.1561. * vmoist - 0.17224               
	               
    !ELSEIF (vmoist < 0.1 ) THEN
		  !litter_gmoist = vmoist
		!ELSE
		  !litter_gmoist = .85
    !END IF
	
	!end original ferrazzoli litter moisture estimation

    biomass: DO ibiomass = 1,1   ! [Kg/m^2]

      ! Soil scattering matrix evaluation  
      !Fresh_Litter_biomass = Dry_Litter_biomass / (1 - litter_gmoist) !original ferrazzoli
      Fresh_Litter_biomass=0.
  	  CALL IEM(sigmaz,lcorr,vmoist,Sgg,RR,litter_gmoist,Fresh_Litter_biomass)
      tau = 0.
		  fm = 2.*PI
		  DO istoki = 1,2
		   DO j = 1,nij
			   ja = nij*(istoki - 1) + j
			   akznn = Kenn(j,istoki)
			   akzb = Keb(j,istoki)
			   !
			   ! VEGETATION TRANSMISSIVITIES
			   !
			   ! Needles transmissivity
			   IF(akznn .GE. 1.) THEN
				   EN(ja) = EXP(-akznn)
			   ELSE
				   EN(ja) = 1. - akznn
			   END IF
			   EEN(ja) = EN(ja)**NSTRATI
			   !
			   ! Branches transmissivity
			   IF(akzb .GE. 1.) THEN
				   EB(ja) = EXP(-akzb)
			   ELSE
				   EB(ja) = 1. -akzb
			   END IF
			   EEB(ja) = EB(ja)**NSTRATI	
			   ! 
			   ! Crown transmissivity
			   IF(akznn + akzb .GE. 1.) THEN
				   EC(ja) = EXP(-akznn - akzb)
			   ELSE
				   EC(ja) = 1. - akznn - akzb
			   END IF
			   EEC(ja) = EC(ja)**NSTRATI		
			   ! 
			   ! Trunks transmissivity
         theta=(j - .5)*deltet
         akatr = SUM(sigmaatr(j,istoki,:)*an0tr(:))
         Abstr(ja) = EXP(-akatr)*EXP(-tau/COS(theta))
			  
			   DO istoks = 1,2
				   DO i = 1,nij
					   ia = nij*(istoks - 1) + i
					   Sv(ia,ja) = (Sb(i,j,1,istoks,istoki) + Snn(i,j,1,istoks,istoki))*fm
					   Tv(ia,ja) = (Sb(i,j,2,istoks,istoki) + Snn(i,j,2,istoks,istoki))*fm
 					   Sg0(ia,ja) = Sgg(i,j,1,istoks,istoki)*fm
					   IF(ia .EQ. ja) THEN
						   Tv(ia,ja) = Tv(ia,ja) + EC(ja)
						   Sg0(ia,ja) = Sg0(ia,ja) + RR(j,istoki)
					   END IF
				   END DO      
			   END DO       
		   END DO      
		  END DO   

		  ! Matrix doubling amond crown sublayers 
		  CALL MATRIX_DOUBLING(Sv,Tv,NSTRATI2)

		  ! Attenuated soil scattering
		  DO ja = 1,NN
		   Sg0(:,ja) = Sg0(:,ja)*Abstr(ja)*Abstr(:)  
		  END DO   

		  ! Matrix doubling between ground and forest
		  CALL MATRIX_DOUBLING(Sv,Tv,Tv,Sg0,Svg1)

		  ! Total vegetation scattering matrix
		  Stot(:,:) = Sv(:,:) + Svg1(:,:)

		  ! Emissivity estimation
		  DO istoki = 1,2
		   DO j = 1,nij
			   Emi(j,istoki) = 1.
			   Emiv(j,istoki) = 1.
			   Emig(j,istoki) = 1.
			   ja = nij*(istoki - 1) + j
			   DO istoks = 1,2
				   DO i = 1,nij
					   ia = nij*(istoks - 1) + i
					   cc = Cot(i)/Cot(j)
					   ! Total emissivity
					   Emi(j,istoki) = Emi(j,istoki) - cc*Stot(ia,ja)
					   ! Crown emissivity
					   Emiv(j,istoki) = Emiv(j,istoki) - cc*(Sv(ia,ja) + Tv(ia,ja))
					   ! Bare soil emissivity
					   Emig(j,istoki) = Emig(j,istoki) - cc*Sg0(ia,ja)    
				   END DO    
			   END DO  
		   END DO   
		  END DO 

		  ! Soil and trunks emissivities, attenuated by the crown
		  DO istoks = 1,2
			  DO i = 1,nij
				  ia = nij*(istoks - 1) + i
				  Emig1(i,istoks) = 0.
				  Emitr(i,istoks) = 0.
				  DO istoki = 1,2
					  DO j = 1,nij
						  ja = nij*(istoki - 1) + j
						  ! Attenuated soil emissivity
						  Emig1(i,istoks) = Emig1(i,istoks) + emig(j,istoki)*tv(ia,ja) 
						  ! Attenuated trunks emissivity
						  Emitr(i,istoks) = Emitr(i,istoks) + (1 - abstr(ja))*tv(ia,ja)
					  END DO   
				  END DO     
			  END DO   
		  END DO     

		  !Original ferrazzoli
		  !WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emi(1:11,1)
		  !WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emiv(1:11,1)
		  !WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emitr(1:11,1)
		  !WRITE(2,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emig1(1:11,1)

		  !WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emi(1:11,2)
		  !WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emiv(1:11,2)
		  !WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emitr(1:11,2)
		  !WRITE(22,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Emig1(1:11,2)
		  !Original ferrazzoli


		  ! j = 6 => 27.5 deg incidence angle
		  WRITE(2,102) ilai,vmoist,Emi(6,1),Emiv(6,1),Emitr(6,1),Emig1(6,1)

		  !WRITE(2,101) ilai,vmoist,Emiv(6,1)
		  !WRITE(2,101) ilai,vmoist,Emitr(6,1)
		  !WRITE(2,101) ilai,vmoist,Emig1(6,1)

		  WRITE(22,102) ilai,vmoist,Emi(6,2),Emiv(6,2),Emitr(6,2),Emig1(6,2)
		  !WRITE(22,101) ilai,Emiv(6,2)
		  !WRITE(22,101) ilai,Emitr(6,2)
		  !WRITE(22,101) ilai,Emig1(6,2)

		  !WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEN(1:11)
		  !WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEB(1:11)
		  !WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Abstr(1:11)		
		  !WRITE(3,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEC(1:11)*Abstr(1:11)

		  !WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEN(19:29)
		  !WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEB(19:29)
		  !WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,Abstr(19:29)
		  !WRITE(33,101) ilai,vmoist,litter_gmoist,Fresh_Litter_biomass,EEC(19:29)*Abstr(19:29)

  101 FORMAT(I2, 2(1x,f4.2), 1x, f5.3, 1x,11(1x, f8.4))
  102 FORMAT(I2, 5(1x,f5.3))

    END DO biomass
	END DO moisture 
END DO Nlai
CLOSE(2)
CLOSE(22)

CONTAINS
!
!******************************************************************************
!
! *********** INPUT
!
! lai	    : LAI values																									[cm^2/cm^2] 
! species : (1)  Hardood,  (2)  Softwood
! name    : Kind of selected forest species
! dbh_min : Minimum value for tree diameter used for the distribution   [cm]
! dbh_man : Maximum value for tree diameter used for the distribution   [cm]
! dbh_mean: Mean value for tree diameter used for the distribution      [cm]
! std     : Gaussian Standard deviation
! thick   : Needles length or leaves thickness								          [cm]
! radius  : Needles or leaves radius																		[cm]
!
! *********** OUTPUT
!
! an0nn       : Leaves surface density					                        [cm^2/cm^2] 
! Trees_group	: Trees number for at a specific dbh									    [#/cm^2] 
! an0br				: Primary branches densities                              [cm^2/cm^2]
! an0br2			: Secondary branches densities                            [cm^2/cm^2]
! htr					: Trees height at a specific dbh						              [cm]
! dbh					: Distribution diameters trees														[cm]
! an0tr				: Trees densities									                        [#/cm^2]
!
!
SUBROUTINE FOREST_GROWTH(lai,name,species,dbh_min,dbh_max,dbh_mean,std,thick,radius,an0nn,Trees_group,htr,dbh,an0tr,Ncyl,an0br,an0br2,Dry_Litter_biomass)
! Dummy argument declaration
!
! INPUT
REAL, INTENT(IN)								:: lai,dbh_min,dbh_max,dbh_mean,std,thick,radius
INTEGER, INTENT(IN)							:: name,species
!
! OUTPUT
REAL, DIMENSION(:,:), POINTER	  :: an0br
REAL, DIMENSION(:), POINTER		  :: htr,dbh,an0tr
INTEGER, DIMENSION(:), POINTER	:: Ncyl

REAL, INTENT(OUT)								:: an0nn,an0br2,Dry_Litter_biomass
INTEGER, INTENT(OUT)						:: Trees_group

! Local variables declaration
REAL, DIMENSION(:), ALLOCATABLE	:: Bio_tot,Bio_tr,Bio_f,Bio_b,Vol_tot,Vol_tr,Vol_f,Vol_b,dmax
REAL, DIMENSION(:), ALLOCATABLE	:: Bio_F_tot,Bio_F_tr,Bio_F_f,Bio_F_b

REAL, DIMENSION(5,2,2) :: biomassa
REAL, DIMENSION(4,2,2) :: component
REAL, DIMENSION(2,3)   :: A,n,bm,b0

REAL		:: area,dbh_step,rhoDB,rhoDT,amoiB,amoiT,inf,sup,rho_b,rho_t,step,db,vol
REAL		:: Trees_Number,rhoDF,amoiF,rho_f,secondary_branches,lbran,Foliage_biomass

INTEGER	:: AllocErr,up_boundary

! Constants values for dry and fresh matter
DATA dbh_step,rhoDB,rhoDT,rhoDF,amoiB,amoiT,amoiF/5., 0.4, 0.4, 0.3, 0.6, 0.5, 0.6/ 

DATA biomassa/ - 2.2094, - 1.9123,  - 2.48,   - 2.0127, 0.,     &   ! 4 Hardwood species
						 & 2.3867, 2.3651,  2.4835, 2.4342, 0.,     &
						 & - 2.0336, - 2.2304,  - 2.5384, - 2.5356, - 2.0773, &   ! 5 Softwood species
						 & 2.2592, 2.4435,  2.4814, 2.4349, 2.3323/

DATA component/ - 4.0813, - 1.6911,  - 2.0129,  - 0.3065, &         ! Hardwood components coefficients 
							& 5.8816, 0.816,   - 1.6805,  - 5.424,  &
							& - 2.9584, - 1.5619,  - 2.098,   - 0.3737, &         ! Softwood components coefficients 
							& 4.4766, 0.6614,  - 1.1432,  - 1.8055/

! Constants for Karam function, which determines the branches orientation
DATA A,n,bm,b0/2.42,2.6,2.12,2.3,1.98,2.01,&
			        &3.66,5.84,3.29,3.99,2.67,3.1,&
					 	  &.47,.39,.58,.58,.57,.6,&
			        &1.,1.,0.,0.,0.,0./  

! Complete leaf area
area = 2*PI*radius*radius*thick*(1/radius + 1/thick)  

! projected leaf area
IF(species  ==  1) THEN
	! Hardwood
  area = area/2
ELSE
  ! Softwood
  area = area/2.5
ENDIF

! Leaves density
an0nn = LAI/area 

! Number of differents consideres diameters values
Trees_group = NINT((dbh_max - dbh_min)/dbh_step)    

! Memory allocation for used pointers
ALLOCATE(an0tr(Trees_group),Bio_tot(Trees_group),Bio_tr(Trees_group),Bio_f(Trees_group),Bio_b(Trees_group),Bio_F_tot(Trees_group),&
				&Bio_F_tr(Trees_group),Bio_F_f(Trees_group),Bio_F_b(Trees_group),Vol_tot(Trees_group),Vol_tr(Trees_group),Vol_b(Trees_group),&
				&Vol_f(Trees_group),htr(Trees_group),dbh(Trees_group),dmax(Trees_group),Ncyl(Trees_group),STAT = AllocErr)

IF(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria per an0tr'

! Diameter trunks step
dbh_step = (dbh_max - dbh_min)/Trees_group

! Integration probability function at each interval
DO i = 1,Trees_group  
	inf = dbh_min + (i - 1) * dbh_step
	sup = inf + dbh_step
	an0tr(i) = INTGAUSS(inf,sup,dbh_mean,std)
	dbh(i) = inf + dbh_step / 2;
	!
	! All the weight are expressed in grams
	!
	! Total dry biomass
	Bio_tot(i) = 1000 * EXP(biomassa(name,1,species) + biomassa(name,2,species) * LOG(dbh(i))) 
  ! Trunks dry biomass 
	Bio_tr(i) = Bio_tot(i) * (EXP(component(3,1,species) + component(3,2,species) / dbh(i)) + &
 						  EXP(component(4,1,species) + component(4,2,species) / dbh(i)))
  ! Foliage dry biomass
	Bio_f(i) = Bio_tot(i) * EXP(component(1,1,species) + component(1,2,species) / dbh(i))
ENDDO  

! Branches dry biomass
Bio_b(:) = Bio_tot(:) - Bio_tr(:) - Bio_f(:)

! Number of tree at each interval
an0tr(:) = an0tr(:) + (1 - SUM(an0tr(:)))/Trees_group

! Fresh biomass
Bio_F_tr(:) = Bio_tr(:) / (1 - amoiT)
Bio_F_b(:) = Bio_b(:) / (1 - amoiB)
Bio_F_f(:) = Bio_f(:) / (1 - amoiF)
Bio_F_tot(:) = Bio_F_tr(:) + Bio_F_b(:) + Bio_F_f(:)

! Efficacious density of fresh matter
rho_b = rhoDB / ((1 - amoiB) + rhoDB*amoiB) 
rho_t = rhoDT / ((1 - amoiT) + rhoDT*amoiT)  	
rho_f = rhoDF / ((1 - amoiF) + rhoDF*amoiF)  	

! Volumes of each dielectric component
Vol_tr(:) = Bio_F_tr(:) / rho_t    
Vol_b(:) = Bio_F_b(:) / rho_b  
Vol_f(:) = Bio_F_f(:) / rho_f
Vol_tot(:) = Vol_tr(:) + Vol_b(:) + Vol_f(:)

! Trunks height
htr(:) = 4 * Vol_tr(:) / (PI * dbh(:) * dbh(:))

! Foliage biomass link with LAI values
IF(species  ==  1) THEN
	Foliage_biomass = LAI / 1.4889                        ! [t/Ha]
  Dry_Litter_biomass = .01 *(1.2463 * Foliage_biomass)  ! [g/cm^2]
ELSE
  Foliage_biomass = LAI / 0.4314                        ! [t/Ha]
  Dry_Litter_biomass = .01 *(0.2498 * Foliage_biomass)  ! [g/cm^2]
ENDIF

! Total number trees
Trees_Number = 1.e-2 * Foliage_biomass / (SUM(Bio_f(:) * an0tr(:))) ! [#/cm^2]

! Trees number for each selected dbh value
an0tr(:) = an0tr(:) * Trees_Number

! Branches maximum value
IF(species  ==  1) THEN
  dmax(:) = (dbh(:) / 4)
ELSE
	dmax(:) = 10**( -0.468 + 0.803 * LOG10(1.0 + 1.2*dbh(:)))
ENDIF

! Print of parameter values
OPEN(192,file = "Biomasse.txt",FORM = 'FORMatted',status = 'unknown')

IF(.true.) THEN

  DO i = 1,Trees_group

    WRITE(192,11)  lai,dbh(i),htr(i)/1e2,NINT(an0tr(i)*1e8),&
					         Bio_tot(i)* an0tr(i) *1e2,Bio_tr(i)* an0tr(i) *1e2,Bio_f(i)* an0tr(i) *1e2,Bio_b(i)* an0tr(i) *1e2,&
					         Bio_F_tot(i)* an0tr(i) *1e2,Bio_F_tr(i)* an0tr(i) *1e2,Bio_F_f(i)* an0tr(i) *1e2,Bio_F_b(i)* an0tr(i) *1e2,&
					         Vol_tot(i)* an0tr(i) *1e2,Vol_tr(i)* an0tr(i) *1e2,Vol_f(i)* an0tr(i) *1e2,Vol_b(i)* an0tr(i) *1e2
  ENDDO

  WRITE(192,*)
  WRITE(192,*)

11 FORMAT(F4.1,1X,2(F6.2,1X),I3,12(F9.3))    
  
  IF(lai  ==  10) CLOSE(192)

ENDIF

step = 0.2

! "secondary_branches * 100" is the percentage component of secondary branches presence
secondary_branches = .3

! Secondary branches density
an0br2 = secondary_branches * SUM(Vol_b(:) * an0tr(:)) / (PI * 0.3 * 0.3 * 20)

Ncyl(:) = NINT(dmax(:) / step)
up_boundary = MAXVAL(Ncyl(:))

ALLOCATE(an0br(Trees_group,up_boundary),STAT = AllocErr)
IF(AllocErr .NE. 0.) STOP 'Errore nell''allocazione di memoria'
an0br(:,:) = 0.  
!
! Primary branches densities for each dbh value
!
DO i = 1,Trees_group
	inf = 0.
	sup = step
	db = (sup + inf)/2.					! dimameter value
  IF(lb/db > 50) THEN
    lbran = db*50.
  ELSE
    lbran = lb
  ENDIF
	vol = pi*db*db*lbran/4.   ! volume cylinder
	DO j = 1,Ncyl(i)	    	
		an0br(i,j) = INTKARAM(inf,sup,A(1,3)/dmax(i),n(1,3),b0(1,3)*dmax(i),bm(1,3)*dmax(i)) *&
	               (1 - secondary_branches) * Vol_b(i) * an0tr(i) / vol												! [#/cm^2]      
	  inf = sup
	  sup = sup + step
	  db = db + step
    IF(lb/db > 50) THEN
      lbran = db*50.
    ELSE
      lbran = lb
    ENDIF
    vol = pi*db*db*lbran/4.
	ENDDO
ENDDO

DEALLOCATE(Bio_tot,Bio_tr,Bio_f,Bio_b,Vol_tot,Vol_tr,Vol_b,Vol_f,dmax,Bio_F_tot,Bio_F_tr,Bio_F_f,Bio_F_b)

RETURN
END SUBROUTINE FOREST_GROWTH
!
! ************ Integral of Gauss function ************
!
REAL FUNCTION INTGAUSS(inf,sup,mean,std)
! Dummy variables declaration
REAL	:: inf,sup,mean,std
! Local variables declaration
REAL	:: b,passo,ck
INTEGER :: i,nint

NINT = 10

INTGAUSS = 0.
passo = (sup - inf)/nint
nint = nint + 1
b = inf

DO i = 1,nint

  IF(i  ==  1 .OR. i  ==  nint) THEN
		ck = 1.
	ELSE
		ck = 2.
	ENDIF

	INTGAUSS = INTGAUSS + ck * EXP(( - (b - mean)**2)/(2*std*std))
	b = b + passo

END DO  

INTGAUSS = (1/(std*SQRT(2*PI))) * INTGAUSS * passo/2.

RETURN
END FUNCTION INTGAUSS
!
! ************ Integral of Karam function ************
!
REAL FUNCTION INTKARAM(inf,sup,A,n,b0,bm)
! Dummy variables declaration
REAL	:: inf,sup,A,n,b0,bm
! Local variables declaration
REAL	:: b,passo,ck,arg
INTEGER :: i,nint

  nint = 10
  INTKARAM = 0.
  passo = (sup - inf)/nint
  nint = nint + 1 
  b = inf

  DO i = 1,nint

		IF(i  ==  1 .OR. i  ==  nint) THEN
			ck = 1.
		ELSE
			ck = 2.
		ENDIF

		arg = ABS((b - bm)/(b0 - bm))
		! Essendo n un numero irrazionale il radicando (arg) deve essere positivo, inoltre per motivi di
		! precisione di calcolo quando l'argomento tende ad 1. il prodotto arg*PI/2 può eccedere PI/2. portando
		! il coseno ad un valore negativo e non regolare per un indice irrazionale.
		IF (arg .LT. 1.) INTKARAM = INTKARAM + ck*A*COS((PI/2.)*arg)**n	   
		b = b + passo

  END DO  

  INTKARAM = INTKARAM*passo/2

RETURN
END FUNCTION INTKARAM

END PROGRAM Emiforest