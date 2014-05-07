! *************************************************************************************
! MOD_GLOBAL_PARAMETER contains the constants definition necessary for the correct 
! execution of the program
! *************************************************************************************
MODULE MOD_GLOBAL_PARAMETERS

PUBLIC

!***************************************************************************************
! GENERIC CONSTANTS INPUT PARAMETER
!***************************************************************************************

REAL, PARAMETER	::  F = 1.4												! Frequency

! ****************************** VEGETATION ********************************
REAL, PARAMETER	::   LB = 25.				,& ! Branches segment length
					 DSDWW = .4				,& ! Wood dry matter density
                     DSDWF = .3             ,& ! Vegetation dry matter density
                     VOLMOIST_TRUNKS = .5   ,& ! Wood water content (Trunks)
                     VOLMOIST_MATRIX = .6      ! Branches and leaves

! Selection of permittivity alghoritm
INTEGER, PARAMETER  :: PERMITTIVITY = 2		! 1-> Matzler 
											! 2-> Ulaby 

! ********************************* SOIL *********************************** 
REAL	::  CORRELATION_LENGTH = 5.		,&  ! Soil correlation length
			SOIL_ROUGHNESS = 1.5		,&	! Soil roughness
            TSOL = 290.					,&  ! Ground temperature
			SAND = 0.8					,&  ! SAND percentage
			CLAY = 0.1					,&  ! CLAY percentage
			ROB = 1.25   

! Selection of Correlation values
INTEGER, PARAMETER	:: CORRELATION = 1	! 1-> Exponential 
										! 2-> Gaussian

! Generic costants
REAL(8), PARAMETER  :: DPI = 3.141592653589793238D0  ! Pi Greek double precision
REAL, PARAMETER     :: PI = 3.1415927								 ! Pi Greek single precision
INTEGER, PARAMETER	:: NMAX0 = 100						       ! Maximum order for Bessel series 

!***************************************************************************************
! SAMPLE FREQUENCY OF BISTATIC SCATTERING FUNCTION
!***************************************************************************************
! NN = 2*(number of discrete intervals of incidence and scattering off - nadir angles);
INTEGER, PARAMETER		:: NN = 36 
! NM = number of Fourier components for dependence on (phi_s - phi) of scattering;
INTEGER, PARAMETER		:: NM = 64	
! NIJ number of intervals in theta and theta_s between 0 and PI/2   
INTEGER, PARAMETER		:: NIJ = NN/2     
! NM1 number of intervals in phi_s - phi between 0 and PI  
INTEGER, PARAMETER		:: NM1 = NM/2 + 1				
	
! NN = Number of Fourier elements
INTEGER, PARAMETER		:: NF = -NINT(LOG(REAL(NM))/LOG(2.))

! Number of sublayer for matix doubling
INTEGER, PARAMETER    :: NSTRATI2 = 11
INTEGER, PARAMETER    :: NSTRATI = 2**nstrati2

!***************************************************************************************
! SAMPLE FREQUENCY FOR LEGENDRE-GAUSS INTEGRATION METHOD 
!***************************************************************************************
! Frequency sampling in integration method, two points o four points sampling points are foreseen.
INTEGER, PARAMETER		:: NS = 2
! Coefficients for two points sampling 
REAL, DIMENSION(NS), PARAMETER	:: OM  = (/1, 1/) ,& ! zeros of Legendre polynomials
								   CHI = (/.5773502692 , -.5773502692/)   ! weights of Gaussian quadrature
! Coefficients for four points sampling 
!REAL, DIMENSION(NS), PARAMETER :: OM  = (/.3478548451, .6521451549, .6521451549, .3478548451/),& ! zeros of Legendre polynomials
!																   CHI = (/.8611363116, .3399810436,-.3399810436,-.8611363116/)   ! weights of Gaussian quadrature


!***************************************************************************************
! EULERIAN ANGLES FOR ALL THE DIELECTRY BODIES
!***************************************************************************************
!
! The values of Select is required in Infinite_Length and Hollow routine
!
! Leaves, used with Raylegh-Gans and Pysical Optics approximation
REAL :: ALPHA1DIS = 15., NALPHADIS = 12, DALPHADIS = 30. ,&
			  BETA1DIS  = 5. , NBETADIS  = 9 , DBETADIS  = 10. ,&
			  GAMMA1DIS = 0. , NGAMMADIS = 1 , DGAMMADIS = 0.  

! Petioles, needles and secondary branches used with Rayleigh-Gans, Infinte length
! and Pysical Optics approximation
!
! Select = 'P'  Petioles
! Select = 'N'  Needles
! Select = 'R'  Secondary Branches
!
REAL :: ALPHA1PET = 15., NALPHAPET = 12, DALPHAPET = 30. ,&
			  BETA1PET  = 5. , NBETAPET  = 9 , DBETAPET  = 10. ,&
			  GAMMA1PET = 0. , NGAMMAPET = 1 , DGAMMAPET = 0. 

! Ears used with Infinite length approximation
!
! Select = 'E'  Ears
!
REAL :: ALPHA1EAR = 15., NALPHAEAR = 1,  DALPHAEAR = 0. ,&
		 	  BETA1EAR  = 0. , NBETAEAR  = 1 , DBETAEAR  = 0. ,&
			  GAMMA1EAR = 0. , NGAMMAEAR = 1 , DGAMMAEAR = 0. 

! Stems used with Infinite length
!
! Select = 'S'  Stem and Trunks
!
REAL :: ALPHA1STE = 15., NALPHASTE = 12, DALPHASTE = 30. ,&
			  BETA1STE  = 2. , NBETASTE  = 2 , DBETASTE  = 3.  ,&
			  GAMMA1STE = 0. , NGAMMASTE = 1 , DGAMMASTE = 0. 

! Used for branches Infinte length
!
! Select = 'B'  Primary branches
!
REAL :: ALPHA1BRA = 15., NALPHABRA = 12, DALPHABRA = 30. ,&
			  BETA1BRA = 5.  , NBETABRA  = 9 , DBETABRA  = 10.  ,&
			  GAMMA1BRA = 0. , NGAMMABRA = 1 , DGAMMABRA = 0. 


END MODULE MOD_GLOBAL_PARAMETERS