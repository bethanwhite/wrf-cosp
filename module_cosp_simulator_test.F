MODULE module_cosp_simulator_test 

 ! TODO: rename when working: module_cosp_wrf
 !-------------------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Model output for comparisons with CloudSat, CALIPSO, MODIS, MISR and ISCCP satellite data
 !
 ! The COSP simulator is (c) copyright by the British Crown / Met Office 2008
 ! See http://cfmip.metoffice.com/cosp/cosp.v0.3/  and
 ! Refer to Met_Office_licencse.text for m^2 / s^2
 !
 ! 
 !------------------------------------------------------------------------------------------------- 
 ! WRF IMPLEMENTATION VERSION 1.0, BETHAN WHITE, DECEMBER 2013 
 ! This is the first implementation of COSP into WRF, written for convection-permitting simulations. 
 ! Thanks go to Mirek Andrejczuk for his help with much of the technical WRF setup, and thanks to 
 ! J. Kay of UCAR for publishing the CAM-COSP interface, which was of great help.
 !
 ! This is not yet a general implementation and requires the following WRF settings:
 ! - Morrison two-moment bulk microphysics scheme
 ! - RRTM radiation scheme
 ! Unfortunately, due to the way that WRF handles nesting, some changes have had to be made to the native
 ! microphysics and radiation code in order to make certain variables available to COSP that are not available
 ! from the WRF data grid structure. In order to run with other mphys/radiation options, similar changes will 
 ! have to be made to those schemes.
 !
 ! The following notes describe what changes have been made:
 ! +++++++++ NOTES ON NATIVE CHANGES MADE TO WRF CODE GO HERE +++++++++++++++++++++++++++++++++++++++
 ! module_mp_morr_two_moment.F:
 ! - added hydrometero effective radii, ERADC, ERADI, ERADS, ERADR, ERADG 
 ! +++++++++ NOTES ON WRF STRUCTURE (WHERE COSP CODE IS) GO HERE ++++++++++++++++++++++++++++++++++++
 !
 !
 ! Obvious points for future updates include:
 ! - moving the logical flags and input parameters for the simluator code to a namelist
 ! - moving some input parameters to the namelist
 ! - an option for running with a convection scheme?
 ! - run code at n*rad_dt if running at every radiation step is too expensive?
 !
 ! TODO BEFORE IMPLEMENTING:
 ! - Check all REALS are set to dp (currently not assigned)
 ! - COSP simulator is currently called from each k-point in the tile - I don't know how this works
 ! - Need to check whether all outputs are point-data or not (ISCCP cloud types, etc?) - how does WRF Registry handle these?
 !-------------------------------------------------------------------------------------------------
 
 ! USE THE MET OFFICE COSP SIMULATOR
 USE mod_cosp_simulator

 ! USE WRF ERROR MODULE?
 !USE module_wrf_error


 ! USE WRF PHYSICS CONSTANTS:
  use module_model_constants, ONLY: &
                                    G, &        ! gravitational acceleration
                                    R => r_d    ! dry gas constant

 ! USE COSP CONSTANTS:
 USE mod_cosp_constants    !, ONLY: R_UNDEF, PARASOL_NREFL, i_cvcice, i_cvcliq, i_lscliq, i_lscice ! NOTE - actually want to read all constants!

 ! USE COSP TYPES TO CREATE DERIVED DATA TYPES
 USE mod_cosp_types

 IMPLICIT NONE

 CONTAINS
  subroutine wrf_cosp_simulator(specified                                & 
                      ,eradc, eradi, erads, eradr, eradg                 &
                      ,ozonemr, cldfra, u10, v10,  tsk, ht               &
                      ,qv, qc, qr, qi, qs, qg, qh                        &
                      ,phb, ph, pb, p, psfc, t, xland                    &
                      ,taucldc, taucldi, emiss, coszen                   &                         
                      ,rainncv, hailncv, snowncv, graupelncv             &
                      ,xlat, xlong                                       &
                      ,ids,ide, jds,jde, kds,kde                         &
                      ,ims,ime, jms,jme, kms,kme                         &
                      ,ips,ipe, jps,jpe, kps,kpe                         &
                      ,i_start,i_end,j_start,j_end,kts,kte               &
                      ,num_tiles, spec_zone, channel_switch              &
                      ,dt, curr_secs)              

   !------------ 1. VARIABLE DECLARATIONS: ------------------------------------------------------------------------

   ! INPUT VARIABLES -----------------------------------------------------------------------------------------------
   ! TODO: check intents - these should just be (IN)?
   LOGICAL,    INTENT(IN   )    :: specified
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradc
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradi
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::erads
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradr
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradg
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::cldfra
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::ozonemr
   REAL,INTENT(INOUT), DIMENSION(ims:ime,jms:jme)::u10
   REAL,INTENT(INOUT), DIMENSION(ims:ime,jms:jme)::v10
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::tsk
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::ht
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qv
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qc
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qr
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qi
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qs
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qg
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::qh
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::phb
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::ph
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::pb
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::p
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::psfc
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::t
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::xland
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::taucldc
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::taucldi
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::emiss
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme)::coszen
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::rainncv
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::hailncv
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::snowncv
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::graupelncv

   REAL, INTENT(IN)     :: xlat ( ims:ime, jms:jme )        ! latitudes on mass-point grid
   REAL, INTENT(IN)     :: xlong( ims:ime, jms:jme )        ! longitudes on mass-point grid

   INTEGER, INTENT(IN)  :: ids,ide, jds,jde, kds,kde                         &
                          ,ims,ime, jms,jme, kms,kme                         &
                          ,ips,ipe, jps,jpe, kps,kpe                         
   INTEGER, INTENT(IN)  :: kts,kte
   INTEGER, INTENT(IN)  :: num_tiles, spec_zone

   INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start,i_end,j_start,j_end    
   LOGICAL, OPTIONAL, INTENT(IN)    :: channel_switch

   REAL, INTENT(IN) :: dt, curr_secs  
  
   
   ! COSP input variables from WRF: most of these are passed from the WRF data grid structure, but 
   ! PLEASE NOTE that for variables not available from the grid, the code currently ONLY gathers these 
   ! from the MORRISON microphysics scheme and the RRTM radiation scheme. Using other physics schemes
   ! will cause the code to FAIL!!! USE AT YOUR PERIL!
   !
   !
   ! eradc, eradi, erads, eradr, eradg - cloud droplet, ice, snow, rain, graupel effective radii, micron (3D)
   ! cldfra - cloud fraction (3D)
   ! u10 - 10m u-wind component, m/s (2D)
   ! v10 - 10m v-wind component, m/s (2D)
   ! tsk - skin temperature, K (2D)
   ! ht - surface height, m (2D)
   ! qv, qc, qr, qi, qs, qg, qh - mass mixing ratios, kg/kg, for water vapour, cloud droplet, rain, ice,  snow, hail, graupel (3D)
   ! phb - base state geopotential, m^2 / s^2 (3D)
   ! ph - gepotential perturbation, m^2 / s^2 (3D)
   ! pb - base state pressure, Pa (3D)
   ! p - pressure perturbation, Pa (3D)
   ! psfc - surface pressure, Pa (2D)
   ! t - temperature, K (3D)
   ! xland - land use mask (2D)
   ! taucldc, taucldi - cloud droplet and ice optical depths (3D)
   ! emiss - surface emissivity (2D)
   ! coszen - the cosine of the solar zenith angle (2D)
   ! rainncv, hailncv, snowncv, graupelncv - sedimentation fluxes, kg/m^2/s, for rain, hail, snow, graupel (3D)
   ! xlat, xlong - latitude and longitude on mass points
   ! dt - time step (s)
   ! curr_secs - current forecast time (in seconds)

   ! OUTPUT VARIABLES -------------------------------------------------------------------------------------------------
   ! Output variable declarations go in here
   !
   !--------------------------------------------------------------------------------------------------------------------


   !--------- 2. WE NEED TO CALCULATE WRF TILE DIMENSIONS IN ORDER TO SET THE DIMENSIONS OF COSP INPUTS -----------------
   !--------- THE FOLLOWING CODE FOLLOWS THE DEFINITION OF ITS, ITE, JTS, JTE FROM THE MICROPHYSICS DRIVER --------------
   INTEGER :: its,ite,jts,jte,ij,sz
   LOGICAL :: channel
   ! Define the size of the variables passed to COSP
   ! Number of gridpoints
   INTEGER :: Npoints
   ! Number of levels
   INTEGER :: Nlevels
   ! Number of columns
   INTEGER :: Ncolumns

   channel = .FALSE.
   IF ( PRESENT ( channel_switch ) ) channel = channel_switch

   IF( specified ) THEN
      sz = spec_zone
   ELSE
      sz = 0
   ENDIF

   #ifndef RUN_ON_GPU
   DO ij = 1 , num_tiles
      IF (channel) THEN
         its = max(i_start(ij),ids)
         ite = min(i_end(ij),ide-1)
      ELSE
         its = max(i_start(ij),ids+sz)
         ite = min(i_end(ij),ide-1-sz)
      ENDIF
         jts = max(j_start(ij),jds+sz)
         jte = min(j_end(ij),jde-1-sz)
   #else
   DO ij = 1 , 1
      IF (channel) THEN
         its = max(ips,ids)
         ite = min(ipe,ide-1)
      ELSE
         its = max(ips,ids+sz)
         ite = min(ipe,ide-1-sz)
      ENDIF
         jts = max(jps,jds+sz)
         jte = min(jpe,jde-1-sz)
   #endif

   ! CALCULATE THE NUMBER OF GRID POINTS IN THE TILE
   ! Npoints = (jte-jts)*(ite-its) ! OK FOR CALLING COSP
   Npoints = 1 ! COSP is called once for each column !TODO: is this the case, or can we call it for the tile?
   ! BEWARE: CODE BELOW DEPENDS ON NPOINTS = 1 BECAUSE WE HAVE CHOSEN TO MOVE DATA ONTO THE COLUMN AND RUN COSP ONCE FOR EACH COLUMN
   ! IF YOU CHANGE NPOINTS, BE SURE TO FOLLOW THROUGH EVERYTHING ELSE CAREFULLY!

   ! Number of levels
   Nlevels = kte - kts

   ! Number of columns
   Ncolumns = Npoints

   !---------- DONE CALCULATING DIMENSIONS -----------------------------------------------------------------------------


   !---------- 3. DECALRE LOCAL VARIABLES:  -----------------------------------------------------------------------------

   !--------- SIMULATOR FLAGS - SET WHICH SIMULATORS TO RUN HERE --------------------------------------------------------
   LOGICAL :: lradar_sim = .false.   ! Radar
   LOGICAL :: llidar_sim = .false.   ! Lidar
   LOGICAL :: lmodis_sim = .false.   ! MODIS
   LOGICAL :: lmisr_sim = .false.    ! MISR
   LOGICAL :: lisccp_sim = .false.   ! ISCCP
   LOGICAL :: lrttov_sim = .false.   ! RRTOV

   !--------- DO YOU WANT SUB-COLUMN OUTPUT ? ---------------------------------------------------------------------------
   LOGICAL :: lfrac_out = .false.

   !--------- OUTPUT FLAGS REQUIRED BY COSP CONFIGURATION ---------------------------------------------------------------	
   !--------- THESE ARE ALL INITIALISED TO 'FALSE', THEN WE OVERWRITE THE VALUES ACCORDING TO WHICH SIMULATORS ARE SET --
   LOGICAL :: lalbisccp = .false.
   LOGICAL :: latb532 = .false.
   LOGICAL :: lboxptopisccp = .false.
   LOGICAL :: lboxtauisccp = .false.
   LOGICAL :: lcfad_dbze94 = .false.
   LOGICAL :: lcfad_lidarsr532 = .false.
   LOGICAL :: lclcalipso = .false.
   LOGICAL :: lclhcalipso = .false.
   LOGICAL :: lclisccp2 = .false.
   LOGICAL :: lcllcalipso = .false.
   LOGICAL :: lclmcalipso = .false.
   LOGICAL :: lcltcalipso = .false.
   LOGICAL :: lctpisccp = .false.
   LOGICAL :: ldbze94 = .false.
   LOGICAL :: lcltradar = .false.
   LOGICAL :: lcltradar2 = .false.
   LOGICAL :: ltauisccp = .false.
   LOGICAL :: ltclisccp = .false.
   LOGICAL :: lparasol_refl = .false.
   LOGICAL :: lclmisr = .false.
   LOGICAL :: lmeantbisccp = .false.
   LOGICAL :: lmeantbclrisccp = .false.
   LOGICAL :: lclcalipso2 = .false.
   LOGICAL :: lcltlidarradar = .false.
   LOGICAL :: lbeta_mol532 = .false.
   LOGICAL :: Llongitude = .false.
   LOGICAL :: Llatitude =.false.
   LOGICAL :: lcltmodis = .false.
   LOGICAL :: lclwmodis = .false.
   LOGICAL :: lclimodis = .false.
   LOGICAL :: lclhmodis = .false.
   LOGICAL :: lclmmodis = .false.
   LOGICAL :: lcllmodis = .false.
   LOGICAL :: ltautmodis = .false.
   LOGICAL :: ltauwmodis = .false.
   LOGICAL :: ltauimodis = .false.
   LOGICAL :: ltautlogmodis = .false.
   LOGICAL :: ltauwlogmodis = .false.
   LOGICAL :: ltauilogmodis = .false.
   LOGICAL :: lreffclwmodis = .false.
   LOGICAL :: lreffclimodis = .false.
   LOGICAL :: lpctmodis = .false.
   LOGICAL :: llwpmodis = .false.
   LOGICAL :: liwpmodis = .false.
   LOGICAL :: lclmodis = .false.
   LOGICAL :: ltbrttov = .false.  !! RTTOV OUTPUT 

   ! OVERWRITE THE OUTPUT FLAGS ACCORDING TO WHICH SIMULATORS ARE TO BE RUN:
   ! RADAR SIMULATOR
   if (lradar_sim) then
      lcfad_dbze94 = .true.
      ldbze94 = .true.
      lcltradar = .true.
      lcltradar2 = .true.
   end if

   ! RADAR AND LIDAR SIMULATORS TOGETHER
   if ((lradar_sim) .and. (llidar_sim)) then
      lclcalipso2 = .true.
      lcltlidarradar = .true.
   end if

   ! LIDAR SIMULATOR
   if (llidar_sim) then
      lcllcalipso = .true.
      lclmcalipso = .true.
      lcltcalipso = .true.
      lclcalipso = .true.
      lclhcalipso = .true.
      lcfad_lidarsr532 = .true.
      latb532 = .true.
      lparasol_refl = .true.
      lbeta_mol532 = .true.
   end if

   ! ISCCP SIMULATOR
   if (lisccp_sim) then
      lalbisccp = .true.
      lboxptopisccp = .true.
      lboxtauisccp = .true.
      lclisccp2 = .true.
      lctpisccp = .true.
      ltauisccp = .true.
      ltclisccp = .true.
      lmeantbisccp = .true.
      lmeantbclrisccp = .true.
   end if

   ! MISR SIMULATOR
   if (lmisr_sim) then
      lclmisr = .true.
   end if

   ! MODIS SIMULATOR
   if (lmodis_sim) then
      lcltmodis = .true.
      lclwmodis = .true.
      lclimodis = .true.
      lclhmodis = .true.
      lclmmodis = .true.
      lcllmodis = .true.
      ltautmodis = .true.
      ltauwmodis = .true.
      ltauimodis = .true.
      ltautlogmodis = .true.
      ltauwlogmodis = .true.
      ltauilogmodis = .true.
      lreffclwmodis = .true.
      lreffclimodis = .true.
      lpctmodis = .true.
      llwpmodis = .true.
      liwpmodis = .true.
      lclmodis = .true.
   end if
   ! DONE SETTING OUTPUT FLAGS

   ! DECLARE LOCAL VARIABLES
   ! These variables are ok on Nlevels only (not (Npoints, Nlevels) ) because we do not pass them directly to the COSP gridbox
   REAL, DIMENSION(Nlevels) :: eradc1d, eradi1d, eradr1d, erads1d, eradg1d   ! droplet, ice, snow, rain, graupel effective radii


   ! DATA ON HALF LEVELS FOR EACH POINT IN THE TILE
   REAL, DIMENSION(Npoints, Nlevels) :: p_half1d, ght_half                  ! pressure and geopotential height on half levels

   ! COLUMNS OF WRF VARIABLES FOR EACH POINT IN THE TILE
   REAL, DIMENSION(Npoints, Nlevels)   :: p1d, ght,                                    & ! pressure and geopotential height on full levels
                                          t1d, rh1d, sh1d,                             & ! temp, relative humidity, specific humidity
                                          ozonemr1d, cldfra1d,                         & ! ozone mixing ratio, cloud fraction, u- and v-wind components
                                          taucldc1d, taucldi1d, taucld1d,              & ! cloud droplet and ice optical depths, total cloud optical depth profile     
                                          qv1d, qc1d, qr1d, qi1d, qs1d, qg1d, qh1d,    & ! water vapour and hydrometeor mass mixing ratios
                                          rainr1d, hailr1d, snowr1d, graupelr1d,       & ! rain, hail, snow and graupel sedimentation rates 
                                          semiss                                         ! longwave cloud emissivity  


   ! TODO: remove the local variables that we don't actually use!

   ! COLUMNS OF 'CONVECTIVE' VARIABLES THAT WE DON'T HAVE (NOT USING A CONVECTION SCHEME) BUT NEED TO POPULATE FOR PASSING TO COSP
   REAL, DIMENSION(Npoints, Nlevels) ::   qv_c, qc_c, qr_c, qi_c, qs_c, qg_c, qh_c   & ! water vapour and hydrometeor mass mixing ratios, convection scheme
                                         ,semiss_c                                   & ! Cloud emissivity from convection scheme
                                         ,taucld_c, cca                              & ! cloud optical depth from convection scheme, convective cloud fraction
                                         ,rainr1d_c, hailr1d_c, snowr1d_c              ! rain, hail and snow sedimentation fluxes from convection scheme


   ! Array containing the effective radii of all the COSP hydrometeor types
   REAL, DIMENSION(Npoints, Nlevels, N_hydro) :: Reff


   ! SURFACE VALUES
   INTEGER, DIMENSION(Npoints) ::  sunlitbox, land   ! is gridbox sunlit? is gridbox a land point?
   REAL, DIMENSION(Npoints) ::    sfcp, skt, sht     ! surface pressure, temp, height
   REAL, DIMENSION(Npoints) :: xlat_ij, xlong_ij     ! latitude and longitude of given point
   REAL, DIMENSION(Npoints)  :: u10_ij, v10_ij        ! surface u-winf and v-wind
 
   ! Parameters for calculating relative humidity
   REAL, PARAMETER :: pq0 = 379.90516
   REAL, PARAMETER :: a2 = 17.2693882
   REAL, PARAMETER :: a3 = 273.16
   REAL, PARAMETER :: a4 = 35.86

   ! COUNTING AND INDEX VARIABLES
   INTEGER  :: i, j, k

   ! COSP DATA TYPES
   type(cosp_config)  :: cfg     ! Configuration options
   type(cosp_gridbox) :: gbx     ! Gridbox information. Input for COSP
   type(cosp_subgrid) :: sgx     ! Subgrid outputs
   type(cosp_sgradar) :: sgradar ! Output from radar simulator
   type(cosp_sglidar) :: sglidar ! Output from lidar simulator
   type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
   type(cosp_modis)   :: modis   ! Output from MODIS simulator
   type(cosp_misr)    :: misr    ! Output from MISR simulator
   type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
   type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
   type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
   if (lrttov_sim) then
      type(cosp_rttov)   :: rttov   ! Output from RTTOV 
   end if
   !---------- END OF VARIABLE DECLARATIONS -----------------------------------------------------------------------------

   !---------- 4. POPULATE COSP CONFIGURATION INPUT VARIABLE ("cfg") WITH SIMULATOR AND OUTPUT FLAGS --------------------
   cfg%Lradar_sim=lradar_sim
   cfg%Llidar_sim=llidar_sim
   cfg%Lisccp_sim=lisccp_sim
   cfg%Lmisr_sim=lmisr_sim
   cfg%Lmodis_sim=lmodis_sim
   cfg%Lrttov_sim=lrttov_sim
   cfg%Lalbisccp=lalbisccp
   cfg%Latb532=latb532
   cfg%Lboxptopisccp=lboxptopisccp
   cfg%Lboxtauisccp=lboxtauisccp
   cfg%Lcfaddbze94=lcfad_dbze94
   cfg%LcfadLidarsr532=lcfad_lidarsr532
   cfg%Lclcalipso=lclcalipso
   cfg%Lclhcalipso=lclhcalipso
   cfg%Lclisccp=lclisccp2
   cfg%Lcllcalipso=lcllcalipso
   cfg%Lclmcalipso=lclmcalipso
   cfg%Lcltcalipso=lcltcalipso
   cfg%Lpctisccp=lctpisccp
   cfg%Ldbze94=ldbze94
   cfg%Lcltradar=lcltradar
   cfg%Lcltradar2=lcltradar2
   cfg%Ltauisccp=ltauisccp
   cfg%Lcltisccp=ltclisccp
   cfg%Llongitude=llongitude
   cfg%Llatitude=llatitude
   cfg%LparasolRefl=lparasol_refl
   cfg%LclMISR=lclMISR
   cfg%Lmeantbisccp=lmeantbisccp
   cfg%Lmeantbclrisccp=lmeantbclrisccp
   cfg%Lclcalipso2=lclcalipso2
   cfg%Lcltlidarradar=lcltlidarradar
   cfg%Lfracout=lfrac_out
   cfg%LlidarBetaMol532=lbeta_mol532
   cfg%Lcltmodis=lcltmodis
   cfg%Lclwmodis=lclwmodis
   cfg%Lclimodis=lclimodis
   cfg%Lclhmodis=lclhmodis
   cfg%Lclmmodis=lclmmodis
   cfg%Lcllmodis=lcllmodis
   cfg%Ltautmodis=ltautmodis
   cfg%Ltauwmodis=ltauwmodis
   cfg%Ltauimodis=ltauimodis
   cfg%Ltautlogmodis=ltautlogmodis
   cfg%Ltauwlogmodis=ltauwlogmodis
   cfg%Ltauilogmodis=ltauilogmodis
   cfg%Lreffclwmodis=lreffclwmodis
   cfg%Lreffclimodis=lreffclimodis
   cfg%Lpctmodis=lpctmodis
   cfg%Llwpmodis=llwpmodis
   cfg%Liwpmodis=liwpmodis
   cfg%Lclmodis=lclmodis
   cfg%Ltbrttov=ltbrttov

   ! make sure COSP radar and lidar stats are calculated
   if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.

   cfg%Lwrite_output=.false.  !! COSP I/O not used (using WRF I/O!), so this should always be .false.

   !--------- END OF COSP CONFIGURATION --------------------------------------------------------------------------------

   ! --------- 5. SET THE COSP INPUTS -----------------------------------------------------------------------------------
   ! TODO - do we want to read these from a namelist...? (Likewise with the configuration and simulator flags)
   !
   LOGICAL :: use_vgrid = .true.  ! Use fixed vertical grid for outputs? (if .true. then you need to define number of levels with Nlr)
   LOGICAL :: csat_vgrid = .true. ! CloudSat vertical grid? (if .true. then the CloudSat standard grid is used for the outputs.
   !                              !  USE_VGRID needs also be .true.)
   INTEGER, PARAMETER ::  Npoints_it = 10000   ! Max # gridpoints to be processed in one iteration (10,000)
   INTEGER :: Nlr = 40   ! Number of levels in statistical outputs 
   DOUBLE PRECISION, PARAMETER :: emsfc_lw = 0.99 ! longwave emissivity of surface at 10.5 microns 
   REAL :: time   ! time since start of run [days], set to 1 because running over single WRF timestep
   DOUBLE PRECISION :: time_bnds(2) ! time boundaries - new in cosp v1.3, set following cosp_test.f90 line 121
   REAL :: time_step, toffset_step, half_time_step

   ! COSP INPUTS: RADAR
   REAL, PARAMETER :: radar_freq = 94.0  ! CloudSat radar frequency (GHz)
   INTEGER, PARAMETER :: surface_radar = 0 ! surface=1, spaceborne=0
   INTEGER, PARAMETER :: use_mie_tables = 0! use a precomputed lookup table? yes=1,no=0
   INTEGER, PARAMETER :: use_gas_abs = 1   ! include gaseous absorption? yes=1,no=0
   INTEGER, PARAMETER :: do_ray = 0        ! calculate/output Rayleigh refl=1, not=0
   INTEGER, PARAMETER :: melt_lay = 0      ! melting layer model off=0, on=1
   INTEGER, PARAMETER :: k2 = -1           ! |K|^2, -1=use frequency dependent default
   LOGICAL :: use_reff = .true. ! True if you want effective radius to be used by radar simulator (always used by lidar)
   LOGICAL :: use_precipitation_fluxes = .true.  ! True if precipitation fluxes are input to the algorithm    

   ! COSP INPUTS: LIDAR
   INTEGER, PARAMETER :: Nprmts_max_hydro = 12 ! Max number of parameters for hydrometeor size distributions
   INTEGER, PARAMETER :: Naero = 1             ! Number of aerosol species (Not used)
   INTEGER, PARAMETER :: Nprmts_max_aero = 1   ! Max number of parameters for aerosol size distributions (Not used)
   INTEGER, PARAMETER :: lidar_ice_type = 0    ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
   INTEGER, PARAMETER :: overlap = 3   !  overlap type: 1=max, 2=rand, 3=max/rand    

   ! COSP INPUTS: ISCCP
   INTEGER, PARAMETER :: isccp_topheight = 1  !  1 = adjust top height using both a computed
                       !  infrared brightness temperature and the visible
                       !  optical depth to adjust cloud top pressure. Note
                       !  that this calculation is most appropriate to compare
                       !  to ISCCP data during sunlit hours.
                       !  2 = do not adjust top height, that is cloud top
                       !  pressure is the actual cloud top pressure
                       !  in the model
                       !  3 = adjust top height using only the computed
                       !  infrared brightness temperature. Note that this
                       !  calculation is most appropriate to compare to ISCCP
                       !  IR only algortihm (i.e. you can compare to nighttime
                       !  ISCCP data with this option)
   INTEGER, PARAMETER ::  isccp_topheight_direction = 2   ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 ! 1 = default setting, and matches all versions of
                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                 ! 2 = experimental setting  
   ! COSP INPUTS: RRTOV
   INTEGER, PARAMETER ::   Platform = 1    ! satellite platform
   INTEGER, PARAMETER ::   Satellite = 15  ! satellite
   INTEGER, PARAMETER ::   Instrument = 0  ! instrument
   INTEGER, PARAMETER ::   Nchannels = 8   ! Number of channels to be computed
   INTEGER, DIMENSION(Nchannels) :: Channels
   DOUBLE PRECISION, DIMENSION(Nchannels) :: Surfem ! Surface emissivity
   REAL, PARAMETER ::  ZenAng = 50.0 ! Satellite Zenith Angle
   ! Mixing ratios of trace gases are below - taken from cosp_test.F90, units in kg/kg
   ! Could we get any of these from the radiation code instead?
   REAL, PARAMETER :: CO2 = 5.241e-04
   REAL, PARAMETER :: CH4 = 9.139e-07
   REAL, PARAMETER :: N2O = 4.665e-07
   REAL, PARAMETER :: CO = 2.098e-07

   ! Populate Channels and Surfem:
   Channels = (/1,3,5,6,8,10,11,13/)             !Channel numbers (please be sure that you supply Nchannels)
   Surfem = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)  ! Surface emissivity (please be sure that you supply Nchannels)

   ! Calculate time offsets of each point from the value in time
   ! If we are running COSP on each column, these parameters correctly make the timeoffset calculation later in gbx evaluate to zero!
   time_step = dt
   time = curr_secs
   toffset_step  = time_step/Npoints
   half_time_step = 0.5*time_step
   time_bnds = (/time-half_time_step,time+half_time_step/)


   ! -------- 6. LOOP OVER THE TILE, CONSTRUCT THE COSP INPUT AND CALL THE SIMULATOR -----------------------------------
      ! LOOP OVER I, J POINTS
      do j=jts,jte      ! j loop (north-south)
      do i=its,ite      ! i loop (east-west)

      ! CALCULATE WHETHER GRID POINT IS SUNLIT OR NOT
         if ( coszen(i,j) .gt. 0.2 ) then
               sunlitbox = 1
         else
               sunlitbox = 0
         end if

       ! only run COSP if the gridbox is sunlit!
       ! TODO: THIS IS ONLY THE CASE FOR SOME SIMULATORS! SOME WE CAN RUN REGARDLESS...
       if ( sunlitbox .eq. 1 ) then
          ! OFF WE GO...!   

         ! INITIALISE POINT VALUES
         sfcp(1:Npoints) = 0.
         skt(1:Npoints) = 0.
         sht(1:Npoints) = 0.
         land(1:Npoints) = 0.
         xlat_ij(1:Npoints) = 0.
         xlong_ij(1:Npoints) = 0.
         u10_ij(1:Npoints) = 0.
         v10_ij(1:Npoints) = 0.

         ! GET POINT VALUES FOR THIS POINT
         sfcp(1:Npoints) = psfc(i,j)      ! surface pressure
         skt(1:Npoints) =  tsk(i,j)       ! skin temperature
         sht(1:Npoints) =  ht(i,j)        ! surface height
         land(1:Npoints) = xland(i,j)     ! is this a land point?
         xlat_ij(1:Npoints) = xlat(i,j)   ! latitude of this point
         xlong_ij(1:Npoints) = xlong(i,j) ! longitude of this point
         u10_ij(1:Npoints) = u10(i,j)     ! surface u-wind
         v10_ij(1:Npoints) = v10(i,j)     ! surface v-wind

         ! INITIALISE 'CONVECTIVE' VALUES (NOT USING A CONVECTION SCHEME)
         semiss_c(1:Npoints,kts:kte) = 0.     ! cloud emissivity
         taucld_c(1:Npoints,kts:kte) = 0.     ! cloud optical depth
         cca(1:Npoints,kts:kte) = 0.          ! cloud fraction
         rainr1d_c(1:Npoints,kts:kte) = 0.    ! rain flux
         hailr1d_c(1:Npoints,kts:kte) = 0.    ! hail flux 
         snowr1d_c(1:Npoints,kts:kte) = 0.    ! snow flux
         graupelr1d_c(1:Npoints,kts:kte) = 0. ! graupel flux
         qv_c(1:Npoints,kts:kte) = 0.         ! water vapour mixing ratio
         qr_c(1:Npoints,kts:kte) = 0.         ! rain mixing ratio
         qs_c(1:Npoints,kts:kte) = 0.         ! snow mixing ratio
         qi_c(1:Npoints,kts:kte) = 0.         ! ice mixing ratio
         qg_c(1:Npoints,kts:kte) = 0.         ! graupel mixing ratio
         qh_c(1:Npoints,kts:kte) = 0.         ! hail mixing ratio 
         qc_c(1:Npoints,kts:kte) = 0.         ! cloud liquid mixing ratio
      
         ! INITIALISE COLUMN PROFILES
         p1d(1:Npoints,kts:kte) = 0.
         ght(1:Npoints,kts:kte) = 0.
         sh1d(1:Npoints,kts:kte) = 0.
         rh1d(1:Npoints,kts:kte) = 0.
         taucld1d(1:Npoints,kts:kte) = 0.
         semiss(1:Npoints,kts:kte) = 0.
         p_half1d(1:Npoints,kts:kte) = 0.
         ght_half(1:Npoints,kts:kte) = 0.
         t1d(1:Npoints,kts:kte) = 0.
         ozonemr1d(1:Npoints,kts:kte) = 0.
         cldfra1d(1:Npoints,kts:kte) = 0.
         eradc1d(kts:kte) = 0.
         eradi1d(kts:kte) = 0.
         eradr1d(kts:kte) = 0.
         erads1d(kts:kte) = 0.
         eradg1d(kts:kte) = 0.
         qv1d(1:Npoints,kts:kte) = 0.
         qc1d(1:Npoints,kts:kte) = 0.
         qr1d(1:Npoints,kts:kte) = 0.
         qi1d(1:Npoints,kts:kte) = 0.
         qs1d(1:Npoints,kts:kte) = 0.
         qg1d(1:Npoints,kts:kte) = 0.
         qh1d(1:Npoints,kts:kte) = 0.
         rainr1d(1:Npoints,kts:kte) = 0.
         hailr1d(1:Npoints,kts:kte) = 0.
         snowr1d(1:Npoints,kts:kte) = 0.
         graupelr1d(1:Npoints,kts:kte) = 0.

         ! INITIALISE REFF
         Reff(1:Npoints, 1:Nlevels, 1:N_Hydro) = 0.


         ! CALCULATE COLUMN PROFILES THAT AREN'T NATIVE TO WRF:
         do k = kts, kte     ! k loop (vertical)
             ! full pressure field
             p1d(1:Npoints,k) = pb(i,k,j) + p(i,k,j)
             ! geopotential height
             ght(1:Npoints,k) = ( phb(i,k,j) + ph(i,k,j) ) / G
             ! specific humidity
             sh1d(1:Npoints,k) = qv(i,k,j) / (1 + qv(i,k,j) )
             ! relative humidity
             rh1d(1:Npoints,k) = ( qv(i,k,j) / ( (pq0 / psfc(i,j)) * exp( a2*( t(i,k,j) - a3) / ( t(i,k,j) - a4) ) ) )*100.
             ! check: if relative humidity is greater than 1, set it to 1!
             IF ( rh1d(1:Npoints,k) .GT. 1. ) rh1d(k) = 1. ! BEWARE: THIS FAILS HORRIBLY IF NPOINTS NE 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT!
             ! cloud optical thickness
             taucld1d(1:Npoints,k) = taucldc(i,k,j) + taucldi(i,k,j)
              ! cloud emissivity = 1-exp(-cloud optical depth)
             semiss(1:Npoints,k) = 1.0 - exp(-taucld1d(1:Npoints,k)) ! BEWARE: THIS FAILS HORRIBLY IF NPOINTS NE 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT!
         end do

         ! LINEAR INTERPOLATION TO CALCULATE FIELDS ON HALF LEVELS:
         do k = kts, kte     ! k loop (vertical)
            if ( k .lt. kte ) then
              ! pressure on half levels
              p_half1d(1:Npoints,k) =  0.5*( p1d(1:Npoints,k) + p1d(1:Npoints,k+1) ) ! BEWARE: THIS FAILS HORRIBLY IF NPOINTS NE 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT!
              ! geopotential height on half levels 
              ght_half(1:Npoints,k) =  0.5*( ght(1:Npoints,k) + ght(1:Npoints,k+1) ) ! BEWARE: THIS FAILS HORRIBLY IF NPOINTS NE 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT! 
            else
              ! For k = kte, need to create the final point on half levels
              ! pressure on half levels
              p_half1d(1:Npoints,k) =  0.5*( 3*p1d(1:Npoints,k) - p1d(1:Npoints,k-1) )
              ! geopotential height on half levels
              ght_half(1:Npoints,k) =  0.5*( 3*ght(1:Npoints,k) - ght(1:Npoints,k-1) )
            end if
         end do

         ! TRANSFER NATIVE 3D ARRAYS INTO 1D COLUMNS FOR INPUT TO COSP:
         do k = kts, kte  
            t1d(1:Npoints,k) = t(i,k,j)                 ! temp
            ozonemr1d(1:Npoints,k) = ozonemr(i,k,j)     ! ozone mixing ratio
            cldfra1d(1:Npoints,k) = cldfra(i,k,j)       ! cloud fraction
            eradc1d(k) = eradc(i,k,j)                   ! liquid cloud effective radius
            eradi1d(k) = eradi(i,k,j)                   ! ice effective radius
            eradr1d(k) = eradr(i,k,j)                   ! rain drop effective radius
            erads1d(k) = erads(i,k,j)                   ! snow effective radius
            eradg1d(k) = eradg(i,k,j)                   ! graupel effective radius
            ! The  mixing ratios can be assigned as follows because we use Npoints = 1
            qv1d(1:Npoints,k) = qv(i,k,j)               ! water vapour mixing ratio 
            qc1d(1:Npoints,k) = qc(i,k,j)               ! liquid cloud mixing ratio
            qr1d(1:Npoints,k) = qr(i,k,j)               ! rain mixing ratio
            qi1d(1:Npoints,k) = qi(i,k,j)               ! ice mixing ratio
            qs1d(1:Npoints,k) = qs(i,k,j)               ! snow mixing ratio
            qg1d(1:Npoints,k) = qg(i,k,j)               ! graupel mixing ratio
            qh1d(1:Npoints,k) = qh(i,k,j)               ! hail mixing ratio  
            rainr1d(1:Npoints,k) = rainncv(i,k,j)       ! rain flux
            hailr1d(1:Npoints,k) = hailncv(i,k,j)       ! hail flux
            snowr1d(1:Npoints,k) = snowncv(i,k,j)       ! snow flux
            graupelr1d(1:Npoints,k) = graupelncv(i,k,j) ! graupel flux
         end do

         ! CONVERT EFFECTIVE RADII FROM MICRONS (WRF) TO METRES (WHAT COSP WANTS)
         ! NOTE: WE USE THE SAME VALUES FOR THE CONVECTIVE AND LARGE-SCALE EFFECTIVE RADII
         ! TODO: FIND OUT WHAT WRF'S FILLVALUE IS - NEED TO RESET IN REFF SINCE HAVE MULTIPLIED BY 1.e-6!
         Reff(1:Npoints,1:Nlevels,I_LSCLIQ) = eradc1d*1.e-6  ! LSCLIQ
         Reff(1:Npoints,1:Nlevels,I_LSCICE) = eradi1d*1.e-6  ! LSCICE
         Reff(1:Npoints,1:Nlevels,I_LSRAIN) = eradr1d*1.e-6  ! LSRAIN
         Reff(1:Npoints,1:Nlevels,I_LSSNOW) = erads1d*1.e-6  ! LSSNOW
         Reff(1:Npoints,1:Nlevels,I_CVCLIQ) = eradc1d*1.e-6  ! CVCLIQ - USE SAME VALUES AS FOR LSCLIQ
         Reff(1:Npoints,1:Nlevels,I_CVCICE) = eradi1d*1.e-6  ! CVCICE - USE SAME VALUES AS FOR LSCICE
         Reff(1:Npoints,1:Nlevels,I_CVRAIN) = eradr1d*1.e-6  ! CVRAIN -  USE SAME VALUES AS FOR LSRAIN
         Reff(1:Npoints,1:Nlevels,I_CVSNOW) = erads1d*1.e-6  ! CVSNOW -  USE SAME VALUES AS FOR LSSNOW
         Reff(1:Npoints,1:Nlevels,I_LSGRPL) = eradg1d*1.e-6  ! LSGRPL     

         !----------------------------------------------------------------------
         ! Allocate memory for gridbox type
         !----------------------------------------------------------------------
         call construct_cosp_gridbox(time,time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
                                    do_ray,melt_lay,k2, &
                                    Npoints,Nlevels,Ncolumns,N_hydro,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                                    lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
                                    use_precipitation_fluxes,use_reff, &
                                    Platform,Satellite,Instrument,Nchannels,ZenAng, &
                                    channels(1:Nchannels),surfem(1:Nchannels),CO2,CH4,N2O,CO,gbx) ! gbx is output
        !-----------------------------------------------------------------------
        ! Code to populate input structure
        !-----------------------------------------------------------------------
        gbx%longitude = xlong_ij 
        gbx%latitude = xlat_ij
        ! Toffset. This assumes that time is the mid-point of the interval.
        ! Sanity check: this does equate to zero for Npoints = 1!
        do n=1,Npoints
          gbx%toffset(n) = -half_time_step + toffset_step*(n-0.5)
        end do
        gbx%p = p1d ! pressure on full levels
        gbx%ph = p_half1d ! pressure on half levels
        gbx%zlev = ght ! height on full levels
        gbx%zlev_half = ght_half ! height on half levels
        gbx%T = t1d ! temperature
        gbx%q = rh1d ! relative humidity
        gbx%sh = sh1d ! specific humidity 
        gbx%cca = cca ! convective cloud amount
        gbx%tca = cldfra1d ! total cloud amount
        gbx%psfc = sfcp ! surface pressure
        gbx%skt  = skt ! skin temperature
        gbx%land = land ! landmask
        gbx%mr_ozone  = ozonemr1d ! ozone mass mixing ratio
        gbx%u_wind  = u10_ij ! surface u-wind
        gbx%v_wind  = v10-ij ! surface v-wind
        gbx%sunlit  = 1 ! we already know we are only running the code if sunlitbox = 1
        gbx%mr_hydro(:,:,I_LSCLIQ) = qc1d ! mixing ratio large-scale cloud
        gbx%mr_hydro(:,:,I_LSCICE) = qi1d ! mixing ratio large-scale ice
        gbx%mr_hydro(:,:,I_CVCLIQ) = qc_c ! mixing ratio convective cloud
        gbx%mr_hydro(:,:,I_CVCICE) = qi_c ! mixing ratio convective ice
        gbx%rain_ls = rainr1d ! large-scale rain flux
        gbx%snow_ls = snowr1d ! large-scale snow flux
        gbx%grpl_ls = graupelr1d ! large-scale graupel flux
        gbx%rain_cv = rainr1d_c ! convective rain flux
        gbx%snow_cv = snowr1d_c ! convective snow flux
        gbx%Reff = Reff

        ! ISCCP simulator
        gbx%dtau_s   = taucld1d ! stratiform cloud optical depth
        gbx%dtau_c   = taucld_c ! convective cloud optical depth
        gbx%dem_s    = semiss   ! stratiform cloud emissivity
        gbx%dem_c    = semiss_c ! convective cloud emissivity

        ! Note: COSP output variables are sgx (sub-grid outputs), sgradar (radar outputs), sglidar (lidar outputs),
        ! isccp (isccp outputs), misr (misr simulator outputs), vgrid (vertical grid info), stradar
        ! (summary statistics radar simulator), stlidar (summary statistics lidar simulator)


        !-----------------------------------------------------------------------
        ! Define new vertical grid
        !-----------------------------------------------------------------------
        call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)

        !-----------------------------------------------------------------------
        ! Allocate memory for other types
        !-----------------------------------------------------------------------
        call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
        call construct_cosp_sgradar(cfg,Npoints,Ncolumns,Nlevels,N_hydro,sgradar)
        call construct_cosp_radarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_hydro,stradar)
        call construct_cosp_sglidar(cfg,Npoints,Ncolumns,Nlevels,N_hydro,PARASOL_NREFL,sglidar)
        call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_hydro,PARASOL_NREFL,stlidar)
        call construct_cosp_isccp(cfg,Npoints,Ncolumns,Nlevels,isccp)
        call construct_cosp_modis(cfg,Npoints,modis)
        call construct_cosp_misr(cfg,Npoints,misr)
        if (lrttov_sim) then
           call construct_cosp_rttov(Npoints,Nchannels,rttov)
        end if

        !-----------------------------------------------------------------------
        ! CALL SIMULATOR
        !-----------------------------------------------------------------------
        ! TODO: DO WE CALL COSP OR cosp_simulator HERE? COSP CALLED FROM THE TEST MODULE, BUT COSP.F90 CALLS THE SIMULATOR...
        if (lrttov_sim) then
           call cosp_simulator(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
        else 
           call cosp_simulator(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
        end if

        ! TODO: TRANSFER OUTPUT TO 3D ARRAYS AND WRITE OUTPUT TO WRF VARIABLES
        ! TODO: DEALLOCATE MEMORY IN DERIVED TYPES

      end if ! only run COSP if gridbox is sunlit TODO: THIS ONLY APPLIES FOR SOME OF THE SIMULATORS! SEE CAM IMPLEMENTATION

      end do ! i-point loop
      end do ! j-point loop

      ! TESTING WRITING TO NEW VARIABLES!!!
      eradc(ims:ime, kms:kme, jms:jme)=-15
      eradi(ims:ime, kms:kme, jms:jme)=-15
      erads(ims:ime, kms:kme, jms:jme)=-15
      eradr(ims:ime, kms:kme, jms:jme)=-15
      eradg(ims:ime, kms:kme, jms:jme)=-15
      cldfra(ims:ime, kms:kme, jms:jme)=-15

     END DO ! tiles loop
   end subroutine wrf_cosp_simulator
END MODULE module_cosp_simulator_test
!END MODULE module_cosp_wrf
