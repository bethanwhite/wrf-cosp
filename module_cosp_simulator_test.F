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
 !
 ! 
 ! - STRUCTURE OF OUTPUTS:
 ! - pass as array to subroutine
 ! - declare as allocatable - because we need to calculate some dimensions (even though these are set in the Registry
 !   for writing out, I don't think we can access that information from here...
 ! - allocate: PROBLEMS HERE! allocate tile dims or memory dims? local or global? where do we calculate kts:kte? etc (see code)
 ! - initialise all to 0
 ! - loop over the tile (its:ite), call simulator and write output
 ! - deallocate allocatable memory for the output before the tile loop re-runs
 ! - UPDATE: kts, kte are passed in from the dynamics solver (solve_em.f) therefore ok to use
 !-------------------------------------------------------------------------------------------------
 
 ! USE THE COSP MODULES
 USE mod_cosp_constants    !, ONLY: R_UNDEF, PARASOL_NREFL, i_cvcice, i_cvcliq, i_lscliq, i_lscice ! NOTE - actually want to read all constants!
 USE mod_cosp_types
 USE mod_cosp_simulator
 USE mod_cosp_modis_simulator


 ! USE WRF ERROR MODULE?
 !USE module_wrf_error


 ! USE WRF PHYSICS CONSTANTS:
  use module_model_constants, ONLY: &
                                    G, &        ! gravitational acceleration
                                    R => r_d    ! dry gas constant

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
                      ,dt, curr_secs                                     &
                      ,COSP_CLLCALIPSO, COSP_CLMCALIPSO, COSP_CLHCALIPSO       & ! COSP OUTPUTS
                      ,COSP_CLTCALIPSO, COSP_CLTLIDARRADAR, COSP_CLTISCCP      &
                      ,COSP_PCTISCCP, COSP_TAUISCCP, COSP_ALBISCCP             &
                      ,COSP_MEANTBCLRISCCP, COSP_MEANTBISCCP, COSP_CLTMODIS    &
                      ,COSP_CLWMODIS, COSP_CLIMODIS, COSP_CLHMODIS             &
                      ,COSP_CLMMODIS, COSP_CLLMODIS, COSP_TAUTMODIS            &
                      ,COSP_TAUWMODIS, COSP_TAUIMODIS, COSP_TAUTLOGMODIS       &
                      ,COSP_TAUWLOGMODIS, COSP_TAUILOGMODIS, COSP_REFFCLWMODIS &
                      ,COSP_REFFCLIMODIS, COSP_PCTMODIS, COSP_LWPMODIS         &
                      ,COSP_IWPMODIS, COSP_TOFFSET, COSP_BOXTAUISCCP           &
                      ,COSP_BOXTOPISCCP, COSP_CLCALIPSO, COSP_CLCALIPSO2       &
                      ,COSP_LIDARBETAMOL532                                    &
                      ,COSP_PARASOLREFL, COSP_DBZE94                           &
                      ,COSP_ATB532, COSP_FRACOUT, COSP_CFADDBZE94              &
                      ,COSP_CFADLIDARSR532, COSP_CLISCCP, COSP_CLMODIS         &
                      ,COSP_CLMISR )


   !------------ 1. VARIABLE DECLARATIONS: ------------------------------------------------------------------------

   ! INPUT VARIABLES -----------------------------------------------------------------------------------------------
   LOGICAL,INTENT(IN)    :: specified
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradc
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradi
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::erads
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradr
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::eradg
   REAL,INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme)::cldfra
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::ozonemr
   REAL,INTENT(IN), DIMENSION(ims:ime,jms:jme)::u10
   REAL,INTENT(IN), DIMENSION(ims:ime,jms:jme)::v10
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::tsk
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::ht
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qv
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qc
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qr
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qi
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qs
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qg
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::qh
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::phb
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::ph
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::pb
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::p
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::psfc
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::t
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::xland
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::taucldc
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::taucldi
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::emiss
   REAL,INTENT(IN), DIMENSION(ims:ime, jms:jme)::coszen
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::rainncv
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::hailncv
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::snowncv
   REAL,INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::graupelncv

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

   ! Number of gridpoints
   ! BEWARE: CODE BELOW DEPENDS ON NPOINTS = 1 BECAUSE WE HAVE CHOSEN TO MOVE DATA ONTO THE COLUMN AND RUN COSP ONCE FOR EACH COLUMN
   ! IF YOU CHANGE NPOINTS, BE SURE TO FOLLOW THROUGH EVERYTHING ELSE CAREFULLY!
   INTEGER, PARAMETER :: Npoints = 1 ! COSP is called once for each column
   ! Number of levels
   INTEGER :: Nlevels
   ! Number of columns
   INTEGER, PARAMETER :: Ncolumns = Npoints
   ! Number of levels in statistical outputs
   INTEGER, PARAMETER :: Nlr = 40   

   ! OUTPUT VARIABLES -------------------------------------------------------------------------------------------------
   ! 2D spatial arrays
   REAL,INTENT(INOUT), DIMENSION(ims:ime, jms:jme) ::                         &
      COSP_CLLCALIPSO, COSP_CLMCALIPSO, COSP_CLHCALIPSO, COSP_CLTCALIPSO,     &
      COSP_CLTLIDARRADAR, COSP_CLTISCCP, COSP_PCTISCCP, COSP_TAUISCCP,        &
      COSP_ALBISCCP, COSP_MEANTBCLRISCCP, COSP_MEANTBISCCP, COSP_CLTMODIS,    &
      COSP_CLWMODIS, COSP_CLIMODIS, COSP_CLHMODIS, COSP_CLMMODIS,             &
      COSP_CLLMODIS, COSP_TAUTMODIS, COSP_TAUWMODIS, COSP_TAUIMODIS,          &
      COSP_TAUTLOGMODIS, COSP_TAUWLOGMODIS, COSP_TAUILOGMODIS,                &
      COSP_REFFCLWMODIS, COSP_REFFCLIMODIS, COSP_PCTMODIS, COSP_LWPMODIS,     &
      COSP_IWPMODIS, COSP_TOFFSET, COSP_BOXTAUISCCP, COSP_BOXTOPISCCP

   ! 3D spatial output arrays on stat levels, Nlr
   REAL,INTENT(INOUT), DIMENSION(ims:ime, Nlr, jms:jme) ::                    &
      COSP_CLCALIPSO, COSP_CLCALIPSO2

   ! 3D spatial arrays on model levels, Nlevels = kms-kme
   REAL,INTENT(INOUT), DIMENSION(ims:ime, 1:kms-kme, jms:jme) ::                &
      COSP_LIDARBETAMOL532, COSP_DBZE94, COSP_ATB532, COSP_FRACOUT

   ! 4D arrays
   REAL,INTENT(INOUT), DIMENSION(ims:ime,jms:jme,PARASOL_NREFL) :: COSP_PARASOLREFL
   REAL,INTENT(INOUT), DIMENSION(ims:ime,1:Nlr,jms:jme,DBZE_BINS) :: COSP_CFADDBZE94
   REAL,INTENT(INOUT), DIMENSION(ims:ime,1:Nlr,jms:jme,SR_BINS) :: COSP_CFADLIDARSR532
   REAL,INTENT(INOUT), DIMENSION(ims:ime,jms:jme,7,7) ::  COSP_CLISCCP, COSP_CLMODIS
   REAL,INTENT(INOUT), DIMENSION(ims:ime,jms:jme,7,MISR_N_CTH) :: COSP_CLMISR


   !-------- LOCAL VARIABLES -------------------------------------------------------------------------------------------
   INTEGER :: its,ite,jts,jte,ij,sz
   LOGICAL :: channel

   !--------- SIMULATOR FLAGS - SET WHICH SIMULATORS TO RUN HERE --------------------------------------------------------
   LOGICAL :: lradar_sim = .true.   ! Radar
   LOGICAL :: llidar_sim = .false.   ! Lidar
   LOGICAL :: lmodis_sim = .false.   ! MODIS
   LOGICAL :: lmisr_sim = .false.    ! MISR
   LOGICAL :: lisccp_sim = .false.   ! ISCCP
   LOGICAL :: lrttov_sim = .false.   ! RRTOV

   !--------- DO YOU WANT SUB-COLUMN OUTPUT ? ---------------------------------------------------------------------------
   LOGICAL :: lfrac_out = .false.

   !--------- COMPUTE STATS? --------------------------------------------------------------------------------------------
   LOGICAL :: Lstats = .false.

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

   ! NEVER USE COSP I/O! (USING WRF I/O). THIS SHOULD ALWAYS BE .FALSE.
   LOGICAL :: Lwrite_output = .false.

   ! TODO: remove the local variables that we don't actually use!
   ! COSP INPUT VARIABLES
   REAL, DIMENSION(1:kms-kme) ::                             &
      eradc1d, eradi1d, eradr1d, erads1d, eradg1d                          ! droplet, ice, snow, rain, graupel effective radii

   REAL, DIMENSION(Npoints, 1:kms-kme) ::                    &
      p_half1d, ght_half, p1d, ght,                                        & ! pressure and geopotential height on half and full levels
      t1d, rh1d, sh1d, ozonemr1d, cldfra1d,                                & ! temp, relative humidity, specific humidity, ozone mixing ratio, cloud fraction
      taucldc1d, taucldi1d, taucld1d,                                      & ! cloud droplet and ice optical depths, total cloud optical depth profile
      qv1d, qc1d, qr1d, qi1d, qs1d, qg1d, qh1d,                            & ! water vapour and hydrometeor mass mixing ratios
      rainr1d, hailr1d, snowr1d, graupelr1d, semiss,                       & ! rain, hail, snow and graupel sedimentation rates, longwave cloud emissivity
      qv_c, qc_c, qr_c, qi_c, qs_c, qg_c, qh_c,                            & ! water vapour and hydrometeor mass mixing ratios, convection scheme
      semiss_c, taucld_c, cca,                                             & ! cloud emissivity, cloud optical depth, convective cloud fraction from convection scheme
      rainr1d_c, snowr1d_c                                                   ! rain and snow sedimentation fluxes from convection scheme

   REAL, DIMENSION(Npoints, 1:kms-kme, N_Hydro) ::  Reff         ! effective radii of all the COSP hydrometeor types

   REAL, DIMENSION(Npoints) ::                               & 
      sunlitbox, land, sfcp, skt, sht, xlat_ij, xlong_ij, u10_ij, v10_ij     ! is gridbox sunlit? is gridbox a land point? surface pressure, temp, height
                                                                             ! latitude and longitude of given point, 10m u-wind and v-wind

   ! Parameters for calculating relative humidity
   REAL, PARAMETER :: pq0 = 379.90516
   REAL, PARAMETER :: a2 = 17.2693882
   REAL, PARAMETER :: a3 = 273.16
   REAL, PARAMETER :: a4 = 35.86

   ! COUNTING AND INDEX VARIABLES
   INTEGER  :: i, j, k, n

   ! COSP DATA TYPES
   type(cosp_config)  :: cfg        ! Configuration options
   type(cosp_gridbox) :: gbx        ! Gridbox input to COSP
   type(cosp_subgrid) :: sgx        ! Subgrid inputs to COSP
   type(cosp_sgradar) :: sgradar    ! Output from radar simulator
   type(cosp_sglidar) :: sglidar    ! Output from lidar simulator
   type(cosp_isccp)   :: isccp      ! Output from ISCCP simulator
   type(cosp_modis)   :: modis      ! Output from MODIS simulator
   type(cosp_misr)    :: misr       ! Output from MISR simulator
   type(cosp_sghydro) :: sghydro    ! Subgrid info for hydrometeors
   type(cosp_vgrid)   :: vgrid      ! Information on vertical grid of stats
   type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
   type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator


   ! COSP INPUT FLAGS AND PARAMETERS
   ! TODO - do we want to read these from a namelist...? (Likewise with the configuration and simulator flags)

   LOGICAL :: use_vgrid = .true.  ! Use fixed vertical grid for outputs? (if .true. then you need to define number of levels with Nlr)
   LOGICAL :: csat_vgrid = .true. ! CloudSat vertical grid? (if .true. then the CloudSat standard grid is used for the outputs.
   !                              !  USE_VGRID needs also be .true.)
   INTEGER, PARAMETER ::  Npoints_it = 10000   ! Max # gridpoints to be processed in one iteration (10,000)
   REAL, PARAMETER :: emsfc_lw = 0.99 ! longwave emissivity of surface at 10.5 microns
   DOUBLE PRECISION :: time   ! time since start of run [days], set to 1 because running over single WRF timestep
   DOUBLE PRECISION :: time_bnds(2) ! time boundaries - new in cosp v1.3, set following cosp_test.f90 line 121
   REAL :: time_step, toffset_step, half_time_step

   ! COSP INPUTS: RADAR
   REAL, PARAMETER :: radar_freq = 94.0  ! CloudSat radar frequency (GHz)
   INTEGER, PARAMETER :: surface_radar = 0 ! surface=1, spaceborne=0
   INTEGER, PARAMETER :: use_mie_tables = 0! use a precomputed lookup table? yes=1,no=0
   INTEGER, PARAMETER :: use_gas_abs = 1   ! include gaseous absorption? yes=1,no=0
   INTEGER, PARAMETER :: do_ray = 0        ! calculate/output Rayleigh refl=1, not=0
   INTEGER, PARAMETER :: melt_lay = 0      ! melting layer model off=0, on=1
   REAL, PARAMETER :: k2 = -1           ! |K|^2, -1=use frequency dependent default
   LOGICAL :: use_reff = .true. ! True if you want effective radius to be used by radar simulator (always used by lidar)
   LOGICAL :: use_precipitation_fluxes = .true.  ! True if precipitation fluxes are input to the algorithm

   ! COSP INPUTS: LIDAR
   INTEGER, PARAMETER :: Nprmts_max_hydro = 12 ! Max number of parameters for hydrometeor size distributions
   INTEGER, PARAMETER :: Naero = 1             ! Number of aerosol species (Not used)
   INTEGER, PARAMETER :: Nprmts_max_aero = 1   ! Max number of parameters for aerosol size distributions (Not used)
   INTEGER, PARAMETER :: lidar_ice_type = 0    ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
   INTEGER, PARAMETER :: overlap = 3   !  overlap type: 1=max, 2=rand, 3=max/rand
   INTEGER, PARAMETER :: Nrefl = 5   ! number of parasol reflectors

   ! COSP INPUTS: ISCCP
   INTEGER, PARAMETER :: isccp_topheight = 1
   ! 1 = adjust top height using both a computed infrared brightness temperature and the visible optical depth to adjust cloud top pressure.
   ! Note that this calculation is most appropriate to compare to ISCCP data during sunlit hours. 
   ! 2 = do not adjust top height, that is cloud top pressure is the actual cloud top pressure in the model
   ! 3 = adjust top height using only the computed infrared brightness temperature. Note that this calculation is most appropriate to compare
   ! to ISCCP IR only algortihm (i.e. you can compare to nighttime ISCCP data with this option)

   INTEGER, PARAMETER ::  isccp_topheight_direction = 2
   ! direction for finding atmosphere pressure level with interpolated temperature equal to the radiance determined cloud-top temperature
   ! 1 = find the *lowest* altitude (highest pressure) level with interpolated temperature equal to the radiance determined cloud-top temperature
   ! 2 = find the *highest* altitude (lowest pressure) level with interpolated temperature equal to the radiance determined cloud-top temperature
   ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
   ! 1 = default setting, and matches all versions of ISCCP simulator with versions numbers 3.5.1 and lower
   ! 2 = experimental setting

   ! COSP INPUT PARAMETERS: RRTOV
   INTEGER, PARAMETER ::   Platform = 1    ! satellite platform
   INTEGER, PARAMETER ::   Satellite = 15  ! satellite
   INTEGER, PARAMETER ::   Instrument = 0  ! instrument
   INTEGER, PARAMETER ::   Nchannels = 8   ! Number of channels to be computed
   INTEGER, PARAMETER ::   Channels(Nchannels) = (/1,3,5,6,8,10,11,13/) !Channel numbers (please be sure that you supply Nchannels)
   REAL, PARAMETER :: Surfem(NChannels) = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)  ! Surface emissivity (please be sure that you supply Nchannels)
   REAL, PARAMETER ::  ZenAng = 50.0 ! Satellite Zenith Angle
   ! Mixing ratios of trace gases are below - taken from cosp_test.F90, units in kg/kg
   ! Could we get any of these from the radiation code instead?
   REAL, PARAMETER :: CO2 = 5.241e-04
   REAL, PARAMETER :: CH4 = 9.139e-07
   REAL, PARAMETER :: N2O = 4.665e-07
   REAL, PARAMETER :: CO = 2.098e-07


   !---------- END OF VARIABLE DECLARATIONS -----------------------------------------------------------------------------

   !---------- CALCULATE DIMENSIONS -----------------------------------------------------------------------------
   ! Calculate Nlevels
    Nlevels = kte - kts

   ! Calculate time offsets of each point from the value in time
   ! If we are running COSP on each column, these parameters correctly make the timeoffset calculation later in gbx evaluate to zero!

   time_step = dt
   time = curr_secs
   toffset_step  = time_step/Npoints
   half_time_step = 0.5*time_step
   time_bnds = (/time-half_time_step,time+half_time_step/)


   !---------- 2. INITIALISE OUTPUT ARRAYS -------------------------------------------------------------------------------
   COSP_CLLCALIPSO(ims:ime,jms:jme) = 0.
   COSP_CLMCALIPSO(ims:ime,jms:jme) = 0.
   COSP_CLHCALIPSO(ims:ime,jms:jme) = 0.
   COSP_CLTCALIPSO(ims:ime,jms:jme) = 0.
   COSP_CLTLIDARRADAR(ims:ime,jms:jme) = 0.
   COSP_CLTISCCP(ims:ime,jms:jme) = 0.
   COSP_PCTISCCP(ims:ime,jms:jme) = 0.
   COSP_TAUISCCP(ims:ime,jms:jme) = 0.
   COSP_ALBISCCP(ims:ime,jms:jme) = 0.
   COSP_MEANTBCLRISCCP(ims:ime,jms:jme) = 0.
   COSP_MEANTBISCCP(ims:ime,jms:jme) = 0.
   COSP_CLTMODIS(ims:ime,jms:jme) = 0.
   COSP_CLWMODIS(ims:ime,jms:jme) = 0.
   COSP_CLIMODIS(ims:ime,jms:jme) = 0.
   COSP_CLHMODIS(ims:ime,jms:jme) = 0.
   COSP_CLMMODIS(ims:ime,jms:jme) = 0.
   COSP_CLLMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUTMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUWMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUIMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUTLOGMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUWLOGMODIS(ims:ime,jms:jme) = 0.
   COSP_TAUILOGMODIS(ims:ime,jms:jme) = 0.
   COSP_REFFCLWMODIS(ims:ime,jms:jme) = 0.
   COSP_REFFCLIMODIS(ims:ime,jms:jme) = 0.
   COSP_PCTMODIS(ims:ime,jms:jme) = 0.
   COSP_LWPMODIS(ims:ime,jms:jme) = 0.
   COSP_IWPMODIS(ims:ime,jms:jme) = 0.
   COSP_TOFFSET(ims:ime,jms:jme) = 0.
   COSP_BOXTAUISCCP(ims:ime,jms:jme) = 0.
   COSP_BOXTOPISCCP(ims:ime,jms:jme) = 0.
   COSP_CLCALIPSO(ims:ime,1:Nlr,jms:jme) = 0.
   COSP_CLCALIPSO2(ims:ime,1:Nlr,jms:jme) = 0.
   COSP_LIDARBETAMOL532(ims:ime,1:Nlevels,jms:jme) = 0.
   COSP_PARASOLREFL(ims:ime,jms:jme,1:PARASOL_NREFL) = 0.
   COSP_DBZE94(ims:ime,1:Nlevels,jms:jme) = 0.
   COSP_ATB532(ims:ime,1:Nlevels,jms:jme) = 0.
   COSP_FRACOUT(ims:ime,1:Nlevels,jms:jme) = 0.
   COSP_CFADDBZE94(ims:ime,1:Nlr,jms:jme,1:DBZE_BINS) = 0.
   COSP_CFADLIDARSR532(ims:ime,1:Nlr,jme:jme,1:SR_BINS) = 0.
   COSP_CLISCCP(ims:ime,jms:jme,1:7,1:7) = 0.
   COSP_CLMODIS(ims:ime,jms:jme,1:7,1:7) = 0.
   COSP_CLMISR(ims:ime,jms:jme,1:7,1:MISR_N_CTH) = 0.


   ! TODO: remove any local input variables that we don't actually use!

   !--------- 3. WE NEED TO CALCULATE WRF TILE DIMENSIONS -----------------
   !--------- THE FOLLOWING CODE FOLLOWS THE DEFINITION OF ITS, ITE, JTS, JTE FROM THE MICROPHYSICS DRIVER --------------
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

   !---------- 4. SET FLAGS --------------------------------------------------------------------------------------------
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

   ! make sure COSP radar and lidar stats are calculated
   if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) Lstats = .true.

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
 

   !---------- 5. POPULATE COSP CONFIGURATION INPUT VARIABLE ("cfg") WITH SIMULATOR AND OUTPUT FLAGS --------------------
   ! NOTE: Ltoffset is missing from this list - we don't set it and don't copy it to cfg structure. Do we need it?
   cfg%Lradar_sim=lradar_sim
   cfg%Llidar_sim=llidar_sim
   cfg%Lisccp_sim=lisccp_sim
   cfg%Lmisr_sim=lmisr_sim
   cfg%Lmodis_sim=lmodis_sim
   cfg%Lrttov_sim=lrttov_sim
   cfg%Lstats = Lstats
   cfg%Lwrite_output = Lwrite_output
   cfg%Lalbisccp=lalbisccp
   cfg%Latb532=latb532
   cfg%Lboxptopisccp=lboxptopisccp
   cfg%Lboxtauisccp=lboxtauisccp
   cfg%Lcfaddbze94=lcfad_dbze94
   cfg%LcfadLidarsr532=lcfad_lidarsr532
   cfg%Lclcalipso=lclcalipso
   cfg%Lclcalipso2=lclcalipso2
   cfg%Lclhcalipso=lclhcalipso
   cfg%Lclisccp=lclisccp2
   cfg%Lcllcalipso=lcllcalipso
   cfg%Lclmcalipso=lclmcalipso
   cfg%Lcltcalipso=lcltcalipso
   cfg%Lcltlidarradar=lcltlidarradar
   cfg%Lpctisccp=lctpisccp
   cfg%Ldbze94=ldbze94
   cfg%Ltauisccp=ltauisccp
   cfg%Lcltisccp=ltclisccp
   cfg%LparasolRefl=lparasol_refl
   cfg%LclMISR=lclmisr
   cfg%Lmeantbisccp=lmeantbisccp
   cfg%Lmeantbclrisccp=lmeantbclrisccp
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

   !--------- END OF COSP CONFIGURATION --------------------------------------------------------------------------------

   ! -------- 6. LOOP OVER THE TILE, CONSTRUCT THE COSP INPUT AND CALL THE SIMULATOR -----------------------------------
      ! LOOP OVER I, J POINTS
      do j=jts,jte      ! j loop (north-south)
      do i=its,ite      ! i loop (east-west)

      ! CALCULATE WHETHER GRID POINT IS SUNLIT OR NOT
      ! VARIABLE 'SUNLITBOX' HAS DIMENSIONS NPOINTS (REQUIRED FOR PASSING TO GBX), BUT WE HAVE NPOINTS = 1
      ! THEREFORE WE CALCULATE SUNLITBOX(1)
         if ( coszen(i,j) .gt. 0.2 ) then
               sunlitbox(1) = 1
         else
               sunlitbox(1) = 0
         end if

       ! only run COSP if the gridbox is sunlit!
       ! TODO: THIS IS ONLY THE CASE FOR SOME SIMULATORS! SOME WE CAN RUN REGARDLESS...
       if ( sunlitbox(1) .eq. 1 ) then
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
         semiss_c(1:Npoints,1:Nlevels) = 0.     ! cloud emissivity
         taucld_c(1:Npoints,1:Nlevels) = 0.     ! cloud optical depth
         cca(1:Npoints,1:Nlevels) = 0.          ! cloud fraction
         rainr1d_c(1:Npoints,1:Nlevels) = 0.    ! rain flux
         snowr1d_c(1:Npoints,1:Nlevels) = 0.    ! snow flux
         qv_c(1:Npoints,1:Nlevels) = 0.         ! water vapour mixing ratio
         qr_c(1:Npoints,1:Nlevels) = 0.         ! rain mixing ratio
         qs_c(1:Npoints,1:Nlevels) = 0.         ! snow mixing ratio
         qi_c(1:Npoints,1:Nlevels) = 0.         ! ice mixing ratio
         qg_c(1:Npoints,1:Nlevels) = 0.         ! graupel mixing ratio
         qh_c(1:Npoints,1:Nlevels) = 0.         ! hail mixing ratio 
         qc_c(1:Npoints,1:Nlevels) = 0.         ! cloud liquid mixing ratio
      
         ! INITIALISE COLUMN PROFILES
         p1d(1:Npoints,1:Nlevels) = 0.
         ght(1:Npoints,1:Nlevels) = 0.
         sh1d(1:Npoints,1:Nlevels) = 0.
         rh1d(1:Npoints,1:Nlevels) = 0.
         taucld1d(1:Npoints,1:Nlevels) = 0.
         semiss(1:Npoints,1:Nlevels) = 0.
         p_half1d(1:Npoints,1:Nlevels) = 0.
         ght_half(1:Npoints,1:Nlevels) = 0.
         t1d(1:Npoints,1:Nlevels) = 0.
         ozonemr1d(1:Npoints,1:Nlevels) = 0.
         cldfra1d(1:Npoints,1:Nlevels) = 0.
         eradc1d(1:Nlevels) = 0.
         eradi1d(1:Nlevels) = 0.
         eradr1d(1:Nlevels) = 0.
         erads1d(1:Nlevels) = 0.
         eradg1d(1:Nlevels) = 0.
         qv1d(1:Npoints,1:Nlevels) = 0.
         qc1d(1:Npoints,1:Nlevels) = 0.
         qr1d(1:Npoints,1:Nlevels) = 0.
         qi1d(1:Npoints,1:Nlevels) = 0.
         qs1d(1:Npoints,1:Nlevels) = 0.
         qg1d(1:Npoints,1:Nlevels) = 0.
         qh1d(1:Npoints,1:Nlevels) = 0.
         rainr1d(1:Npoints,1:Nlevels) = 0.
         hailr1d(1:Npoints,1:Nlevels) = 0.
         snowr1d(1:Npoints,1:Nlevels) = 0.
         graupelr1d(1:Npoints,1:Nlevels) = 0.

         ! INITIALISE REFF
         Reff(1:Npoints, 1:Nlevels, 1:N_Hydro) = 0.

         ! CALCULATE COLUMN PROFILES THAT AREN'T NATIVE TO WRF:
         do k = 1, Nlevels     ! k loop (vertical)
             ! full pressure field
             p1d(1:Npoints,k) = pb(i,k,j) + p(i,k,j)
             ! geopotential height
             ght(1:Npoints,k) = ( phb(i,k,j) + ph(i,k,j) ) / G
             ! specific humidity
             sh1d(1:Npoints,k) = qv(i,k,j) / (1 + qv(i,k,j) )
             ! relative humidity
             rh1d(1:Npoints,k) = ( qv(i,k,j) / ( (pq0 / psfc(i,j)) * exp( a2*( t(i,k,j) - a3) / ( t(i,k,j) - a4) ) ) )*100.
             ! check: if relative humidity is greater than 1, set it to 1!
             IF ( rh1d(1,k) .GT. 1. ) rh1d(1,k) = 1. ! BEWARE: ONLY WORKS IF NPOINTS = 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT!
             ! cloud optical thickness
             taucld1d(1:Npoints,k) = taucldc(i,k,j) + taucldi(i,k,j)
              ! cloud emissivity = 1-exp(-cloud optical depth)
             semiss(1:Npoints,k) = 1.0 - exp(-taucld1d(1:Npoints,k)) ! BEWARE: THIS FAILS HORRIBLY IF NPOINTS NE 1!!! WE HAVE SET NPOINTS = 1, BUT WATCH OUT!
         end do

         ! LINEAR INTERPOLATION TO CALCULATE FIELDS ON HALF LEVELS:
         do k = 1, Nlevels     ! k loop (vertical)
            if ( k .lt. Nlevels ) then
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
         do k = 1, Nlevels   
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
         ! NOTE: The following code only works for Npoints = 1, since Reff(Npoints, Nlevels, N_Hydro)
         ! TODO: Do we need to generalise this for the case where Npoints > 1?  
         !       Probably not if we are running COSP on each column...?
         Reff(1,1:Nlevels,I_LSCLIQ) = eradc1d*1.e-6  ! LSCLIQ
         Reff(1,1:Nlevels,I_LSCICE) = eradi1d*1.e-6  ! LSCICE
         Reff(1,1:Nlevels,I_LSRAIN) = eradr1d*1.e-6  ! LSRAIN
         Reff(1,1:Nlevels,I_LSSNOW) = erads1d*1.e-6  ! LSSNOW
         Reff(1,1:Nlevels,I_CVCLIQ) = eradc1d*1.e-6  ! CVCLIQ - USE SAME VALUES AS FOR LSCLIQ
         Reff(1,1:Nlevels,I_CVCICE) = eradi1d*1.e-6  ! CVCICE - USE SAME VALUES AS FOR LSCICE
         Reff(1,1:Nlevels,I_CVRAIN) = eradr1d*1.e-6  ! CVRAIN -  USE SAME VALUES AS FOR LSRAIN
         Reff(1,1:Nlevels,I_CVSNOW) = erads1d*1.e-6  ! CVSNOW -  USE SAME VALUES AS FOR LSSNOW
         Reff(1,1:Nlevels,I_LSGRPL) = eradg1d*1.e-6  ! LSGRPL     

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
        gbx%v_wind  = v10_ij ! surface v-wind
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
        call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,N_hydro,sghydro)
        call construct_cosp_radarstats(cfg,Npoints,Ncolumns,Nlevels,N_hydro,stradar)
        call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,Nlevels,N_hydro,PARASOL_NREFL,stlidar)

        !-----------------------------------------------------------------------
        ! CALL SIMULATOR
        !-----------------------------------------------------------------------

           call cosp_simulator(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)


        ! TRANSFER COSP OUTPUT TO WRF VARIABLES (following cosp_io.F90)
        ! NEED TO BE CAREFUL WITH DIMENSIONS HERE! Lots of arrays are multi-dimensional so for writing COSP output
        ! we have to make sure we have added the correct dimensions to registry.dimspec, then specified the size of 
        ! the output array in the relevant Registry file, as well as taking care of the dimensions locally and in the 
        ! call to this code from the solver.

        ! 2D VARIABLES (i,j)
        COSP_CLLCALIPSO(i,j) = stlidar%cldlayer(1,1) ! CALIPSO low level cloud fraction (profile)
        COSP_CLMCALIPSO(i,j) = stlidar%cldlayer(1,2) ! CALIPSO mid level cloud fraction (profile)
        COSP_CLHCALIPSO(i,j) = stlidar%cldlayer(1,3) ! CALIPSO high level cloud fraction (profile)
        COSP_CLTCALIPSO(i,j) = stlidar%cldlayer(1,4) ! CALIPSO total cloud fraction (profile)
        COSP_CLTLIDARRADAR(i,j) = stradar%radar_lidar_tcc(1) ! Lidar and radar total cloud fraction
        COSP_CLTISCCP(i,j) = isccp%totalcldarea(1) ! ISCCP total cloud fraction
        COSP_PCTISCCP(i,j) = isccp%meanptop(1) ! ISCCP mean cloud top pressure
        COSP_TAUISCCP(i,j) = isccp%meantaucld(1) ! ISCCP mean optical depth
        COSP_ALBISCCP(i,j) = isccp%meanalbedocld(1) ! ISCCP mean cloud albedo
        COSP_MEANTBCLRISCCP(i,j) = isccp%meantbclr(1) ! ISCCP mean clear-sky brightness temp
        COSP_MEANTBISCCP(i,j) = isccp%meantb(1) ! ISCCP mean all-sky brightness temp
        COSP_CLTMODIS(i,j) = modis%Cloud_Fraction_Total_Mean(1) ! MODIS total cloud fraction
        COSP_CLWMODIS(i,j) = modis%Cloud_Fraction_Water_Mean(1) ! MODIS liquid cloud fraction
        COSP_CLIMODIS(i,j) = modis%Cloud_Fraction_Ice_Mean(1) ! MODIS ice cloud fraction
        COSP_CLHMODIS(i,j) = modis%Cloud_Fraction_High_Mean(1) ! MODIS high level cloud fraction
        COSP_CLMMODIS(i,j) = modis%Cloud_Fraction_Mid_Mean(1) ! MODIS mid level loud fraction
        COSP_CLLMODIS(i,j) = modis%Cloud_Fraction_Low_Mean(1) ! MODIS low level cloud fraction
        COSP_TAUTMODIS(i,j) = modis%Optical_Thickness_Total_Mean(1) ! MODIS total cloud optical thickness
        COSP_TAUWMODIS(i,j) = modis%Optical_Thickness_Water_Mean(1) ! MODIS liquid cloud optical thickness
        COSP_TAUIMODIS(i,j) = modis%Optical_Thickness_Ice_Mean(1) ! MODIS ice cloud optical thickness
        COSP_TAUTLOGMODIS(i,j) = modis%Optical_Thickness_Total_LogMean(1) ! MODIS total cloud optical thickness (log10 mean)
        COSP_TAUWLOGMODIS(i,j) = modis%Optical_Thickness_Water_LogMean(1) ! MODIS liquid cloud optical thickness (log10 mean)
        COSP_TAUILOGMODIS(i,j) = modis%Optical_Thickness_Ice_LogMean(1) ! MODIS ice cloud optical thickness (log10 mean)
        COSP_REFFCLWMODIS(i,j) = modis%Cloud_Particle_Size_Water_Mean(1) ! MODIS liquid cloud particle size
        COSP_REFFCLIMODIS(i,j) = modis%Cloud_Particle_Size_Ice_Mean(1) ! MODIS ice cloud particle size
        COSP_PCTMODIS(i,j) = modis%Cloud_Top_Pressure_Total_Mean(1) ! MODIS cloud top pressure
        COSP_LWPMODIS(i,j) = modis%Liquid_Water_Path_Mean(1) ! MODIS cloud liquid water path
        COSP_IWPMODIS(i,j) = modis%Ice_Water_Path_Mean(1) ! MODIS cloud ice water path
        COSP_TOFFSET(i,j) = gbx%toffset(1) ! time difference
        ! These are 2D variables because we use Ncolumns = 1 for WRF implementation:
        COSP_BOXTAUISCCP(i,j) = isccp%boxtau(1,1) ! Optical depth in each column as calculated by ISCCP simulator
        COSP_BOXTOPISCCP(i,j) = isccp%boxptop(1,1) ! Cloud top pressure in each column as calculated by ISCCP simulator

        ! 3D variables: some are on stat grid i.e. (i,l,j)
        ! registry.dimspec: Nlr = 40 (= number of levels in statistical output)
        COSP_CLCALIPSO(i,1:Nlr,j) = stlidar%lidarcld(1,:) ! CALIPSO cloud area fraction
        COSP_CLCALIPSO2(i,1:Nlr,j) = stradar%lidar_only_freq_cloud(1,:) ! CALIPSO cloud fraction undetected by CloudSat
        ! LIDARBETAMOL532 is on model grid i.e. (i,k,j)
        COSP_LIDARBETAMOL532(i,1:Nlevels,j) = sglidar%beta_mol(1,:) ! Lidar molecular backscatter (532 nm)
        
        ! dims of PARASOLREFL(i,j,p)
        COSP_PARASOLREFL(i,j,1:Nrefl) = stlidar%parasolrefl(1,:) ! PARASOL reflectance         

        ! 4D variables are actually 3D: (Npoints, Ncolumns, Nlevels) = (i,k,j) because we have Ncolumns = 1 for WRF implementation...
        ! The following arrays are on model levels (1:Nlevels) not the statistical grid (Nlr = 40)
        COSP_DBZE94(i,1:Nlevels,j) = sgradar%Ze_tot(1,1,:) ! CloudSat radar reflectivity
        COSP_ATB532(i,1:Nlevels,j) = sglidar%beta_tot(1,1,:) ! Lidar attenuated total backscatter (532 nm)
        COSP_FRACOUT(i,1:Nlevels,j) = sgx%frac_out(1,1,:) ! Subcolumn output from SCOPS

        ! The following have non-spatial 4th (and some 3rd) dimensions:
        ! CFADDBZE94(i,l,j,b); DBZE_BINS is read in from cosp_constants
        COSP_CFADDBZE94(i,1:Nlr,j,1:DBZE_BINS) = stradar%cfad_ze(1,:,:)
        ! CFADLIDARSR532: (i,l,j,s); SR_BINS is read in from cosp_constants
        COSP_CFADLIDARSR532(i,1:Nlr,j,1:SR_BINS) = stlidar%cfad_sr(1,:,:) ! CALIPSO scattering ratio CFAD
        ! CLISCCP and CLMODIS: CLISCCP(i,j,d,e); CLMODIS(i,j,da,db)
        COSP_CLISCCP(i,j,1:7,1:7) = isccp%fq_isccp(1,:,:) ! ISCCP cloud area fraction
        COSP_CLMODIS(i,j,1:7,1:7) = modis%Optical_Thickness_vs_Cloud_Top_Pressure(1,:,:) ! MODIS cloud area fraction"
        !CLMISR(i,j,d,m)
        COSP_CLMISR(i,j,1:7,1:MISR_N_CTH) = misr%fq_MISR(1,:,:) ! MISR cloud fraction


      end if ! only run COSP if gridbox is sunlit TODO: THIS ONLY APPLIES FOR SOME OF THE SIMULATORS! SEE CAM IMPLEMENTATION

        ! DEALLOCATE MEMORY IN DERIVED TYPES: free_cosp_X is accessed from cosp_types, this is within the i,j loop
        call free_cosp_gridbox(gbx)
        call free_cosp_vgrid(vgrid)
        call free_cosp_subgrid(sgx)
        call free_cosp_sgradar(sgradar)
        call free_cosp_radarstats(stradar)
        call free_cosp_sglidar(sglidar)
        call free_cosp_lidarstats(stlidar)
        call free_cosp_isccp(isccp)
        call free_cosp_misr(misr)
        call free_cosp_modis(modis)


      end do ! i-point loop
      end do ! j-point loop

      ! TESTING WRITING TO NEW VARIABLES!!!
      ! eradc(ims:ime, kms:kme, jms:jme)=-15
      ! eradi(ims:ime, kms:kme, jms:jme)=-15
      ! erads(ims:ime, kms:kme, jms:jme)=-15
      ! eradr(ims:ime, kms:kme, jms:jme)=-15
      ! eradg(ims:ime, kms:kme, jms:jme)=-15
      ! cldfra(ims:ime, kms:kme, jms:jme)=-15

     END DO ! tiles loop
   end subroutine wrf_cosp_simulator
 END MODULE module_cosp_simulator_test
!END MODULE module_cosp_wrf
