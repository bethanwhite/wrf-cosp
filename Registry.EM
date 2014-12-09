#  Registry file, EM

#INCLUDES
include registry.dimspec
include Registry.EM_COMMON
include registry.io_boilerplate

include registry.fire
include registry.avgflx
include registry.stoch
include registry.les
include registry.cam
include registry.clm
include registry.ssib
include registry.diags
include registry.bdy_perturb

# added to output 5 for ESMF
state   real   landmask            ij    misc          1     -     i0125rh05d=(interp_fcnm_imask)u=(copy_fcnm)   "LANDMASK"      "LAND MASK (1 FOR LAND, 0 FOR WATER)"  ""

# Masked SST interpolation from the CG
#state    real   SST              ij    misc        1         -     i01245rh05d=(interp_mask_field:lu_index,iswater)f=(p2c_mask:lu_index,tslb,num_soil_layers,iswater)   "SST"              "SEA SURFACE TEMPERATURE" "K"
# Simple SST interpolation from the CG
#state    real   SST              ij    misc        1         -     i01245rh05d=(interp_mask_field:lu_index,iswater)f=(p2c)   "SST"              "SEA SURFACE TEMPERATURE" "K"
state    real   SST              ij     misc        1         -     i01245rh05d=(interp_mask_field:lu_index,iswater)   "SST"              "SEA SURFACE TEMPERATURE" "K"


# Registry entries that are exclusive to Registry.EM

# Chem Scalars
rconfig   integer     chem_opt            namelist,physics	max_domains    0       rh       "chem_opt"              ""      ""
state   real    -              ikjftb   chem        1         -     -    -                                       
#
# Tracer Scalars: If you want full resolved and non-resolved dispersion, compile WRF-Chem
#
state   real    -          ikjftb  tracer        1         -     -    -                                       
state   real    tr17_1     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_1"         "tr17_1"          -
state   real    tr17_2     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_2"         "tr17_2"          -
state   real    tr17_3     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_3"         "tr17_3"          -
state   real    tr17_4     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_4"         "tr17_4"          -
state   real    tr17_5     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_5"         "tr17_5"          -
state   real    tr17_6     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_6"         "tr17_6"          -
state   real    tr17_7     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_7"         "tr17_7"          -
state   real    tr17_8     ikjftb  tracer        1         -     irhusdf=(bdy_interp:dt)    "tr17_8"         "tr17_8"

package   tracer_test1  tracer_opt==2       -             tracer:tr17_1,tr17_2,tr17_3,tr17_4,tr17_5,tr17_6,tr17_7,tr17_8


# COSP input variables 
state    real   OZONEMR	        ikj     misc        1       -       rh		"OZONEMR"			    "COSP OZONE MIXING RATIO" "" 	
state	 real	ERADC		    ikj		misc		1		-		rh		"ERADC"				     "CLOUD DROPLET EFFECTIVE RADIUS" ""
state    real   ERADI           ikj     misc        1       -       rh       "ERADI"                  "ICE EFFECTIVE RADIUS" ""
state    real   ERADS           ikj     misc        1       -       rh       "ERADS"                  "SNOW EFFECTIVE RADIUS" ""
state    real   ERADR           ikj     misc        1       -       rh       "ERADR"                  "RAIN DROPLET EFFECTIVE RADIUS" ""
state    real   ERADG           ikj     misc        1       -       rh       "ERADG"                  "GRAUPEL EFFECTIVE RADIUS" ""
state    real   PRECR           ikj     misc        1       -       r       "PRECR"                  "LARGE-SCALE RAIN SEDIMENTATION FLUX" ""
state    real   PRECI           ikj     misc        1       -       r       "PRECI"                  "LARGE-SCALE ICE SEDIMENTATION FLUX" ""
state    real   PRECS           ikj     misc        1       -       r       "PRECS"                  "LARGE-SCALE SNOW SEDIMENTATION FLUX" "" 
state    real   PRECG           ikj     misc        1       -       r       "PRECG"                  "LARGE-SCALE GRAUPEL SEDIMENTATION FLUX" ""

# COSP output variables
# if variables only in 1D (curtain mode)
state    real  TOFFSET         ij       misc        1         -      h        "TOFFSET"              "time difference in days between time of profile and recorded time"     ""
# ISCCP simulator
state    real  ALBISCCP        ik      misc        1         -      h        "ALBISCCP"             "ISCCP mean cloud albedo"  					      ""
state    real  BOXTOPISCCP     ikj     misc        1         -      h        "BOXTOPISCCP"          "Cloud top pressure in each column as calculated by ISCCP simulator"    "Pa"
state    real  BOXTAUISCCP     ikj     misc        1         -      h        "BOXTAUISCCP"          "Optical depth in each column as calculated by ISCCP simulator"         ""
# check CLISCCP dims: should be (i,j,7,7)?
state    real  CLISCCP         ij      misc        1         -      h        "CLISCCP"	            "ISCCP cloud area fraction"                                             ""
state    real  CLTISCCP        ij      misc        1         -      h        "CLTISCCP"             "ISCCP total cloud fraction"                                            ""
state    real  MEANTBCLRISCCP  ij      misc        1         -      h        "MEANTBCLRISCCP"       "Mean clear-sky 10.5 micron brightness temp from the ISCCP simulator"   "K"
state    real  MEANTBISCCP     ij      misc        1         -      h        "MEANTBISCCP"          "Mean all-sky 0.5 micron brightness temp from the ISCCP simulator"      "K"
state    real  TAUISCCP        ij      misc        1         -      h        "TAUISCCP"             "Mean optical depth as calculated by the ISCCP simulator"               ""
state    real  PCTISCCP        ij      misc        1         -      h        "PCTISCCP"             "ISCCP mean cloud top pressure"					      "Pa"
# MISR simulator
state    real  CLMISR          ik      misc        1         -      h        "CLMISR"               "MISR cloud fraction"                                                   ""
# CALIPSO simulator
state    real  LATB532         ik      misc        1         -      h       "LATB532"              "Lidar attenuated total backscatter (532 nm)"                           ""
# dims of (i,k,j,SR_bins)?
state    real  CFADLIDARSR532  ik      misc        1         -      h        "CFADLIDARSR532"       "CALIPSO scattering ratio CFAD"                                         ""
state    real  CLCALIPSO       ikj     misc        1         -      h        "CLCALIPSO"            "CALIPSO cloud area fraction" 					      ""
state    real  CLHCALIPSO	   ij      misc        1         -      h       "CLHCALIPSO"           "CALIPSO high level cloud fraction"                                     ""
state    real  CLLCALIPSO      ij      misc        1         -      h        "CLLCALIPSO"           "CALIPSO low level cloud fraction"                                      ""
state    real  CLMCALIPSO      ij      misc        1         -      h        "CLMCALIPSO"           "CALIPSO mid level cloud fraction" 	                              ""
state    real  CLTCALIPSO      ij      misc        1         -      h        "CLTCALIPSO"           "CALIPSO total cloud fraction"                                          ""
# dims of PARASOLREFL(i,j,Nrefl)?
state    real  PARASOLREFL     ij      misc        1         -      h        "PARASOLREFL"          "PARASOL reflectance"                                                   ""
# CloudSat simulator  
# dims of (i,k,j,DBZE_BINS)?
state	 real  CFADDBZE94      ikj      misc   	    1         -      h	     "CFADDBZE94"           "CloudSat radar reflectivity CFAD"				      "dBZ"
state    real  DBZE94          ikj		misc	    1 	      -	     h       "DBZE9"                "CloudSat radar reflectivity"					      "dBZ"
# Use lidar and radar
state    real  CLTLIDARRADAR   ik      misc        1         -      h        "CLTLIDARRADAR"        "Lidar and radar total cloud fraction"                                  ""
state    real  CLCALIPSO2      ikj      misc        1         -      h        "CLCALIPSO2"           "CALIPSO cloud fraction undetected by CloudSat"                         ""
# MODIS simulator
state    real  CLHMODIS        ij      misc        1         -      h        "CLHMODIS"             "MODIS high level cloud fraction"                                       ""
state    real  CLIMODIS        ij      misc        1         -      h        "CLIMODIS"             "MODIS ice cloud fraction"                                              ""
state    real  CLLMODIS        ij      misc        1         -      h        "CLLMODIS"             "MODIS low level cloud fraction"                                        ""
state    real  CLMMODIS        ij      misc        1         -      h        "CLMMODIS"             "MODIS mid level loud fraction"                                         ""
# required? dims?
state    real  CLMODIS         ik      misc        1         -      h        "CLMODIS"              "MODIS cloud area fraction"                                             ""
state    real  CLTMODIS        ij      misc        1         -      h        "CLTMODIS"             "MODIS total cloud fraction"                                            ""
state    real  CLWMODIS        ij      misc        1         -      h        "CLWMODIS"             "MODIS liquid cloud fraction"                                           ""
state    real  IWPMODIS        ij      misc        1         -      h        "IWPMODIS"             "MODIS cloud ice water path"                                            ""
state    real  LWPMODIS        ij      misc        1         -      h        "LWPMODIS"             "MODIS cloud liquid water path"                                         ""
state    real  PCTMODIS        ij      misc        1         -      h        "PCTMODIS"             "MODIS cloud top pressure"                                              ""
state    real  REFFCLIMODIS    ij      misc        1         -      h        "REFFCLIMODIS"         "MODIS ice cloud particle size"                                         "um"
state    real  REFFCLWMODIS    ij      misc        1         -      h        "REFFCLWMODIS"         "MODIS liquid cloud particle size"                                      "um"
state    real  TAUILOGMODIS    ij      misc        1         -      h        "TAUILOGMODIS"         "MODIS ice cloud optical thickness (log10 mean)"                        ""
state    real  TAUIMODIS       ij      misc        1         -      h        "TAUIMODIS"            "MODIS ice cloud optical thickness"                                     ""
state    real  TAUTLOGMODIS    ij      misc        1         -      h        "TAUTLOGMODIS"         "MODIS total cloud optical thickness (log10 mean)"                      ""
state    real  TAUTMODIS       ij      misc        1         -      h        "TAUTMODIS"            "MODIS total cloud optical thickness"                                   ""
state    real  TAUWLOGMODIS    ij      misc        1         -      h        "TAUWLOGMODIS"         "MODIS liquid cloud optical thickness (log10 mean)"                     ""
state    real  TAUWMODIS       ij      misc        1         -      h        "TAUWMODIS"            "MODIS liquid cloud optical thickness"                                  ""
# RRTOV simulator
state    real  TBRTTOV         ikj     misc        1         -      h        "TBRTTOV"              "Mean clear-sky brightness temp as calculated by RRTOV"                 "K"
# Other COSP outputs
state    real  FRACOUT         ik      misc        1         -      h        "FRACOUT"              "Subcolumn output from SCOPS"                                           ""
state    real  LIDARBETAMOL532 ik      misc        1         -      h        "LIDARBETAMOL532"      "Lidar molecular backscatter (532 nm)"                                  ""
state    real  COSP_SUNLIT      ik      misc        1         -      h       "COSP_SUNLIT"           "COSP sunlit fraction"                                                  ""
state    real  COSP_SUNLIT_AV   ik      misc        1         -      h       "COSP_SUNLIT_AV"        "COSP average sunlit fraction"                                          ""

 
