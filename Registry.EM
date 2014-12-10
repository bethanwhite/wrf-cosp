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
state    real  COSP_TOFFSET         ij       misc        1         -      h        "TOFFSET"              "time difference in days between time of profile and recorded time"     ""
# ISCCP simulator
state    real  COSP_ALBISCCP        ij      misc        1         -      h        "ALBISCCP"             "ISCCP mean cloud albedo"  					      ""
state    real  COSP_BOXTOPISCCP     ij      misc        1         -      h        "BOXTOPISCCP"          "Cloud top pressure in each column as calculated by ISCCP simulator"    "Pa"
state    real  COSP_BOXTAUISCCP     ij      misc        1         -      h        "BOXTAUISCCP"          "Optical depth in each column as calculated by ISCCP simulator"         ""
state    real  COSP_CLISCCP         ijde    misc        1         -      h        "CLISCCP"	            "ISCCP cloud area fraction"                                             ""
state    real  COSP_CLTISCCP        ij      misc        1         -      h        "CLTISCCP"             "ISCCP total cloud fraction"                                            ""
state    real  COSP_MEANTBCLRISCCP  ij      misc        1         -      h        "MEANTBCLRISCCP"       "Mean clear-sky 10.5 micron brightness temp from the ISCCP simulator"   "K"
state    real  COSP_MEANTBISCCP     ij      misc        1         -      h        "MEANTBISCCP"          "Mean all-sky 0.5 micron brightness temp from the ISCCP simulator"      "K"
state    real  COSP_TAUISCCP        ij      misc        1         -      h        "TAUISCCP"             "Mean optical depth as calculated by the ISCCP simulator"               ""
state    real  COSP_PCTISCCP        ij      misc        1         -      h        "PCTISCCP"             "ISCCP mean cloud top pressure"					      "Pa"
# MISR simulator
state    real  COSP_CLMISR          ijdm    misc        1         -      h        "CLMISR"               "MISR cloud fraction"                                                   ""
# CALIPSO simulator
state    real  COSP_ATB532          ikj     misc        1         -      h       "LATB532"              "Lidar attenuated total backscatter (532 nm)"                           ""
state    real  COSP_CFADLIDARSR532  iljs    misc        1         -      h        "CFADLIDARSR532"       "CALIPSO scattering ratio CFAD"                                         ""
state    real  COSP_CLCALIPSO       ilj     misc        1         -      h        "CLCALIPSO"            "CALIPSO cloud area fraction" 					      ""
state    real  COSP_CLHCALIPSO	    ij      misc        1         -      h       "CLHCALIPSO"           "CALIPSO high level cloud fraction"                                     ""
state    real  COSP_CLLCALIPSO      ij      misc        1         -      h        "CLLCALIPSO"           "CALIPSO low level cloud fraction"                                      ""
state    real  COSP_CLMCALIPSO      ij      misc        1         -      h        "CLMCALIPSO"           "CALIPSO mid level cloud fraction" 	                              ""
state    real  COSP_CLTCALIPSO      ij      misc        1         -      h        "CLTCALIPSO"           "CALIPSO total cloud fraction"                                          ""
# dims of PARASOLREFL(i,j,Nrefl)?
state    real  COSP_PARASOLREFL     ijp     misc        1         -      h        "PARASOLREFL"          "PARASOL reflectance"                                                   ""
# CloudSat simulator  
state	 real  COSP_CFADDBZE94      iljb    misc   	    1         -      h	     "CFADDBZE94"           "CloudSat radar reflectivity CFAD"				      "dBZ"
state    real  COSP_DBZE94          ikj		misc	    1 	      -	     h       "DBZE9"                "CloudSat radar reflectivity"					      "dBZ"
# Use lidar and radar
state    real  COSP_CLTLIDARRADAR   ik      misc        1         -      h        "CLTLIDARRADAR"        "Lidar and radar total cloud fraction"                                  ""
state    real  COSP_CLCALIPSO2      ilj     misc        1         -      h        "CLCALIPSO2"           "CALIPSO cloud fraction undetected by CloudSat"                         ""
# MODIS simulator
state    real  COSP_CLHMODIS        ij      misc        1         -      h        "CLHMODIS"             "MODIS high level cloud fraction"                                       ""
state    real  COSP_CLIMODIS        ij      misc        1         -      h        "CLIMODIS"             "MODIS ice cloud fraction"                                              ""
state    real  COSP_CLLMODIS        ij      misc        1         -      h        "CLLMODIS"             "MODIS low level cloud fraction"                                        ""
state    real  COSP_CLMMODIS        ij      misc        1         -      h        "CLMMODIS"             "MODIS mid level loud fraction"                                         ""
# required? dims?
state    real  COSP_CLMODIS         ijde    misc        1         -      h        "CLMODIS"              "MODIS cloud area fraction"                                             ""
state    real  COSP_CLTMODIS        ij      misc        1         -      h        "CLTMODIS"             "MODIS total cloud fraction"                                            ""
state    real  COSP_CLWMODIS        ij      misc        1         -      h        "CLWMODIS"             "MODIS liquid cloud fraction"                                           ""
state    real  COSP_IWPMODIS        ij      misc        1         -      h        "IWPMODIS"             "MODIS cloud ice water path"                                            ""
state    real  COSP_LWPMODIS        ij      misc        1         -      h        "LWPMODIS"             "MODIS cloud liquid water path"                                         ""
state    real  COSP_PCTMODIS        ij      misc        1         -      h        "PCTMODIS"             "MODIS cloud top pressure"                                              ""
state    real  COSP_REFFCLIMODIS    ij      misc        1         -      h        "REFFCLIMODIS"         "MODIS ice cloud particle size"                                         "um"
state    real  COSP_REFFCLWMODIS    ij      misc        1         -      h        "REFFCLWMODIS"         "MODIS liquid cloud particle size"                                      "um"
state    real  COSP_TAUILOGMODIS    ij      misc        1         -      h        "TAUILOGMODIS"         "MODIS ice cloud optical thickness (log10 mean)"                        ""
state    real  COSP_TAUIMODIS       ij      misc        1         -      h        "TAUIMODIS"            "MODIS ice cloud optical thickness"                                     ""
state    real  COSP_TAUTLOGMODIS    ij      misc        1         -      h        "TAUTLOGMODIS"         "MODIS total cloud optical thickness (log10 mean)"                      ""
state    real  COSP_TAUTMODIS       ij      misc        1         -      h        "TAUTMODIS"            "MODIS total cloud optical thickness"                                   ""
state    real  COSP_TAUWLOGMODIS    ij      misc        1         -      h        "TAUWLOGMODIS"         "MODIS liquid cloud optical thickness (log10 mean)"                     ""
state    real  COSP_TAUWMODIS       ij      misc        1         -      h        "TAUWMODIS"            "MODIS liquid cloud optical thickness"                                  ""
# RRTOV simulator
state    real  COSP_TBRTTOV         ijn     misc        1         -      h        "TBRTTOV"              "Mean clear-sky brightness temp as calculated by RRTOV"                 "K"
# Other COSP outputs
state    real  COSP_FRACOUT         ikj     misc        1         -      h        "FRACOUT"              "Subcolumn output from SCOPS"                                           ""
state    real  COSP_LIDARBETAMOL532 il      misc        1         -      h        "LIDARBETAMOL532"      "Lidar molecular backscatter (532 nm)"                                  ""
#state    real  COSP_SUNLIT      ij      misc        1         -      h       "COSP_SUNLIT"           "COSP sunlit fraction"                                                  ""
#state    real  COSP_SUNLIT_AV   ij      misc        1         -      h       "COSP_SUNLIT_AV"        "COSP average sunlit fraction"                                          ""

 
