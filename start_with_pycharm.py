from argodmqc_owc import pyowc as owc
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

FLOAT_NAMES = ["3901960", "OTHER FLOAT"] # add float numbers here
USER_CONFIG = owc.configuration.load() # fetch the default configuration and parameters
USER_CONFIG['MAP_USE_PV'] = 0 # Possibly tune options
print(owc.configuration.print_cfg(USER_CONFIG))

#owc.configuration.set_calseries("/", FLOAT_NAME, USER_CONFIG)

# only need salinity mapping

for float in FLOAT_NAMES:
    owc.calibration.update_salinity_mapping("/", float, USER_CONFIG)


#owc.calibration.calc_piecewisefit("/", FLOAT_NAME, USER_CONFIG)
#owc.dashboard("/", FLOAT_NAME, USER_CONFIG)
