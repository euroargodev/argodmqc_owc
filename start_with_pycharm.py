from argodmqc_owc import pyowc as owc
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import multiprocessing
from functools import partial

FLOAT_NAMES = ["3901960", "FLOAT NAME"] # add float names here
USER_CONFIG = owc.configuration.load() # fetch the default configuration and parameters
USER_CONFIG['MAP_USE_PV'] = 0 # Possibly tune options
print(owc.configuration.print_cfg(USER_CONFIG))

#owc.configuration.set_calseries("/", FLOAT_NAME, USER_CONFIG)

pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
func = partial(owc.calibration.update_salinity_mapping, "/", USER_CONFIG)
pool.map(func, FLOAT_NAMES)
pool.close()
pool.join()

#owc.calibration.calc_piecewisefit("/", FLOAT_NAME, USER_CONFIG)
#owc.dashboard("/", FLOAT_NAME, USER_CONFIG)
