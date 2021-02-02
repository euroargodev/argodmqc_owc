from argodmqc_owc import pyowc as owc
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
import multiprocessing
from functools import partial
import time

if __name__ == '__main__':

    FLOAT_NAMES = ["3901960"]  # add float names here
    USER_CONFIG = owc.configuration.load()  # fetch the default configuration and parameters
    #USER_CONFIG['MAP_USE_PV'] = 0  # Possibly tune options
    print(owc.configuration.print_cfg(USER_CONFIG))

    start = time.time()

    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    func = partial(owc.calibration.update_salinity_mapping, "/", USER_CONFIG)
    pool.map(func, FLOAT_NAMES)
    pool.close()
    pool.join()

    end = time.time()

    print("\nTOTAL TIME ELAPSED: ", end - start)

    # loop for sequential run
    for float in FLOAT_NAMES:
        owc.configuration.set_calseries("/", float, USER_CONFIG)
        owc.calibration.calc_piecewisefit("/", float, USER_CONFIG)
        owc.dashboard("/", float, USER_CONFIG)
        mid = time.time()
        print("Time for float: ", mid - start)

    # FLOAT_NAME = "3901960"
    # owc.configuration.set_calseries("/", FLOAT_NAME, USER_CONFIG)
    # owc.calibration.calc_piecewisefit("/", FLOAT_NAME, USER_CONFIG)
    # owc.dashboard("/", FLOAT_NAME, USER_CONFIG)
