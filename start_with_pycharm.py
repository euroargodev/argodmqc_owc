""" starting code
"""
import multiprocessing
from functools import partial
import time
import warnings
import pyowc as owc

warnings.filterwarnings("ignore", category=RuntimeWarning)

if __name__ == '__main__':

    FLOAT_NAMES = ["3901960"]  # add float names here
    USER_CONFIG = owc.configuration.load()  # fetch the default configuration and parameters
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
    for flt in FLOAT_NAMES:
        owc.configuration.set_calseries("/", flt, USER_CONFIG)
        owc.calibration.calc_piecewisefit("/", flt, USER_CONFIG)
        owc.dashboard("/", flt, USER_CONFIG)
        mid = time.time()
        print("Time for float: ", mid - start)

