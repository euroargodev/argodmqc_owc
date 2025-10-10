"""starting code"""

import multiprocessing
import time
import warnings
from functools import partial

import pyowc as owc

warnings.filterwarnings("ignore", category=RuntimeWarning)


def main() -> None:
    """Entry point for processing."""
    FLOAT_NAMES = ["3901960"]
    config_file_location = "owc_config.json"

    USER_CONFIG = owc.utilities.load_configuration_from_json_file(config_file_location)
    print("App configuration is as follows: \n")
    for key, value in USER_CONFIG.items():
        print(key, value)

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

if __name__ == "__main__":
    main()
