"""starting code"""

import multiprocessing
import sys
import time
import warnings
from functools import partial

import pyowc as owc

warnings.filterwarnings("ignore", category=RuntimeWarning)


def main() -> None:
    """Entry point for processing."""
    # Handle any input args
    args = sys.argv[1:]
    headless = "--headless" in args
    plot_file_suffix = ""
    try:
        append_index = args.index("-appendRef")
        if len(args) > append_index + 1:
            plot_file_suffix = args[append_index + 1]
            if not plot_file_suffix.startswith("_"):
                raise Exception("appendRef must start with an underscore")
    except ValueError:
        pass
    FLOAT_NAMES = ["3901960"]
    config_file_location = "owc_config.json"

    USER_CONFIG = owc.utilities.load_configuration_from_json_file(config_file_location)
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
        owc.dashboard("/", flt, USER_CONFIG, headless=headless, file_suffix=plot_file_suffix)
        mid = time.time()
        print("Time for float: ", mid - start)

if __name__ == "__main__":
    main()
