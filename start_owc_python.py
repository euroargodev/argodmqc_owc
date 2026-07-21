"""starting code"""

import multiprocessing
import sys
import time
import warnings
from functools import partial

import click

import pyowc as owc

warnings.filterwarnings("ignore", category=RuntimeWarning)


@click.command
@click.option("--headless", is_flag=True)
@click.option("--floats", required=True)
@click.option("-appendRef", "plot_file_suffix")
def main(headless, plot_file_suffix, floats) -> None:
    """Entry point for processing."""
    if plot_file_suffix and not plot_file_suffix.startswith("_"):
        raise Exception("appendRef must start with an underscore")

    FLOAT_NAMES = floats.split(",")
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
        start = time.time()
        owc.configuration.set_calseries("/", flt, USER_CONFIG)
        owc.calibration.calc_piecewisefit("/", flt, USER_CONFIG)
        owc.dashboard("/", flt, USER_CONFIG, headless=headless, file_suffix=plot_file_suffix)
        mid = time.time()
        print("Time for float: ", mid - start)

if __name__ == "__main__":
    main()
