"""A script to allow profiling of the update_salinity_mapping function.

Run this code as follows from the repository root directory:

python -m cProfile -o baseline.stats utils/profile_update_salinity_mapping.py

This will export the results to the ``baseline.stats`` file, which can be read using:

python -m pstats baseline.stats
"""
import tempfile

import pyowc

if __name__ == "__main__":

    FLOAT_NAME = "3901960"
    USER_CONFIG = pyowc.configuration.load()

    print(pyowc.configuration.print_cfg(USER_CONFIG))

    with tempfile.TemporaryDirectory() as tmp_dir:
        # override the output directory for this call
        tmp_config = {**USER_CONFIG, "FLOAT_MAPPED_DIRECTORY": tmp_dir}

        pyowc.calibration.update_salinity_mapping("/", tmp_config, FLOAT_NAME)
