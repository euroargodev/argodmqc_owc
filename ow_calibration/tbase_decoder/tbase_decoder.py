"""
-----tbase_decoder-----

Written by: Edward Small
When: 28/04/2020

The old matlab version of this uses an old .int file from NOAA which contains
5 arcminute data of global terrain. Whilst other more complete data sets now
exists, the data files for them are very large, so this file can be used for now,
and perhaps we can update to better data when we either move the code online, or
find a way to susbet it/compress it reasonably.

The .int file is very weird, and stores 16bit integers in binary. The below
code opens the file and converts the old binary back to the numbers we expect.
It then finds global terrain over a specified area before closing the file

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import struct
import numpy as np


# pylint: disable=too-many-locals
def get_topo_grid(min_long, max_long, min_lat, max_lat):
    """
    Find depth grid over given area using tbase.int file
    :param min_long: minimum longitudinal value for grid
    :param max_long: maximum longitudinal value for grid
    :param min_lat: minimum latidunal value for grid
    :param max_lat: maximum latidunal value for grid
    :return: matrices containing a uniform grid of latitudes
    and longitudes, along with the depth at these points
    """
    # manipulate input values to match file for decoding
    blat = int(np.max((np.floor(min_lat * 12), -90 * 12 + 1)))
    tlat = int(np.ceil(max_lat * 12))
    llong = int(np.floor(min_long * 12))
    rlong = int(np.ceil(max_long * 12))

    # use these values to form the grid
    lgs = np.arange(llong, rlong + 1, 1) / 12
    lts = np.flip(np.arange(blat, tlat + 1, 1) / 12, axis=0)

    if rlong > 360 * 12 - 1:
        rlong = rlong - 360 * 12
        llong = llong - 360 * 12

    if llong < 0:
        rlong = rlong + 360 * 12
        llong = llong + 360 * 12

    decoder = [llong, rlong, 90 * 12 - blat, 90 * 12 - tlat]

    # get the amount of elevation values we need
    nlat = int(round(decoder[2] - decoder[3])) + 1
    nlong = int(round(decoder[1] - decoder[0])) + 1

    # initialise matrix to hold z values
    topo = np.zeros((nlat, nlong))

    # Open the binary file
    elev_file = open("data/constants/tbase.int", "rb")

    # decode the file, and get values
    for i in range(nlat):
        elev_file.seek((i + decoder[3]) * 360 * 12 * 2 + decoder[0] * 2)

        for j in range(nlong):
            topo[i, j] = struct.unpack('h', elev_file.read(2))[0]

    # make the grid
    longs, lats = np.meshgrid(lgs, lts)

    # close the file
    elev_file.close()

    return topo, longs, lats
