import numpy as np


def wrap_longitude(grid_long):
    """
    Makes sure that, if the longitude is near the boundary (0 or 360) that we wrap the values beyond
    360 so it appears nicely on a map
    :param grid_long: array of longitude values
    :return: wrray of longitude values that can extend past 360
    """
    
    neg_long = np.argwhere(grid_long < 0)
    grid_long[neg_long] = grid_long[neg_long] + 360

    # if we have data close to upper boundary (360), then wrap some of the data round
    # so it appears on the map
    top_long = np.argwhere(grid_long >= 320)
    if top_long.__len__() != 0:
        bottom_long = np.argwhere(grid_long <= 40)
        grid_long[bottom_long] = 360 + grid_long[bottom_long]

    return grid_long
