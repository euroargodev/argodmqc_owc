import numpy as np

def covarxy_pv(input_coords, coords, long, lat, phi, use_pv):
    """
    Returns a matrix for the horizontal covariance
    :param input_coords: the input coordinates of the the float profile
    :param coords: coordinates for all the float profiles
    :param long: longitude scale
    :param lat: latitude scale
    :param phi: potential gradient
    :param use_pv: whether or not to use potential vorticity
    :return: horizontal covariance matrix
    """

    # Get data dimensions
    input_coord_rows = input_coords.shape[0]
    coord_rows = coords.shape[0]

    # Set up matrices
    cov = np.zeros((input_coord_rows, coord_rows))
    one = np.ones((coord_rows, 1))

    # Derive the planetary vorticity at each point

    # Get the depth for each data point
    z_input_coords = input_coords[2]
    z_coords = coords[:, 2]

    print(input_coords[0])
    print(np.sin(input_coords[0]))
    print(z_input_coords)

    pv_input_coords = (2 * 7.292 * 10 ** -5 * np.sin(input_coords[0] * np.pi/180))/z_input_coords
    print(pv_input_coords)
