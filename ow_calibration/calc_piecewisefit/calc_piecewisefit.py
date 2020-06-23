"""
-----Calculate Piece wise fit-----

Written by: Annie Wong
When: xx/10/2008
Converted to python by: Edward Small
When: 16/06/2020

Calculate the fit of each break and calibrate salinities

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


def calc_piecewisefit(float_dir, float_name, system_config):
    """
    calibrate salinities
    :param float_dir: float directory name
    :param float_name: name of float
    :param system_config: configuration parameter set up
    :return: Nothing, save output
    """
    print(float_dir)
    print(float_name)
    print(system_config)