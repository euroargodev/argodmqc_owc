"""
-----Map Data to Grid-----

Written by: Breck Owens
When: xx/11/2007
Converted to python by: Edward Small
When: 22/11/2019

An optimal mapping routine, taking data measured in arbitrary geographic locations and
mapping these data onto a more regular grid. As should happen in every mapping problem
the data are both mapped onto the prescribed grid, and the data locations (so that
the mapped field can be checked with the original data to ensure that the statistics
are valid and consistent).

Before objectively mapping the data, a mean, using the correlation scales to define the
weights for the sum (see Bretherton, etal, 1975) is removed.  The error estimate includes
the contributions from both the mapping and the mean estimation.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""