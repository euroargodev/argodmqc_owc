"""
-----Covariance of 2D space and time (and potential vorticity)-----

Written by: A. Wong
When: 29/04/2005
Converted to python by: Edward Small
When: 13/11/2019

Calculates covariance of each point against every other point using the
Squared Exponential (SE) covariance function:

SE(x,y) = exp(-(x-y)^2/2l) where (x-y) is the difference between points (could be distance,
time, etc), and l is the characteristic length scale (how close points have to be
to influence each other significantly).

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

x = 1

for i in x:
    print(i)