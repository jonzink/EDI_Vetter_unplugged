#!/usr/bin/python
"""
EDI-Vetter_unplugged 

identifies false positive transit signals using TLS information and has been
simplified from the full EDI-Vetter algorithm for easy implementation with the TLS output.

Created by Jon Zink  (jzink@astro.ucla.edu)
Developed with python 3.7.1 

If you make use of this code, please cite: 
J. K. Zink et al. 2020a


"""

from EDIunplugged.main import parameters
from EDIunplugged.main import Go
import EDIunplugged.metrics