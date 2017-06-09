"""
pymockgal is a package written primarily in Python. It is designed to mimic
observed galaxies from hydro-dynamical simulations.

More details of the functions and models.
"""

__author__ = 'Weiguang Cui'
__email__ = 'cuiweiguang@gmail.com'
__version__ = 'beta'

# import numpy as np  # For modern purposes
from pymgal.SSP_models import SSP_models
from pymgal.load_data import load_data
from pymgal.filters import filters
from pymgal import utils
from pymgal.projection import projection
from pymgal import dusts

# ezsps = ezsps.ezsps
# model = ezsps
# sspmodel = SSP_models.SSP_model
# loaddata = load_data.load_data
# filters = filters.filters
# # ezsps_light = ezsps_light.ezsps_light
# dust = dusts.charlot_fall
