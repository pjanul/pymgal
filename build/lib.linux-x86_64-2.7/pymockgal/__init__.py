"""
pymockgal is a package written primarily in Python. It is designed to mimic
observed galaxies from hydro-dynamical simulations.

More details of the functions and models.
"""

__author__ = 'Weiguang Cui'
__email__ = 'cuiweiguang@gmail.com'
__ver__ = '1.0'

# import numpy as np  # For modern purposes
from SSP_models import SSP_models
from load_data import load_data
from filters import filters
import utils
# import astro_filter
# import wrapper
# import sfhs
# import weight
import dusts

# ezsps = ezsps.ezsps
# model = ezsps
# sspmodel = SSP_models.SSP_model
# loaddata = load_data.load_data
# filters = filters.filters
# # ezsps_light = ezsps_light.ezsps_light
# dust = dusts.charlot_fall
