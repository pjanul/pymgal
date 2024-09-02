"""
pymgal is a package written primarily in Python. It is designed to mimic
observed galaxies from hydro-dynamical simulations.

More details of the functions and models.
"""

# import numpy as np  # For modern purposes
from pymgal.SSP_models import SSP_models
from pymgal.load_data import load_data
from pymgal.filters import filters
from pymgal import utils
from pymgal.projection import projection
from pymgal import dusts
from pymgal import __version__

__author__ = ['Weiguang Cui', 'Patrick Janulewicz']
__email__ = ['cuiweiguang@gmail.com', 'patrick.janulewicz@gmail.com']
__version__ = __version__.__version__
__all__ = ["pymgal"]
# ezsps = ezsps.ezsps
# model = ezsps
# sspmodel = SSP_models.SSP_model
# loaddata = load_data.load_data
# filters = filters.filters
# # ezsps_light = ezsps_light.ezsps_light
# dust = dusts.charlot_fall
