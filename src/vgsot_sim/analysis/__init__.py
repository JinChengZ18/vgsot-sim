"""Analysis sub-package: experiment / Monte-Carlo / fitting utilities.

Modules here re-implement the closed-form and statistical tools previously
scattered across the 04 / 06 / 07 / 08 chapter scripts as a single coherent
API. The scripts now delegate the math to these modules and only handle
data loading and plotting.

Public surface (re-exported below):

- nb_fit:       Néel-Brown log-linear fits from hysteresis loops
- sigmoid_fit:  4-parameter Sigmoid fit + η_c calibration vs NB
- variability:  Brinkman-decomposed PDK variability budget + D2D MC sweep
- sampling:     Wilson interval, exact binomial coverage, MC sampling sensitivity
- materials:    T-dependent M_s/K_i/η helpers (thin re-export of material_temperature)
"""
from . import nb_fit, sigmoid_fit, variability, sampling
from .. import material_temperature as materials

__all__ = ["nb_fit", "sigmoid_fit", "variability", "sampling", "materials"]
