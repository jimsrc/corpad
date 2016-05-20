#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
from scipy.interpolate import (
    splrep,  # spline interpolation
    splev    # spline evaluator
)

fnm_lam = './fig01_bavassano_etal00_lambda.dat'
fnm_rad = './fig01_bavassano_etal00_radius.dat'

tl, l = np.loadtxt(fnm_lam).T
tr, r = np.loadtxt(fnm_rad).T


time = np.linspace(95., 97.0, 100)
#--- interp lambda
tck = splrep(tl, l, s=0)
lam = splev(time, tck, der=0) # der: order of derivative to compute (less than "k")
#--- interp radius
tck = splrep(tr, r, s=0)
rad = splev(time, tck, der=0)


