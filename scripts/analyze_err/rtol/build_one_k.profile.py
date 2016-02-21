#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from calc_k import kdiff
import numpy as np
import os

PLAS    = os.environ['PLAS']
dir_post     = '../../../post'
dir_plot    = '.'

Nm      = 128
sig     = 1e0
Lc_2d   = 3e-3
Lc_slab = 3e-2
Ek      = 6e5       # [eV]
perc_slab   = 0.2

dir_info    = '%s/output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e/info' % (PLAS, Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)

kd = kdiff(Ek, dir_info, dir_post, dir_plot)

t_decr      = 800 #TDECRs[i]  #200#200 #500 # 1000
kd.fit_kperp(t_decr)
kd.plot_kperp()

kd.fit_kparall(t_decr=0.)
kd.plot_kparall()
#EOF
