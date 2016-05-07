#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
- genera ASCII .dat de los perfiles de k(t)
- *no* genera figuras
"""
from calc_k_vs_t import k_vs_t
import os
from os.path import isdir, isfile

Ek	 = 1e9 #6e5	# [eV]
Bo   = 5e-5
atol = 1e-6

#atol = 1e-6
PLAS = os.environ['PLAS']
dir_data = '%s/output/output_Ek.1e9eV_atol.1e-6_nB.50' % PLAS
#dir_data = '%s/output/Ek.%1.1eeV_rtol.%1.0e' % (PLAS, Ek, rtol)
dir_out  = '%s/post/atol/Ek.%1.1eeV_atol.%1.0e' % (PLAS, Ek, atol)
if not(isdir(dir_out)):
    os.mkdir(dir_out)
"""
Nm		= 128
perc_slab	= 0.2 #0.02 #0.05 #0.10 # 0.00, 0.20, 0.40, 0.60, 1.00
sig		= 1e0
Lc_2d		= 1e-3
Lc_slab		= 1e-2
dir_data	= '%s/output' % PLAS +\
'/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e' % (Ek, Nm, perc_slab, sig) +\
'/Lc2D.%1.1e_LcSlab.%1.1e' % (Lc_2d, Lc_slab)
dir_out  = '%s/post/rtol/Ek.%1.1eeV_rtol.%1.0e' % (PLAS, Ek, rtol)
"""
assert isdir(dir_data) and isdir(dir_out), \
    " NO EXISTEN??: \n"+ dir_data + '\n' + dir_out

kt = k_vs_t(Ek, dir_data)
kt.calc_k_versus_t_ii(Bo, dir_out)

#EOF
