#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
- generate plots of the k(t) profiles from ASCII .dat files
- fit the kperp(t) & kparall(t) profiles (hyperbola) to 
  determine kdiff.
"""
from calc_k import kdiff
import numpy as np
import os
from funcs import calc_mfp, larmor
import console_colors as ccl
from os.path import isdir, isfile

PLAS    = os.environ['PLAS']

dir_plot    = './test'

Nm      = 128       # [1] n_modos
sig     = 1e0       # [1] (sigma/Bo)^2
Lc_2d   = 1e-2      # [AU]
Lc_slab = 1e-2      # [AU]
perc_slab = 0.2

Ek      = 1e8 #1e6 #6e5       # [eV] kinetic energy
atol    = 1e-6 #1e-5
dir_post = '%s/post/atol/Ek.%1.1eeV_atol.%1.0e' % (PLAS, Ek, atol)
#dir_info = '%s/output/Ek.%1.1eeV_rtol.%1.0e/info' % (PLAS, Ek, rtol)
dir_info = '%s/output/output_Ek.1e8eV_atol.1e-6_nB.50/info' % PLAS
assert isdir(dir_post) and isdir(dir_info), \
    " revisar nombres de directorios!!!"
#dir_info = '%s/output/Ek.%1.1eeV' % (PLAS, Ek) +\
#'/Nm%03d/slab%1.2f/sig.%1.1e' % (Nm, perc_slab, sig) +\
#'/Lc2D.%1.1e_LcSlab.%1.1e/info' % (Lc_2d, Lc_slab)

kd = kdiff(Ek, dir_info, dir_post, dir_plot)
t_decr      = 800 #TDECRs[i]  #200#200 #500 # 1000

kd.fit_kperp(t_decr)
kd.plot_kperp()
#kperp[i] = kd.kperp

kd.fit_kparall(t_decr=0.)
kd.plot_kparall()
#kparall[i] = kd.kparall

"""
#---------- Mean Free Paths:
MFPs_perp	= calc_mfp(kperp, Ek)	# [cm]
MFPs_parall	= calc_mfp(kparall, Ek)	# [cm]

data_out	= np.array([SLABs, kperp, kparall, MFPs_perp, MFPs_parall]).T
dir_out		= './test' #'../../post'
fname_out	= '%s/diff.pars_Nm%03d_Ek.%1.1e_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e__.dat' % (dir_out, Nm, Ek, sig, Lc_2d, Lc_slab)
np.savetxt(fname_out, data_out, fmt='%5.5e')
print ccl.Rn + " ------> hemos creado: " + fname_out + ccl.W
"""
#EOF
