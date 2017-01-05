#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import plot, close, show, figure
import argparse
import src_Bmodel.Bmodel as bmt # turbulence B-model

#--- read parameters into pd={}
pd = {
'Nm_slab'   : 128,
'Nm_2d'     : 128,
'lmin_s'    : 5e-3,
'lmax_s'    : 1e2,
'lmin_2d'   : 5e-3,
'lmax_2d'   : 1e2,
'Lc_slab'   : 1.,
'xi'        : 1.,
'sigma_Bo_ratio' : 0.3,
'ratio_slab': 0.5,
#TODO: the seeds are not in this .h5 :(
'sem_slab0' : 17,
'sem_slab1' : 101,
'sem_slab2' : 33,
'sem_two0'  : 14,
'sem_two1'  : 79,
}

#--- load it into the B-model module
Bm = bmt.Bmodel()
Bm._build_pturb(pd=pd)
Bm._build_par(nB=0)

#--- read the generated spectra (Bk_SLAB, Bk_2D) && plot
pnms = ['k_s', 'k_2d', 'Bk_SLAB', 'Bk_2D', 'dk_s', 'dk_2d', 'Lc_slab', 'Lc_2d']
p = {}
for nm in pnms:
    p[nm] = Bm.read_param(nm)

Bk_2D   = p['Bk_2D']
Bk_SLAB = p['Bk_SLAB']
k_s     = p['k_s']
k_2d    = p['k_2d']
dk_s    = p['dk_s']
dk_2d   = p['dk_2d']
Lc_slab = p['Lc_slab']
Lc_2d  = p['Lc_2d']


fig = figure(1, figsize=(4,6))
ax  = fig.add_subplot(111)

ax.loglog(k_2d, Bk_2D, label='2D')
ax.loglog(k_s, Bk_SLAB, label='SLAB')

ax.grid(True)
ax.legend(loc='best')
show()
close(fig)

#EOF
