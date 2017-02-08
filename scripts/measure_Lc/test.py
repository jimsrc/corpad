#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia

#AUincm = 1.5e13     # [cm]
xyz = np.zeros(3, dtype=np.float32)

#--- set B-turbulence model
pd = {
'Nm_slab'       : 128,
'Nm_2d'         : 128,
'lmin_s'        : 5e-3, #[lmin_s/Rl] 
'lmax_s'        : 1e+2,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-3, #[lmin_2d/Rl] 
'lmax_2d'       : 1e+2,  #[lmax_2d/Rl] 
'Lc_slab'       : 1.0,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': 0.3, # [1] fluctuation energy
'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
#--- seeds
'sem_slab0'     : 17,
'sem_slab1'     : 101,
'sem_slab2'     : 33,
'sem_two0'      : 14,
'sem_two1'      : 79,
}

import funcs 
fl = funcs.LcMeasure(pd)

dr   = np.linspace(0., 5.*pd['Lc_slab'], 128) # all positions displacements
nB_total = 256 # number of turbulence realizations
Rxx_perp, Ryy_perp = np.zeros((2,nB_total,dr.size))
Rxx_para, Ryy_para = np.zeros((2,nB_total,dr.size))

for nB in range(nB_total):
    print nB
    Rxx_perp[nB,:], Ryy_perp[nB,:] = fl.one_R_realiz(Nro=20, dr=dr, nB=nB, direcc='perp')
    Rxx_para[nB,:], Ryy_para[nB,:] = fl.one_R_realiz(Nro=20, dr=dr, nB=nB, direcc='parall')

#--- normalize the perpendicular R-functions (i.e. R(r_perp)
Rxx_perp /= Rxx_perp[:,0].mean()
Ryy_perp /= Ryy_perp[:,0].mean()
#--- normalize the parallel R-functions (i.e. R(r_parallel))
Rxx_para /= Rxx_para[:,0].mean()
Ryy_para /= Ryy_para[:,0].mean()

cc = dr<2.0
x = dr[cc]
y = Rxx_perp.mean(axis=0)[cc]

m, b = np.polyfit(x, np.log(y), deg=1, cov=False)
print m, b

#print fl.m.Bxyz(xyz)
#EOF
