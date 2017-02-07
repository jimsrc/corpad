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

dr   = np.linspace(0., 6.*pd['Lc_slab'], 1000)
BiBj = fl.one_R_realiz(Nro=20, dr=dr, ij=(0,0), pd=pd, nB=0)

#print fl.m.Bxyz(xyz)
#EOF
