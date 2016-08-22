#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import norm
from src_Bmodel import Bmodel
#from Bparker.Bparker import return_B as Bparker_vector
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia
from numpy.random import RandomState
#--- parameters
#from params import (
#    nB, pd, psim, pother, mu, ph, AU_in_cm
#)

nB = 0 # id of B-realization
AUincm = 1.5e13     # [cm]

#xyz = np.zeros(3, dtype=np.float32)

m  = Bmodel.Bmodel() # model slab/2D
ro = 1.0                   # [AU]
Lc_slab = Lc_memilia(r=ro) # [AU]
Bo      = Bo_parker(r=ro)  # [G]
Rl = calc_Rlarmor(
    rigidity=1e9,       # [V] rigidity
    Bo=Bo,              # [G] Bo-field
    )/AUincm            # [AU]

##--- sample of positions
#Xs, Ys = 200.*RS.random_sample((2,n))
#Zs     = 3500.*RS.random_sample(n)

#--- set B-turbulence model
pd = {
'Nm_slab'       : 256,
'Nm_2d'         : 256,
'lmin_s'        : 5e-5/Rl, #[lmin_s/Rl] 
'lmax_s'        : 1.0/Rl,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-5/Rl, #[lmin_2d/Rl] 
'lmax_2d'       : 1.0/Rl,  #[lmax_2d/Rl] 
'Lc_slab'       : Lc_slab/Rl,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': 1.0, # [1] fluctuation energy
'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
'sem_slab0'     : 17,
'sem_slab1'     : 101,
'sem_slab2'     : 33,
'sem_two0'      : 14,
'sem_two1'      : 79,
}
m.set_Bmodel(pdict=pd, nB=nB)
#print m.Bxyz(xyz)

nhist = 100
hc = np.zeros(nhist, dtype=np.int32)
th_lo, th_hi = 0., 90.
hx = np.linspace(th_lo, th_hi, nhist)

RS = RandomState(seed=1)
n = 100000
#for xyz in zip(Xs,Ys,Zs):
for i in range(n):
    x,y = 200.*RS.random_sample(2)
    z   = 3500.*RS.random_sample(1)
    Bx,By,Bz = m.Bxyz([x,y,z])
    bz = Bz/norm([Bx,By,Bz])
    th = np.arccos(bz)*180./np.pi
    hc += np.histogram(th, bins=nhist, range=(th_lo,th_hi))[0]


"""
Rs = np.arange(0.2, 5., 0.05)
x[1] = np.pi/2. # eliptica

Bm = []
for x[0] in Rs:
    B = bp.return_B(x)
    Bm += [ np.square(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]) ]

Bm = np.array(Bm)
"""
#EOF
