#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia
import argparse, h5py

#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-o', '--out',
type=str,
default='./test.h5',
help='HDF5 output for R-function profiles',
)
parser.add_argument(
'-Nrlz', '--Nrlz',
type=int,
default=256,
help='number of turbulence realizations.',
)
parser.add_argument(
'-Lc', '--Lc',
type=float,
default=1.0,
help='correlation lenght (the same for slab & 2D models)',
)
parser.add_argument(
'-sig', '--sigma',
type=float,
default=0.3,
help='turbulence energy; i.e. (dB/Bo)^2.',
)
parser.add_argument(
'-slab', '--slab',
type=float,
default=0.2,
help='slab energy; i.e. (dB_slab/dB)^2.',
)
pa = parser.parse_args()

#AUincm = 1.5e13     # [cm]
xyz = np.zeros(3, dtype=np.float32)

#--- set B-turbulence model
pd = {
'Nm_slab'       : 128,
'Nm_2d'         : 128,
'lmin_s'        : 5e-3*pa.Lc, #[lmin_s/Rl] 
'lmax_s'        : 1e+2*pa.Lc,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-3*pa.Lc, #[lmin_2d/Rl] 
'lmax_2d'       : 1e+2*pa.Lc,  #[lmax_2d/Rl] 
'Lc_slab'       : pa.Lc,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': pa.sigma, # [1] fluctuation energy
'ratio_slab'    : pa.slab, # [1] (energy_slab)/(energy_total)
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
nB_total = pa.Nrlz # number of turbulence realizations
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


fo = h5py.File(pa.out, 'w')
for nm, var in pd.iteritems():
    fo['psim/'+nm] = var
#--- more parameters
fo['psim/Nrlz'] = pa.Nrlz

# now the physical results
fo['dr']       = dr
fo['Rxx_perp'] = Rxx_perp
fo['Ryy_perp'] = Ryy_perp
fo['Rxx_para'] = Rxx_para
fo['Ryy_para'] = Rxx_para
fo.close()

#EOF
