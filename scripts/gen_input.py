#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#import src.cython_wrapper as cw
import shared.funcs as sf
from pylab import pause
#--- parameters
#from params import (
#    nB, pd, psim, pother, mu, ph, AU_in_cm
#)
#from mpi4py import MPI
from h5py import File as h5
from os.path import isfile, isdir
import os, sys
from glob import glob
from shared.Bparker.Bparker import calc_Rlarmor

#--- globals
AUincm = 1.5e13             # [cm]

def write(f, dc, name):
    ostr = str(dc[name])+'\t'+name+'\n'
    f.write(ostr)

class fwrapper(object):
    def __init__(self,fname, mode, dc):
        assert mode=='w' or mode=='r+', " only write mode!"
        self.f = open(fname, mode)
        self.dc = dc

    def write(self, name):
        ostr = str(self.dc[name])+'\t'+name+'\n'
        self.f.write(ostr)

    def __del__(self):
        print " we wrote: "+self.f.name
        self.f.close()


"""
(Bo = 5nT; omega=omega_{rel})
Ek/eV   Rigidity/V   Rl/AU          beta
1e6     4.33306E+07  1.929785E-04   4.613215E-02
1e7     1.37352E+08  6.117125E-04   1.448440E-01
1e8     4.44583E+08  1.980009E-03   4.281955E-01
1e9     1.69604E+09  7.553521E-03   8.750257E-01  (*)
1e10    1.0898E+10   4.853544E-02   9.963142E-01
--- Bparker model
(Ek=1e9eV)
r[AU]    B[nT]       Rl[AU]         Lc[AU]      Rl/Lc   Rl/(5e-5AU)
0.2      91.0542857  4.147812E-04   0.0044549   0.09    8.30
0.3      41.4298571  9.116035E-04   0.0053034   0.17    18.23
0.4      24.041      1.570966E-03   0.0060017   0.26    31.42
0.5      15.972      2.364613E-03   0.0066061   0.36    47.29
0.7      8.897       4.244982E-03   0.0076345   0.56    84.90
0.8      7.14642857  5.284822E-03   0.0080857   0.65    105.70
1.0      5.0         7.553521E-03   0.0089      0.85    151.07
2.0      1.99653571  1.891657E-02   0.0119904   1.58    378.33
"""
ro = 0.5 # (0.2, 0.3, 0.5, 0.7, 0.9)
lc = sf.Lc_memilia(r=ro)   # [AU]
Lc_slab = lc
Rl = calc_Rlarmor(
    rigidity=1.69604E+09,    # [V]
    Bo=sf.Bo_parker(r=ro)    # [Gauss]
    )/AUincm                 # [AU] Larmor radii
#--- set B-turbulence model
pturb = {
'Nm_slab'       : 128,
'Nm_2d'         : 256,
'lmin_s'        : 5e-5/Rl, #[lmin_s/Rl] 
'lmax_s'        : ro/Rl,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-5/Rl, #[lmin_2d/Rl] 
'lmax_2d'       : ro/Rl,  #[lmax_2d/Rl] 
'Lc_slab'       : Lc_slab/Rl,  # in units of Larmor-radii
'xi'            : 1.0, # [1] xi=Lc_2d/Lc_slab 
'sigma_Bo_ratio': 0.3, # [1] fluctuation energy
'ratio_slab'    : 0.2, # [1] (energy_slab)/(energy_total)
#--- seeds
'sem.slab0'     : 17,
'sem.slab1'     : 101,
'sem.slab2'     : 33,
'sem.two0'      : 14,
'sem.two1'      : 79,
}
#--- corregimos input
eps_o = 1e-4 #3.33e-6 #3.33e-5 #1.0e-4 #3.3e-6 #4e-5 # ratio: (error-step)/(lambda_min)
lmin = np.min([pturb['lmin_s'], pturb['lmin_2d']]) # [1] smallest turb scale
ppla = {
'frac_gyroperiod'   : 5e-2,
'tmax'              : 4e4,
'npoints'           : 20,
'atol'              : lmin*eps_o,  # [1]
'rtol'              : 0.0, #1e-6
'nB'                : 50,
'time_scale'        : 'mixed',
'tmaxHistTau'       : 150.0,
'nHistTau'          : 150,
'nThetaColl'        : 150,
}

#---- create turb-file
ft = fwrapper('../inputs/turb.in', 'w', pturb)
ft.write('Nm_slab')
ft.write('Nm_2d')
ft.write('lmin_s')
ft.write('lmax_s')
ft.write('lmin_2d')
ft.write('lmax_2d')
ft.write('Lc_slab')
ft.write('xi')
ft.write('ratio_slab')
ft.write('sigma_Bo_ratio')
ft.write('sem.slab0')
ft.write('sem.slab1')
ft.write('sem.slab2')
ft.write('sem.two0')
ft.write('sem.two1')
del ft # close

fp = fwrapper('../inputs/plas.in', 'w', ppla)
fp.write('frac_gyroperiod')
fp.write('tmax')
fp.write('npoints')
fp.write('atol')
fp.write('rtol')
fp.write('nB')
fp.write('time_scale')
fp.write('tmaxHistTau')
fp.write('nHistTau')
fp.write('nThetaColl')
del fp

#EOF
