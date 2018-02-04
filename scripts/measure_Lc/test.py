#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
import argparse, h5py
import funcs 

#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-bidim', '--bidim',
action='store_true',
default=False,
help='whether to do the bidimensional map of R-function'+\
' (ignored by default).',
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
'-nth', '--nth',
type=int,
default=7,
help='number theta values (angle w/ respect to z) for R(r,th).',
)
parser.add_argument(
'-ndr', '--ndr',
type=int,
default=32,
help='number dr values (distance from origin ro) for R(ro,ro+dr). The more values, the more resolution for a fixed maximum max(dr); see also option --ndrLc.',
)
parser.add_argument(
'-nLc', '--nLc',
type=int,
default=20,
help='th "origin" position xo will go until max(xo)=nLc*Lc.',
)
parser.add_argument(
'-ndrLc', '--ndrLc',
type=float,
default=2.5,
help='the displacements dr of R(ro,ro+dr) will go until max(dr)=ndrLc*Lc_slab. For more resolution, see option --ndr.',
)
parser.add_argument(
'-Lc', '--Lc',
type=float,
default=1.0,
help='correlation lenght (the same for slab & 2D models)',
)
parser.add_argument(
'-sigma', '--sigma',
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

# compiled Cython libraries
# NOTE: these must be compiled by hand (make clean && make) in 
# their respective directories.
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia

#AUincm = 1.5e13     # [cm]
xyz = np.zeros(3, dtype=np.float32)

#--- set B-turbulence model
pd = {
'Nm_slab'       : 128,
'Nm_2d'         : 128,
'lmin_s'        : 5e-5, #[lmin_s/Rl] 
'lmax_s'        : 1.0,  #[lmax_s/Rl] 
'lmin_2d'       : 5e-5, #[lmin_2d/Rl] 
'lmax_2d'       : 1.0,  #[lmax_2d/Rl] 
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

fl = funcs.LcMeasure(pd, pa.bidim)

dr   = np.linspace(0., pa.ndrLc*pd['Lc_slab'], pa.ndr) # all positions displacements
#nB_total = pa.Nrlz # number of turbulence realizations

# intialize R-functions
if not(pa.bidim):
    Rxx_perp, Ryy_perp = np.zeros((2,pa.Nrlz,dr.size))
    Rxx_para, Ryy_para = np.zeros((2,pa.Nrlz,dr.size))
    #--- run
    for nB in range(pa.Nrlz):
        print nB
        Rxx_perp[nB,:], Ryy_perp[nB,:] = fl.one_R_realiz(Nro=20, dr=dr, nB=nB, direcc='perp')
        Rxx_para[nB,:], Ryy_para[nB,:] = fl.one_R_realiz(Nro=20, dr=dr, nB=nB, direcc='parall')

    #--- normalize the perpendicular R-functions (i.e. R(r_perp)
    Rxx_perp /= Rxx_perp[:,0].mean()
    Ryy_perp /= Ryy_perp[:,0].mean()
    #--- normalize the parallel R-functions (i.e. R(r_parallel))
    Rxx_para /= Rxx_para[:,0].mean()
    Ryy_para /= Ryy_para[:,0].mean()

    #--- save to file
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

else: # pa.bibim==True
    dth  = np.linspace(0., np.pi/2, pa.nth) # polar angle (respect to z-axis)
    _Rxx, _Ryy = np.zeros((2,pa.Nrlz,dr.size,dth.size))
    #--- run
    for nB in range(pa.Nrlz):
        print nB
        _Rxx[nB,:,:], _Ryy[nB,:,:] = fl.one_R_realiz(Nro=pa.nLc,dr=dr,dth=dth,nB=nB)

    #--- normalize the perpendicular R-functions (i.e. R(r_perp)
    Rxx = _Rxx/_Rxx[:,0,:].mean()
    Ryy = _Ryy/_Ryy[:,0,:].mean()

    #--- save to file
    fo = h5py.File(pa.out, 'w')
    for nm, var in pd.iteritems():
        fo['psim/'+nm] = var
    #--- more parameters
    fo['psim/Nrlz'] = pa.Nrlz

    # now the physical results
    fo['dr']        = dr
    fo['dth']       = dth
    fo['Rxx']       = Rxx
    fo['Ryy']       = Ryy
    fo.close()
    #import pdb; pdb.set_trace()

print('\n [+] we generated: %s\n' % pa.out)

#EOF
