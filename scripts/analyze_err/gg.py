#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from h5py import File as f5
import os

def read_pla(iB=0, ipla=0):
    path = 'B%02d/pla%03d' % (iB, ipla)
    t   = f[path+'/t'][...]
    mu  = f[path+'/mu'][...]
    err = f[path+'/err'][...]
    return t, mu, err

PLAS     = os.environ['PLAS']
Ek       = 1e9 #1e8 #1e7 #1e6
slab     = 0.2
Nm       = 128
sig      = 1.0
Lc2d, LcSlab = 1e-2, 1e-2

name_Ek = 'Ek.%1.1eeV' % Ek
name_Nm = 'Nm%03d' % Nm
name_slab = 'slab%2.2f' % slab
name_sig  = 'sig.%1.1e' % sig
name_Lc   = 'Lc2D.%1.1e_LcSlab.%1.1e' % (Lc2d, LcSlab)

dir_src  = '%s/output' % PLAS + '/' + name_Ek
dir_src  += '/' + name_Nm + '/' + name_slab + '/' + name_sig + '/' + name_Lc
fname_inp_base = '%s__%s__%s__%s__%s.h5' % (name_Ek, name_Nm, name_slab, name_sig, name_Lc)
fname_inp = '%s/%s' % (dir_src, fname_inp_base)

#"""
dir_src = '../../output/Ek.1e6eV_rtol.1e-6'
dir_dst = dir_src
fname_inp = '%s/Ek.1e6eV_rtol.1e-6.h5' % dir_src
#"""

f = f5(fname_inp, 'r')

iB, ipla = 0, 0
t, mu, err = read_pla(iB, ipla)
#EOF
