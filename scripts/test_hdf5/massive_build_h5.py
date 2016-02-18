#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from lib_build_h5 import build_h5
import os
"""
PLAS     = os.environ['PLAS']
dir_src  = '%s/output/Ek.1.0e+10eV_ii' % PLAS
dir_dst  = dir_src #'.'
ntot_B   = 50
ntot_pla = 122
fname_out_base = 'Ek.1.0e+10eV_ii.h5'
#
PLAS     = os.environ['PLAS']
dir_src  = '%s/output/Ek.1.0e+06eV/Nm128/slab0.20/sig.1.0e+00/Lc2D.1.0e-02_LcSlab.1.0e-02' % PLAS
dir_dst  = dir_src #'.'
ntot_B   = 50
ntot_pla = 242 #122
fname_out_base = 'Ek.1.0e+06eV__Nm128__slab0.20__sig.1.0e+00__Lc2D.1.0e-02_LcSlab.1.0e-02.h5'
#"""
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

PLAS     = os.environ['PLAS']
dir_src  = '%s/output' % PLAS + '/' + name_Ek
dir_src  += '/' + name_Nm + '/' + name_slab + '/' + name_sig + '/' + name_Lc
dir_dst  = dir_src #'.'
ntot_B   = 50
ntot_pla = 122 #242 #122
fname_out_base = '%s__%s__%s__%s__%s.h5' % (name_Ek, name_Nm, name_slab, name_sig, name_Lc)

print " dir_src: " + dir_src
if raw_input('proceed? ([y]/n): ')=='n': raise SystemExit

fname_out, nok, nbad = build_h5(ntot_B, ntot_pla, dir_src, dir_dst, fname_out_base)
#EOF
