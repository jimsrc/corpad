#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
- genera ASCII .dat de los perfiles de k(t)
- *no* genera figuras
"""
from calc_k_vs_t import mfp_vs_t
import os
from os.path import isdir, isfile
import argparse

# retrieve args
parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', '--dir_src', 
    type=str, 
    default='out', 
    help='source dir, relative to repo root'
)
parser.add_argument(
    '-l', '--label', 
    type=str,
    default='__dummy__',
    help='basename for .h5 output file'
)
parser.add_argument(
    '-o', '--dir_dst', 
    type=str, 
    default='post/vs_r', 
    help='output dir, relative to repo root'
)
pa = parser.parse_args()

PLAS    = os.environ['PLAS']
dir_src = PLAS+'/'+pa.dir_src #'%s/out/r.0.30_NmS.128_Nm2d.256' % PLAS
dir_out = PLAS+'/'+pa.dir_dst #'%s/post/vs_r' % PLAS
if not(isdir(dir_out)):
    os.mkdir(dir_out)

assert isdir(dir_src) and isdir(dir_out), \
    " NO EXISTEN??: \n"+ dir_src + '\n' + dir_out

mt = mfp_vs_t(dir_src)
mt.calc_mfp_profile(dir_out, pa.label, moreinfo=True)

#EOF
