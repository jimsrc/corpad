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
import funcs 

# retrieve args
parser = argparse.ArgumentParser()
parser.add_argument(
    '-di', '--dir_src', 
    type=str, 
    default='out', 
    help='source dir, relative to repo root'
)
#parser.add_argument(
#    '-l', '--label', 
#    type=str,
#    default='__dummy__',
#    help='basename for .h5 output file'
#)
parser.add_argument(
    '-do', '--dir_dst', 
    type=str, 
    default='post/vs_r', 
    help='output dir, relative to repo root'
)
parser.add_argument(
    '-df', '--dir_fig', 
    type=str, 
    default='plots', 
    help='output dir (relative to repo root)',
)
parser.add_argument(
    '-sm_pe', '--seed_m_perp',
    type=float,
    default=190.,
    help='seed value \'m\' for hyperbola fit-function',
)
parser.add_argument(
    '-sb_pe', '--seed_b_perp',
    type=float,
    default=0.1,
    help='seed value \'b\' for hyperbola fit-function',
)
parser.add_argument(
    '-sm_pa', '--seed_m_para',
    type=float,
    default=-1.6e4,
    help='seed value \'m\' for hyperbola fit-function',
)
parser.add_argument(
    '-sb_pa', '--seed_b_para',
    type=float,
    default=30.,
    help='seed value \'b\' for hyperbola fit-function',
)
parser.add_argument(
    '-td', '--t_decr',
    type=float,
    default=1e3,
    help='value \'t_decr\' for hyperbola fit-function. It gives '+\
         'the low threshols value, from which the data is fitted.',
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
# build mean-free-path time profile
mt.calc_mfp_profile(moreinfo=True)
# build total-histogram of backscattering times
mt.build_TauHist()
# save all to .h5
mt.save2file(dir_out)

#--- fits && plots
mfp = funcs.mfp_mgr(
    dir_fig=pa.dir_fig, 
    fname_inp=mt.fname_out
)
#print mfp.psim

mfp.fits_and_plots(pa.t_decr, pa.seed_b_perp, pa.seed_m_perp, pa.seed_m_para, pa.seed_m_para)

#EOF
