#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
- genera ASCII .dat de los perfiles de k(t)
- *no* genera figuras
"""
import os
from os.path import isdir, isfile
import argparse

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


# NOTE: import this after the argparse, so we don't need to compile the Cython
# modules before check the help of this script.
from calc_k_vs_t import mfp_vs_t
import funcs 

#--- build 'post.h5'
mt = mfp_vs_t(pa.dir_src)
# build mean-free-path time profile
mt.calc_mfp_profile(moreinfo=True)
# build global-histogram of backscattering times
mt.build_TauHist()
# build global-histogram of backscattering orientations
mt.build_ThetaHist()
# histos for runtimes
mt.build_RuntimeHist(nbin=100)
# save all to .h5
mt.save2file(pa.dir_dst)


#--- fits && plots && append to .h5
mfp = funcs.mfp_mgr(
    dir_fig=pa.dir_fig, 
    fname_inp=mt.fname_out
)
mfp.fits_and_plots(pa.t_decr, pa.seed_b_perp, pa.seed_m_perp, pa.seed_m_para, pa.seed_m_para)
# append fit results to .h5 file
mfp.save2file(fname_out=mt.fname_out)

#EOF
