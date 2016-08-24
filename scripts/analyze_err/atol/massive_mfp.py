#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import argparse
import funcs as ff
import os


# retrieve args
parser = argparse.ArgumentParser()
parser.add_argument(
    '-ds', '--dir_src', 
    type=str, 
    default='out', 
    help='source dir (relative to repo root)'
)
parser.add_argument(
    '-df', '--dir_fig', 
    type=str, 
    default='plots', 
    help='output dir (relative to repo root)'
)
parser.add_argument(
    '-f', '--fname_inp', 
    type=str,
    default='post.h5',
    help='file name for .h5 output file (excluding the path)'
)
parser.add_argument(
    '-o', '--olabel', 
    type=str,
    default='test_output',
    help='output basename for .png output file'
)
parser.add_argument(
    '-sm_pe', '--seed_m_perp',
    type=float,
    default=1.0,
    help='seed value \'m\' for hyperbola fit-function'
)
parser.add_argument(
    '-sb_pe', '--seed_b_perp',
    type=float,
    default=0.1,
    help='seed value \'b\' for hyperbola fit-function'
)
parser.add_argument(
    '-sm_pa', '--seed_m_para',
    type=float,
    default=1.0,
    help='seed value \'m\' for hyperbola fit-function'
)
parser.add_argument(
    '-sb_pa', '--seed_b_para',
    type=float,
    default=0.1,
    help='seed value \'b\' for hyperbola fit-function'
)
parser.add_argument(
    '-td', '--t_decr',
    type=float,
    default=100.,
    help='value \'t_decr\' for hyperbola fit-function. It gives '+\
         'the low threshols value, from which the data is fitted.'
)

pa = parser.parse_args()

dir_src = os.environ['PLAS']+'/'+ pa.dir_src
mfp = ff.mfp_mgr(pa.dir_fig, dir_src, pa.fname_inp, pa.olabel)
print mfp.psim

mfp.fit_perp(pa.t_decr, pa.seed_b_perp, pa.seed_m_perp)
mfp.plot_perp()
mfp.fit_parall(pa.t_decr, pa.seed_b_para, pa.seed_m_para)
mfp.plot_parall()

#EOF
