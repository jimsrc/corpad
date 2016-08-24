#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
execute as:
./massive_mfp.py -- -p h. -id 1 -sm_pe 190. -sb_pe 0.1 -sm_pa -16000 -sb_pa 30. -td 1000.
"""
import argparse
import funcs as ff
import os


# retrieve args
parser = argparse.ArgumentParser()
parser.add_argument(
    '-df', '--dir_fig', 
    type=str, 
    default='None', 
    help='output dir (relative to repo root)',
)
parser.add_argument(
    '-f', '--fname_inp', 
    type=str,
    default='None',
    help='file name for .h5 output file (excluding the path)',
)
parser.add_argument(
    '-p', '--prefix',
    type=str,
    default='DummyPrefix.',
    help='prefix for output-filename',
)
parser.add_argument(
    '-id', '--myid', 
    type=int,
    default=999,
    help='id-number for the output',
)
parser.add_argument(
    '-sm_pe', '--seed_m_perp',
    type=float,
    default=1.0,
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
    default=1.0,
    help='seed value \'m\' for hyperbola fit-function',
)
parser.add_argument(
    '-sb_pa', '--seed_b_para',
    type=float,
    default=0.1,
    help='seed value \'b\' for hyperbola fit-function',
)
parser.add_argument(
    '-td', '--t_decr',
    type=float,
    default=100.,
    help='value \'t_decr\' for hyperbola fit-function. It gives '+\
         'the low threshols value, from which the data is fitted.',
)

pa = parser.parse_args()

#dir_src = os.environ['PLAS']+'/'+ pa.dir_src
mfp = ff.mfp_mgr(pa.dir_fig, pa.fname_inp, pa.prefix, pa.myid)
print mfp.psim

mfp.fits_and_plots(pa.t_decr, pa.seed_b_perp, pa.seed_m_perp, pa.seed_m_para, pa.seed_m_para)

#EOF
