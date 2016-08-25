#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from lib_build_h5 import build_hdf5
import os
from os.path import isfile, isdir
import argparse

# retrieve the IDentifiers :-)
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir_src', type=str)
#parser.add_argument('-l', '--legend', type=str)
parser.add_argument('-o', '--fname_out', type=str, default='out.h5')
pa = parser.parse_args()
#PLAS     = os.environ['PLAS']
dir_src  = pa.dir_src #'../../output/r.0.20_NmS.128_Nm2d.256'
dir_dst  = dir_src #'.'
#ntot_B   = 50
#ntot_pla = 122 #242 #122
fname_out_base = pa.fname_out

#
print " dir_src: " + dir_src
if raw_input('proceed? ([y]/n): ')=='n': raise SystemExit

#fname_out, nok, nbad = build_h5_ii(ntot_B, ntot_pla, dir_src, dir_dst, fname_out_base)
bh5 = build_hdf5(dir_src, dir_dst, fname_out_base)
#bh5.get_trajs() # converts all .dat to a single .h5
bh5.grab_all() # grabs trajectories coordinates, tau-histos && theta-histos, and convert to .h5
#EOF
