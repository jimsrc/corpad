#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from h5py import File as h5
import os


def read_params(fname):
    f = open(fname_misc, 'r')
    par = {} #output
    for i in range(10): # esta dentro de las primeras 10 lineas
        l = f.readline().split()
        if l[0]=='#####':
            break

        name = l[0][:-1] # -1 para comernos el ":"
        value = np.float(l[1])
        par[name] = value

    return par


dir_src = '../Ek.1.0e+10eV_ii'
dir_dst = '.'
fname_out = '%s/test.h5' % dir_dst
fo      = h5(fname_out, 'w')
nok, nbad = 0, 0

ntot_B   = 50
ntot_pla = 122
for iB in range(ntot_B):
    for ipla in range(ntot_pla):
        fname_traj = '%s/B%02d_pla%03d_traj.dat' % (dir_src, iB, ipla)
        ok1 = os.path.isfile(fname_traj)
        fname_misc = '%s/B%02d_pla%03d_misc.dat' % (dir_src, iB, ipla)
        ok2 = os.path.isfile(fname_misc)
        if (ok1 and ok2):
            t, x, y, z, mu, err = np.loadtxt(fname_traj).T
            path = 'B%02d/pla%03d' % (iB, ipla)
            fo[path+'/t'] = t
            fo[path+'/xyz'] = np.array([x, y, z])
            fo[path+'/mu'] = mu
            fo[path+'/err'] = err

            tau, nreb = np.loadtxt(fname_misc, skiprows=5).T
            fo[path+'/hist_tau.coll'] = np.array([tau, nreb])

            par = read_params(fname_misc)
            for name in par.keys():
                fo[path+'/'+name] = par[name]

            nok += 1
            print " --> " + fname_traj

        else:
            assert not(ok1 or ok2), " ----> Existe traj o misc: RIDICULO!!!... wierd."
            nbad += 1

fo.close()
