#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from h5py import File as f5
import os
from subprocess import check_output

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

PLAS    = os.environ['PLAS']
dir_src = '%s/output/Ek.1.0e+10eV_ii' % PLAS
dir_dst = '.'
fname_out = '%s/test.h5' % dir_dst
fo      = f5(fname_out, 'w')
nok, nbad = 0, 0

ntot_B   = 50
ntot_pla = 122
for iB in range(ntot_B):
    for ipla in range(ntot_pla):
        fname_traj = '%s/B%02d_pla%03d_traj.dat' % (dir_src, iB, ipla)
        ok1 = os.path.isfile(fname_traj) # trajectory file
        fname_misc = '%s/B%02d_pla%03d_misc.dat' % (dir_src, iB, ipla)
        ok2 = os.path.isfile(fname_misc) # misc file
        if (ok1 and ok2):
            # info on the trajecctory (see funcs.cc for units)
            t, x, y, z, mu, err = np.loadtxt(fname_traj).T
            path = 'B%02d/pla%03d' % (iB, ipla)
            fo[path+'/t']   = t                   # [sec] 
            fo[path+'/xyz'] = np.array([x, y, z]) # [AU]
            fo[path+'/mu']  = mu                  # [1]
            fo[path+'/err'] = err                 # [1] gamma/scl.gamma - 1.

            # histogram on scatterings
            tau, nreb = np.loadtxt(fname_misc, skiprows=5).T
            fo[path+'/hist_tau.coll'] = np.array([tau, nreb])

            # save misc parameters
            par = read_params(fname_misc)
            for name in par.keys():
                if name=='trun/min': name_='trun-in-min' #runtime [min]
                else: name_ = name
                fo[path+'/'+name_] = par[name]

            # save modification time
            command  = 'date -r ' + fname_traj + ' -u +%d-%m-%Y\ %H:%M:%S'
            mod_time = check_output(command, shell=True)
            fo[path+'/modification-time'] = mod_time[:-1] # -1 come el '\n'

            nok += 1
            print " --> " + fname_traj

        else:
            assert not(ok1 or ok2), " ----> Existe traj o misc: RIDICULO!!!... wierd."
            nbad += 1

fo['nok']   = nok   # total nmbr of simulated plas
fo['nbad']  = nbad  # plas that were *not* simulated (but were suposed to)
fo.close() 
print "\n dir_src: " + dir_src
print " dir_dst: " + dir_dst + '\n'
assert nok>0, "\n ---> NO LEIMOS NADA!!!\n"  # just a check
print "\n ---> generated: " + fname_out + '\n'
#EOF
