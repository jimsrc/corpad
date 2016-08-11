#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from h5py import File as h5
import os
from os.path import isfile, isdir
from subprocess import check_output
from glob import glob

def read_params(fname):
    """ done for first header.
    TODO: we need to grab the 2nd header too!!
    """
    f = open(fname, 'r')
    par = {} #output
    for i in range(10): # esta dentro de las primeras 10 lineas
        l = f.readline().split()
        #print " ---> ", l
        number = u'%s' % l[-1] # presumably a number
        if not number.replace('.','').replace('-','').isnumeric():
            if l[0]=='#####':
                break
            else:
                continue # we proceed ONLY IF this is numeric string
        #print ' FIRST: ', l[0]
        if l[0]=='#####':
            #print "IM I HERE????"
            break # end of header

        name = l[1][:-1] # l[0] es '#', y -1 para comernos el ":"
        value = np.float(l[2]) # l[2] es el valor
        par[name] = value

    return par

def nans(sh, dt=np.float64):
    return np.nan*np.ones(sh, dtype=dt)


"""
PLAS    = os.environ['PLAS']
dir_src = '%s/output/Ek.1.0e+10eV_ii' % PLAS
dir_dst = '.'
ntot_B   = 50
ntot_pla = 122
"""

def build_h5_ii(ntot_B, ntot_pla, dir_src, dir_dst, fname_out_base):
    nskip_misc = 10 #5 # lineas para obviar en los misc-files
    assert isdir(dir_dst) or isdir(dir_src), \
        " ---> paths don't exist! :\n --> %s\n --> %s" % (dir_dst, dir_src)
    flist = glob(dir_src+'/B*traj*.dat')
    assert len(flist)>0, " ---> NO hay trayectorias! :("
    with open(flist[0], 'r') as ftemp:
        nt = len(ftemp.readlines()) # number of times

    flist = glob(dir_src+'/B*misc*.dat')
    with open(flist[0], 'r') as ftemp:
        ntau = len(ftemp.readlines())-nskip_misc # N of tau-histogram bins
    
    fname_out = '%s/%s' % (dir_dst, fname_out_base) #'%s/test.h5' % dir_dst
    fo        = h5(fname_out, 'w')
    nok, nbad = 0, 0
    for iB in range(ntot_B):
        o_x, o_y, o_z = nans((nt,ntot_pla)), nans((nt,ntot_pla)), nans((nt,ntot_pla))
        o_mu, o_err = nans((nt, ntot_pla)), nans((nt, ntot_pla))
        #o_tau = nans((ntau, ntot_pla))
        o_nreb = nans((ntau, ntot_pla))
        # misc parameters
        temp = read_params(flist[0])
        misc = {}
        for name in temp.keys():
            if name=='trun/min': 
                name_ = 'trun-in-min' #runtime [min]
            else: 
                name_ = name
            misc[name_] = nans(ntot_pla)
        misc['modification-time'] = np.empty(ntot_pla, dtype='S80')

        for ipla in range(ntot_pla):
            fname_traj = '%s/B%02d_pla%03d_traj.dat' % (dir_src, iB, ipla)
            ok1 = isfile(fname_traj) # trajectory file
            fname_misc = '%s/B%02d_pla%03d_misc.dat' % (dir_src, iB, ipla)
            ok2 = isfile(fname_misc) # misc file
            if (ok1 and ok2):
                # info on the trajecctory (see funcs.cc for units)
                t, x, y, z, mu, err = np.loadtxt(fname_traj).T
                o_x[:,ipla]   = x    # [AU]
                o_y[:,ipla]   = y    # [AU]
                o_z[:,ipla]   = z    # [AU]
                o_mu[:,ipla]  = mu   # [1]
                o_err[:,ipla] = err  # [1] gamma/scl.gamma-1.0

                # histogram of scatterings
                tau, nreb = np.loadtxt(fname_misc, skiprows=5).T
                if ipla==0:
                    o_tau = tau
                    time  = t

                o_nreb[:,ipla] = nreb
                path = 'B%02d/misc_params/pla%03d' % (iB, ipla)
                # save misc parameters
                par = read_params(fname_misc)
                for name in par.keys():
                    if name=='trun/min': 
                        name_='trun-in-min' #runtime [min]
                    else: 
                        name_ = name
                    #fo[path+'/'+name_] = par[name]
                    misc[name_][ipla] = par[name]

                # save modification time
                command  = 'date -r ' + fname_traj + ' -u +%d-%m-%Y\ %H:%M:%S'
                mod_time = check_output(command, shell=True)
                misc['modification-time'][ipla] = mod_time[:-1] # -1 come el '\n'
                nok += 1
                print " --> " + fname_traj

            else:
                assert not(ok1 or ok2), " ----> Existe traj o misc: RIDICULO!!!... wierd."
                nbad += 1

        path = 'B%02d' % iB
        fo[path+'/x'] = o_x
        fo[path+'/y'] = o_y
        fo[path+'/z'] = o_z
        fo[path+'/mu']  = o_mu
        fo[path+'/err'] = o_err
        for name in misc.keys():
            fo[path+'/misc/'+name] = misc[name]

    fo['time']  = time
    fo['nok']   = nok   # total nmbr of simulated plas
    fo['nbad']  = nbad  # plas that were *not* simulated (but were suposed to)
    fo['ntot_B'] = ntot_B # total nmbr of B realizations
    fo['ntot_pla'] = ntot_pla # total nmbr of plas
    fo['ntimes'] = nt
    fo['ntau'] = ntau
    fo.close() 
    print "\n dir_src: " + dir_src
    print " dir_dst: " + dir_dst + '\n'
    assert nok>0, "\n ---> NO LEIMOS NADA!!!\n"  # just a check
    print " nok/nbad: %d/%d" % (nok, nbad)
    print " ---> generated: " + fname_out

    return (fname_out, nok, nbad)


def build_h5(ntot_B, ntot_pla, dir_src, dir_dst, fname_out_base):
    assert os.path.isdir(dir_dst) or os.path.isdir(dir_src), \
            " ---> paths don't exist! :\n --> %s\n --> %s" % (dir_dst, dir_src)

    fname_out = '%s/%s' % (dir_dst, fname_out_base) #'%s/test.h5' % dir_dst
    fo        = h5(fname_out, 'w')
    nok, nbad = 0, 0

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
    fo['ntot_B'] = ntot_B # total nmbr of B realizations
    fo['ntot_pla'] = ntot_pla # total nmbr of plas
    fo.close() 
    print "\n dir_src: " + dir_src
    print " dir_dst: " + dir_dst + '\n'
    assert nok>0, "\n ---> NO LEIMOS NADA!!!\n"  # just a check
    print " nok/nbad: %d/%d" % (nok, nbad)
    print " ---> generated: " + fname_out

    return (fname_out, nok, nbad)


class build_hdf5(dir_src, dir_dst, fname_out_base):
    def __init__():
        assert isdir(dir_dst) or isdir(dir_src), \
            " ---> paths don't exist! :\n --> %s\n --> %s"%(dir_dst, dir_src)
        self.dir_src = dir_src
        self.dir_dst = dir_dst
        self.fname_out_base = fname_out_base
        flist = glob(dir_src+'/B*traj*.dat')
        assert len(flist)>0, " ---> NO hay trayectorias! :("
        self.psim = {}

    def extract_block(self, block_name='TRAJECTORY', fname_traj):
        lines = []
        Read  = False
        for line in open(fname_traj,'r'):
            if line.startswith('#BEGIN '+block_name):
                Read = True; continue
            if Read and line.startswith('#END'):
                Read = False; continue
            if flag:
                lines += [line]
        return lines

    def get_traj(self, fname_traj):
        block = self.extract_block('TRAJECTORY', fname_traj)
        #--- read trajectory
        i = 0
        for line in block:
            if not line.startswith('#'):
                #t[i],x[i],y[i],z[i],mu[i],err[i]=map(float,line.split(' '))
                num_trj += [ map(float,line.split(' ')) ]
        trj = np.array(num_trj)
        return trj

    def get_trajs(self):
        nB = len(glob(self.dir_src+'/B*_traj000.dat'))
        nP = len(glob(self.dir_src+'/B00_traj*.dat'))
        flist = glob(dir_src+'/B*traj*.dat')
        #TODO: hacer un 'assert' para chekear q:
        # - en todas las realizac B, haya la misma cantidad
        #   particulas.
        # - para cada ID de pla, haya la misma cantidad de 
        #   realizaciones B.
        fname_out = self.dir_dst+'/'+self.fname_out_base
        fo = h5(fname_out,'w')
        for iB in range(nB):
            flist = glob(self.dir_src+'B%02d_traj*.dat'%iB)
            for fname_traj in flist:
                t,x,y,z,mu,err = self.get_traj(fname_traj)

            


#EOF
