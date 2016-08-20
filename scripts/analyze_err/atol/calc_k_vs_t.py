#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from funcs import *
from funcs import (
    value, count_lines_in_file, 
    sqr_deviations_ii
)
from numpy import (
    nan, ones, zeros, empty, 
    mean, median, nanmean, nanmedian, 
    array, savetxt, log10
)
import os
from h5py import File as h5


def nans(n):
    return nan*ones(n)

def sqr_dsplmts(DATA, time):
    """
    desplazamientos cuadraticos entre un tiempo y el consecutivo
    """
    nfil    = len(time)
    #n  = nfil/every + (nfil%every!=0)
    dx2     = nans(nfil)
    dy2     = nans(nfil)
    dz2     = nans(nfil)
    x2      = nans(nfil)
    y2      = nans(nfil)
    z2      = nans(nfil)
    tt      = nans(nfil)
    dt      = nans(nfil)

    tt[0] = dt[0] = nan
    # data de *todas* las plas
    t       = DATA[:,0]
    x       = DATA[:,1]
    y       = DATA[:,2]
    z       = DATA[:,3]

    #j  = 0
    print " ---> calculando <(x-x')^2>(t),  <(y-y')^2>(t),  <(z-z')^2>(t)..."
    for i in range(1, nfil):    
        cc0     = t==time[i-1]      # tiempo t-dt
        cc1     = t==time[i]        # tiempo t
        # delta de tiempo
        dt[i]   = time[i] - time[i-1]
        # posiciones en t-dt
        x0      = x[cc0]
        y0      = y[cc0]
        z0      = z[cc0]
        # posiciones en t
        x1      = x[cc1]
        y1      = y[cc1]
        z1      = z[cc1]

        dx2[i]  = mean( (x1-x0)**2 )    # promedio (over plas) del dezplazam cuadrat
        dy2[i]  = mean( (y1-y0)**2 )
        dz2[i]  = mean( (z1-z0)**2 )
        
        x2[i]   = mean(x1*x1)             # [AU^2]
        y2[i]   = mean(y1*y1)             # [AU^2]
        z2[i]   = mean(z1*z1)             # [AU^2]
        
        tt[i]   = time[i]               # [seg]

    out         = {}
    out['t']    = tt
    out['dt']   = dt
    out['dx2']  = dx2
    out['dy2']  = dy2
    out['dz2']  = dz2
    out['x2']   = x2
    out['y2']   = y2
    out['z2']   = z2

    out['x']    = x
    out['y']    = y
    out['z']    = z
    out['txyz'] = t

    return out
    #return t, dt, dx2, dy2, dz2

class mfp_vs_t(object):
    def __init__(s, dir_src):
        s.dir_src   = dir_src
        dir_info    = '%s/info' % dir_src
        s.fname_orient = '%s/orientations.in' % dir_info
        s.fname_plas  = '%s/plas.in' % dir_info
        s.fname_turb  = '%s/turb.in' % dir_info
        s.NPLAS       = count_lines_in_file(s.fname_orient)
        #--- simulation parameters
        s.par = {
        'nplas'  : s.NPLAS,
        'nB'     : int(value(s.fname_plas, 'nB')),
        'Nm_s'   : int(value(s.fname_turb, 'Nm_slab')),
        'Nm_2d'  : int(value(s.fname_turb, 'Nm_2d')),
        'sig'    : value(s.fname_turb, 'sigma_Bo_ratio'),
        #'perc_2d': value(s.fname_turb, 'percent_2d'),
        'perc_slab': value(s.fname_turb, 'ratio_slab'),
        'Lc_slab': value(s.fname_turb, 'Lc_slab'),
        'xi'     : value(s.fname_turb, 'xi'),
        'lmin_s' : value(s.fname_turb, 'lmin_s'),
        'lmin_2d': value(s.fname_turb, 'lmin_2d'),
        'lmax_s' : value(s.fname_turb, 'lmax_s'),
        'lmax_2d': value(s.fname_turb, 'lmax_2d'),
        }
        #s.par.update({
        #'Lc_2d' : s.par['Lc_slab']*s.par['xi'],
        #})

    def calc_mfp_profile(s, moreinfo=False):
        """
        build 'self.profile' dict that contains
        statistics on the particles trajectories
        """
        s.fname_inp = s.dir_src+'/out.h5'
        fi   = h5(s.fname_inp,'r')
        assert fi['ntot_pla'].value==s.NPLAS, \
            ' > problem with self.NPLAS!'
        print " > calculating mfp ..."
        SQR  = sqr_deviations_ii(s.fname_inp, moreinfo)
        # tiempos 
        t       = SQR['time'] #SQR['t_dim']     # [seg]
        # promedios sobre realizaciones
        x2_avr  = SQR['x2_mean']
        y2_avr  = SQR['y2_mean']
        z2_avr  = SQR['z2_mean']
        # desv standard correspondientes
        x2_std  = SQR['x2_std']
        y2_std  = SQR['y2_std']
        z2_std  = SQR['z2_std']
        # *otra* desv standard correspondientes
        #x2_std2 = SQR['x2_std2']*AUincm**2  # [cm2]
        #y2_std2 = SQR['y2_std2']*AUincm**2  # [cm2]
        #z2_std2 = SQR['z2_std2']*AUincm**2  # [cm2]

        nt       = t.size 
        kxx      = x2_avr/(2.*t)    # [1]
        kyy      = y2_avr/(2.*t)    # [1]
        kzz      = z2_avr/(2.*t)    # [1]
        kxx_std  = x2_std/(2.*t)    # [1]
        kyy_std  = y2_std/(2.*t)    # [1]
        kzz_std  = z2_std/(2.*t)    # [1]

        Lc_slab = s.par['Lc_slab']  # [1]
        s.profile = {
        'time'   : t,               # [1]
        'lxx'    : 3.*kxx/Lc_slab,     # [1]
        'lyy'    : 3.*kyy/Lc_slab,     # [1]
        'lzz'    : 3.*kzz/Lc_slab,     # [1]
        'lxx_std': 3.*kxx_std/Lc_slab, # [1]
        'lyy_std': 3.*kyy_std/Lc_slab, # [1]
        'lzz_std': 3.*kzz_std/Lc_slab, # [1]
        #'kxx_std2': kxx_std2,      # [cm2/s]
        #'kyy_std2': kyy_std2,      # [cm2/s]
        #'kzz_std2': kzz_std2,      # [cm2/s]
        'nB'     : SQR['nB'],
        'npla'   : SQR['npla'],
        }
        if moreinfo:
            s.profile.update({
            'x2': SQR['x2'],
            'y2': SQR['y2'],
            'z2': SQR['z2'],
            })

    def save2file(s, dir_out, label):
        fname_out = dir_out+'/'+label+'.h5'
        print " ---> saving: "+fname_out+'\n'
        fo = h5(fname_out, 'w')
        #--- time profiles
        for name in s.profile.keys():
            fo[name] = s.profile[name]
        #--- simulation parameters
        for pname, pval in s.par.iteritems():
            fo['psim/'+pname] = pval
        fo.close()

        print " > placing a link to output-file in: "+s.dir_src
        os.system('ln -sf {fname_out} {symlink}'.format(
            fname_out=fname_out,
            symlink=s.dir_src+'/post.h5')
        )

#++++++++++++++++++++++++++++++++++++++++++++++++++
"""
if __name__=='__main__':
"""

##
