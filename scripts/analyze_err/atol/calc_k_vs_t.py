#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from funcs import *
from numpy import (
    nan, ones, zeros, empty, 
    mean, median, nanmean, nanmedian, 
    array, savetxt
)
import os
from h5py import File as h5


def nans(n):
    return nan*ones(n)


"""
desplazamientos cuadraticos entre un tiempo y el consecutivo
"""
def sqr_dsplmts(DATA, time):
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


class diff_mgr:
    def __init__(s, Ek, dir_data, dir_plots, dir_out):
        s.dir_data  = dir_data
        s.dir_plots = dir_plots
        s.dir_out   = dir_out
        s.dir_info  = '%s/info' % dir_data
        s.fname_orient  = '%s/orientations.in' % s.dir_info
        s.fname_plas    = '%s/plas.in' % s.dir_info
        s.fname_turb    = '%s/turb.in' % s.dir_info
        s.NPLAS         = count_lines_in_file(s.fname_orient)
        # orden del nro max de giroperiodos
        s.order_tmax  = int(log10(value(s.fname_plas, 'nmax_gyroperiods')))
        # nro de filas x archivo (nro de puntos q le pedi a la simulacion)
        s.nfil        = int(s.order_tmax*value(s.fname_plas, 'npoints') + 1)
        s.ncol        = 6         # nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
        s.NBs       = int(value(s.fname_plas, 'nro_Bfield_realizations'))
        s.rigidity  = value(s.fname_plas, 'rigidity')
        #--------------
        s.Nm        = int(value(s.fname_turb, 'n_modos'))
        s.Bo        = value(s.fname_turb, 'Bunif')
        s.Sig       = value(s.fname_turb, 'sigma_Bo_ratio')
        s.perc_2d   = value(s.fname_turb, 'percent_2d')
        s.perc_slab = value(s.fname_turb, 'percent_slab')
        s.Lc_2d     = value(s.fname_turb, 'Lc_2d')
        s.Lc_slab   = value(s.fname_turb, 'Lc_slab')
        s.lambda_min = value(s.fname_turb, 'lambda_min')
        s.lambda_max = value(s.fname_turb, 'lambda_max')


    def k_vs_t(s):
        #---------------------
        # nok   : nro de files q existe Y tienen data
        # nbad  : nro de files q solicite y no existen
        # time  : grilla temporal
        DATA, time, nok, nbad = load_trajectories(s.NBs, s.NPLAS, s.nfil, s.ncol, s.dir_data)

        print " nro de plas: ", s.NPLAS
        print " nro de B-realizations: ", s.NBs
        print " nro de ptos por trayectoria: %d\n" % s.nfil
        print " nro de archivos q existe c/data: %d/%d " % (nok, nok+nbad)
        print " nro de archivos q pedi y NO existen: %d/%d " % (nbad, nok+nbad)
        #---------------------
        s.DATA  = DATA
        s.time  = time

        s.out   = sqr_dsplmts(DATA, time)


class mfp_vs_t:
    def __init__(s, dir_src):
        s.dir_src   = dir_src
        dir_info    = '%s/info' % dir_src
        s.fname_orient = '%s/orientations.in' % dir_info
        s.fname_plas  = '%s/plas.in' % dir_info
        s.fname_turb  = '%s/turb.in' % dir_info
        s.NPLAS       = count_lines_in_file(s.fname_orient)

        # simulation parameters
        s.par = {
        'nplas'  : s.NPLAS,
        'nB'     : int(value(s.fname_plas, 'nB')),
        #'rigidity': value(s.fname_plas, 'rigidity'),
        'Nm_s'   : int(value(s.fname_turb, 'Nm_slab')),
        'Nm_2d'  : int(value(s.fname_turb, 'Nm_2d')),
        #'Bo': value(s.fname_turb, 'Bunif'),
        'sig'    : value(s.fname_turb, 'sigma_Bo_ratio'),
        #'perc_2d': value(s.fname_turb, 'percent_2d'),
        'perc_slab': value(s.fname_turb, 'ratio_slab'),
        'Lc_slab': value(s.fname_turb, 'Lc_slab'),
        #'xi'     : value(s.fname_turb, 'xi'),
        #'Lc_slab': value(s.fname_turb, 'Lc_slab'),
        'lmin_s' : value(s.fname_turb, 'lmin_s'),
        'lmin_2d': value(s.fname_turb, 'lmin_2d'),
        'lmax_s' : value(s.fname_turb, 'lmax_s'),
        'lmax_2d': value(s.fname_turb, 'lmax_2d'),
        }
        #s.par.update({
        #'Lc_2d' : s.par['Lc_slab']*s.par['xi'],
        #})

    def calc_mfp_profile(s, dir_out, label, moreinfo=False):
        fname_out = dir_out+'/'+label+'.h5'
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

        print " ---> saving: "+fname_out+'\n'
        fo = h5(fname_out, 'w')
        for name in s.profile.keys():
            fo[name] = s.profile[name]
        fo.close()
        print " > placing a link to output-file in: "+s.dir_src
        os.system('ln -sf {fname_out} {symlink}'.format(
            fname_out=fname_out,
            symlink=s.dir_src+'/post.h5')
        )



def generate_k_vs_t(Ek, dir_data):
    dir_info    = '%s/info' % dir_data
    fname_orient    = '%s/orientations.in' % dir_info
    fname_plas  = '%s/plas.in' % dir_info
    fname_turb  = '%s/turb.in' % dir_info
    NPLAS       = count_lines_in_file(fname_orient)
    # orden del nro max de giroperiodos
    order_tmax  = int(log10(value(fname_plas, 'nmax_gyroperiods')))
    # nro de filas x archivo (nro de puntos q le pedi a la simulacion)
    nfil        = int(order_tmax*value(fname_plas, 'npoints') + 1)
    ncol        = 6         # nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
    NBs     = int(value(fname_plas, 'nro_Bfield_realizations'))
    rigidity    = value(fname_plas, 'rigidity')
    #--------------
    Nm      = int(value(fname_turb, 'n_modos'))
    Bo      = value(fname_turb, 'Bunif')
    Sig     = value(fname_turb, 'sigma_Bo_ratio')
    perc_2d     = value(fname_turb, 'percent_2d')
    perc_slab   = value(fname_turb, 'percent_slab')
    Lc_2d       = value(fname_turb, 'Lc_2d')
    Lc_slab     = value(fname_turb, 'Lc_slab')
    lambda_min  = value(fname_turb, 'lambda_min')
    lambda_max  = value(fname_turb, 'lambda_max')

    print " ------> Ek [eV]: %g" % Ek
    calc_k_versus_t(dir_data, Ek, Sig, NPLAS, NBs, nfil, ncol, Bo, 
            Lc_2d, Lc_slab, Nm, perc_slab)


def calc_k_versus_t(dir_data, Ek, Sig, NPLAS, NBs, nfil, ncol, Bo,
        Lc_2d, Lc_slab, Nm, perc_slab):
    dir_plots = '../../plots'
    #dir_data= '../../output/Ek.%1.1eeV/sig%d' % (Ek, Sig)
    """dir_out    = '../../post/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
    try: os.system('mkdir %s' % dir_out)
    except: print ""
    """
    dir_out   = '../../post'
    fname_out = '%s/k_vs_t_Ek.%1.1eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
    #---------------------
    # nok   : nro de files q existe Y tienen data
    # nbad  : nro de files q solicite y no existen
    # time  : grilla temporal
    DATA, time, nok, nbad = load_trajectories(NBs, NPLAS, nfil, ncol, dir_data)

    print " nro de plas: ", NPLAS
    print " nro de B-realizations: ", NBs
    print " nro de ptos por trayectoria: %d\n" % nfil
    print " nro de archivos q existe c/data: %d/%d " % (nok, nok+nbad)
    print " nro de archivos q pedi y NO existen: %d/%d " % (nbad, nok+nbad)
    #---------------------
    every   = 1         # no en c/tiempo, sino cada 'every'
    tt, x2, y2, z2 = sqr_deviations(DATA, time, every)

    AUinkm  = 1.5e8
    AUincm  = AUinkm*1e5        # [cm]
    r2  = x2 + y2
    r2  = r2*AUincm**2      # [cm^2]
    x2  = x2*AUincm**2      # [cm^2]
    y2  = y2*AUincm**2      # [cm^2]
    z2  = z2*AUincm**2      # [cm^2]
    wc  = calc_omega(Bo, Ek) #4.781066E-01 #4.325188E-01 #Ek=1e8eV  #4.735689E-01 # Ek=1e7eV #4.781066E-01 # Ek=1e6eV
    print " wc[s-1]: ", wc
    tt  = tt*wc             # [1]
    #-------------------
    kxx = x2/(2.*tt/wc)     # [cm2/s]
    kyy = y2/(2.*tt/wc)     # [cm2/s]
    kzz = z2/(2.*tt/wc)     # [cm2/s]
    #-- guarda data kxx(t)
    data_out    = array([tt, kxx, kyy, kzz]).T
    data_out    = data_out[1:]              # el 1er tiempo no lo guardo xq es division por zero 1/2t
    print " ---> guardando: %s" % fname_out
    print ""
    savetxt(fname_out, data_out, fmt='%12.2f')


#++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__=='__main__':
    Ek          = 6e5   # [eV]
    Nm          = 128
    perc_slab   = 0.01 #0.02 #0.05 #0.10 # 0.00, 0.20, 0.40, 0.60, 1.00
    sig         = 1e0
    Lc_2d       = 3e-3
    Lc_slab     = 3e-2
    dir_data    = '../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)
    dir_plots   = '../../plots'
    dir_out     = '../../post'
    dm          = diff_mgr(Ek, dir_data, dir_plots, dir_out)
    print "---> leamos trayectorias..."
    dm.k_vs_t()


##
