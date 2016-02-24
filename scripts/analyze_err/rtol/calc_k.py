#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from funcs import *
"""
Puedo comparar estos resultdos con los de Qin. Ver figura 3.4 de Shalchi 2009 (libro).
Ahi muestra k(t). Son perfiles muy parecidos.
"""
class kdiff:
    def __init__(self, Ek, dir_info, dir_post, dir_plots):
        self.dir_info   = dir_info
        self.dir_plots  = dir_plots
        self.dir_post   = dir_post
        self.Ek         = Ek

        fname_plas      = '%s/plas.in' % dir_info
        fname_turb      = '%s/turb.in' % dir_info
        fname_orient    = '%s/orientations.in' % dir_info
        NPLAS           = count_lines_in_file(fname_orient)
        order_tmax      = int(log10(value(fname_plas, 'nmax_gyroperiods')))
        nfil            = int(order_tmax*value(fname_plas, 'npoints') + 1)
        ncol            = 6

        self.par  = par  = {}
        par_names = ('n_modos', 'Bunif', 'sigma_Bo_ratio', 'percent_2d', 'percent_slab', \
                    'Lc_2d', 'Lc_slab', 'lambda_min', 'lambda_max')
        for pname in par_names:
            par[pname] = value(fname_turb, pname)
            if pname=='n_modos':
                par[pname] = int(par[pname]) # special treatment :\

        fname_inp    = '%s/k_vs_t_Ek.%1.1eeV' % (dir_post, Ek)
        fname_inp   += '_Nm%03d' % par['n_modos']
        fname_inp   += '_slab%1.2f' % par['percent_slab']
        fname_inp   += '_sig.%1.1e' % par['sigma_Bo_ratio']
        fname_inp   += '_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (par['Lc_2d'], par['Lc_slab'])
        self.fname_inp = fname_inp

        #t, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
        t_dim, t_adim, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
        self.inp = {
            't_dim' : t_dim,
            't_adim': t_adim,
            'kxx'   : kxx,
            'kyy'   : kyy,
            'kzz'   : kzz
        }

        self.wc    = calc_omega(par['Bunif'], Ek)
        self.k = {}

    def fit_kperp(self, t_decr):
        t = self.inp['t_dim']
        kxx = self.inp['kxx']
        kyy = self.inp['kyy']
        wc  = self.wc
        #-------------------------------------------- fit k_perp
        #cc = t>5000
        cc  = t>t_decr
        b   = 1.8e19/wc
        m   = 3.6e22/wc
        xo  = 0.
        px  = make_fit_kxx([(t/wc)[cc], kxx[cc]], [b, m, xo])
        py  = make_fit_kxx([(t/wc)[cc], kyy[cc]], [b, m, xo])

        self.k['kxx_fit']    = px[0]
        self.k['kxx_fit_p1'] = px[1]
        self.k['kxx_fit_p2'] = px[2]

        self.k['kyy_fit']    = py[0]
        self.k['kyy_fit_p1'] = py[1]
        self.k['kyy_fit_p2'] = py[2]
        
        self.t_decr_kperp = t_decr
        self.kperp = 0.5*(self.k['kxx_fit'] + self.k['kyy_fit'])


    def plot_kperp(self):
        # we need this parameter to make a
        # reasonable plot
        assert hasattr(self, 't_decr_kperp'), \
            "\n ----> antes de plotear, hay q fitear kperp!!\n"

        t    = self.inp['t_dim']
        cc   = t>self.t_decr_kperp

        kxx  = self.inp['kxx']
        kyy  = self.inp['kyy']
        k    = self.k
        wc   = self.wc
        par  = self.par

        # build fitted curves
        fitted_kxx = fun_hyperbola(k['kxx_fit'], k['kxx_fit_p1'], k['kxx_fit_p2'], t[cc]/wc)
        fitted_kyy = fun_hyperbola(k['kyy_fit'], k['kyy_fit_p1'], k['kyy_fit_p2'], t[cc]/wc)

        fig1 = figure(1, figsize=(6, 4))
        ax1  = fig1.add_subplot(111)
        # plot *all* simulation data
        ax1.scatter(t, kxx, edgecolor='none', c='red', alpha=.3)
        ax1.scatter(t, kyy, edgecolor='none', c='blue', alpha=.3)
        # plot *only* fitted data
        ax1.scatter(t[cc], kxx[cc], edgecolor='none', c='red', label='kxx')
        ax1.scatter(t[cc], kyy[cc], edgecolor='none', c='blue', label='kyy')
        # plot fitted curves
        ax1.plot(t[cc], fitted_kxx, c='red')
        ax1.plot(t[cc], fitted_kyy, c='blue')


        FIT_RESULTS = 'fit: y=b+m/(x-xo)' +\
        '\nKxx: %1.1e   Kyy:%1.1e'%(k['kxx_fit'],k['kyy_fit']) +\
        '\nKperp: %1.2e' % self.kperp

        SIMULAT_PARAMS  = 'Ek[eV]:%1.1e    perc_slab:%1.2f' % (self.Ek, par['percent_slab']) +\
        '\nNm:%d    $(\sigma/Bo)^2$:%1.1e' % (par['n_modos'], par['sigma_Bo_ratio']) +\
        '\n$Lc^{2D}$[AU]:%1.1e    $Lc^{slab}$[AU]:%1.1e' % (par['Lc_2d'], par['Lc_slab'])

        TITLE   = '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)

        ax1.set_title(TITLE)
        ax1.set_xlabel('time [$\Omega^{-1}$]')
        ax1.set_ylabel('[cm2/s]')
        ax1.grid()
        #--------------------------------------------
        fname_fig_perp = '%s/kperp_asymptotic.fit' % self.dir_plots +\
        '_Ek.%1.2eeV' % self.Ek +\
        '_Nm%03d' % par['n_modos'] +\
        '_slab%1.2f' % par['percent_slab'] +\
        '_sig.%1.1e' % par['sigma_Bo_ratio'] +\
        '_Lc2d.%1.1e' % par['Lc_2d'] +\
        '_LcSlab.%1.1e.png' % par['Lc_slab']

        print " ---> generando: %s" % fname_fig_perp
        fig1.savefig(fname_fig_perp, format='png', dpi=135, bbox_inches='tight')
        close(fig1)


    def fit_kparall(self, t_decr):
        t   = self.inp['t_dim']
        kzz = self.inp['kzz']
        wc  = self.wc
        #----------------------------------- fit k_parall
        #cc = t>5000
        cc  = t>t_decr
        b   = 3.5e22*(.2/wc)
        m   = -7.4e24*(.2/wc)
        xo  = 0.
        pz  = make_fit_kxx([(t/wc)[cc], kzz[cc]], [b, m, xo])

        self.k['kzz_fit']    = pz[0]
        self.k['kzz_fit_p1'] = pz[1]
        self.k['kzz_fit_p2'] = pz[2]

        self.t_decr_kparall = t_decr
        self.kparall = self.k['kzz_fit']


    def plot_kparall(self):
        # we need this parameter to make a
        # reasonable plot
        assert hasattr(self, 't_decr_kparall'), \
            "\n ----> antes de plotear, hay q fitear kparall!!\n"

        t    = self.inp['t_dim']
        cc   = t>self.t_decr_kparall

        kzz  = self.inp['kzz']
        k    = self.k
        wc   = self.wc
        par  = self.par

        # build fitted curves
        fitted_kzz = fun_hyperbola(k['kzz_fit'], k['kzz_fit_p1'], k['kzz_fit_p2'], t[cc]/wc)

        fig1 = figure(1, figsize=(6, 4))
        ax1  = fig1.add_subplot(111)
        # plot *all* simulation data
        ax1.scatter(t, kzz, edgecolor='none', c='black', alpha=.3)
        # plot *only* fitted data
        ax1.scatter(t[cc], kzz[cc], edgecolor='none', c='black', label='kzz')
        # plot fitted curves
        ax1.plot(t[cc], fitted_kzz, c='green', lw=3, alpha=.6)

        FIT_RESULTS = 'fit: y=b+m/(x-xo)' +\
        '   Kparall: %1.1e' % k['kzz_fit']

        SIMULAT_PARAMS  = 'Ek[eV]:%1.1e    perc_slab:%1.2f' % (self.Ek, par['percent_slab']) +\
        '\nNm:%d    $(\sigma/Bo)^2$:%1.1e' % (par['n_modos'], par['sigma_Bo_ratio']) +\
        '\n$Lc^{2D}$[AU]:%1.1e    $Lc^{slab}$[AU]:%1.1e' % (par['Lc_2d'], par['Lc_slab'])

        TITLE   = '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)

        ax1.set_title(TITLE)
        ax1.set_xlabel('time [$\Omega^{-1}$]')
        ax1.set_ylabel('[cm2/s]')
        ax1.grid()
        #--------------------------------------------
        fname_fig_zz = '%s/kzz_asymptotic.fit' % self.dir_plots +\
        '_Ek.%1.2eeV' % self.Ek +\
        '_Nm%03d' % par['n_modos'] +\
        '_slab%1.2f' % par['percent_slab'] +\
        '_sig.%1.1e' % par['sigma_Bo_ratio'] +\
        '_Lc2d.%1.1e' % par['Lc_2d'] +\
        '_LcSlab.%1.1e.png' % par['Lc_slab']

        print " ---> generando: %s" % fname_fig_zz
        fig1.savefig(fname_fig_zz, format='png', dpi=135, bbox_inches='tight')
        close(fig1)


#+++++ run example +++++
if __name__=='__main__':
    import os
    PLAS    = os.environ['PLAS']
    Nm      = 128
    sig     = 1e0
    Lc_2d   = 3e-3
    Lc_slab = 3e-2
    Ek      = 6e5       # [eV]
    perc_slab   = 0.2

    dir_post     = '../../../post'
    dir_plot    = '.'
    dir_info    = '%s/output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e/info' % (PLAS, Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)

    kd = kdiff(Ek, dir_info, dir_post, dir_plot)

    kd.fit_kperp(t_decr=800.) #TDECRs[i]  #200#200 #500 # 1000
    kd.plot_kperp()

    kd.fit_kparall(t_decr=0.)
    kd.plot_kparall()

#EOF
