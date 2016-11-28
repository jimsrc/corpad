#
from pylab import *
from numpy import *
from lmfit import minimize, Parameters, Parameter, report_errors
from h5py import File as h5
from numpy import (
    zeros, empty, ones, nan, 
    savetxt, loadtxt, mean, median, std,
    nanmean, nanmedian, log10
)
from os.path import isfile, isdir
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import ticker # handles ticks
from numpy import power as pow
M_PI = np.pi
#

class mfp_mgr(object):
    """
    - read the processed .h5 file 'fname_inp'
    - append calculated parameters into the file
    - generate figures of mfp(t)
    """
    def __init__(self, dir_fig='plots', fname_inp='None'):
        """
        fname_inp: full path of the post.h5 file
        dir_fig: full path to directory for the .pdf output
        """
        #self.prefix  = prefix
        #self.myid    = myid
        #--- directory for figures
        self.dir_fig = dir_fig
        #--- fullpath of input file
        self.fname_inp = fname_inp

        assert isfile(self.fname_inp),\
        ' ---> doesn\'t exist: '+self.fname_inp+\
        '\n Aborting...'
        self.f = h5(self.fname_inp,'r')
        self.psim = {} # dict of sim-parameters
        for nm, v in self.f['psim'].iteritems():
            self.psim[nm] = v.value

        self.t  = self.f['time'][...]
        self.nB = self.f['nB'].value
        self.nP = self.f['npla'].value

    def fits_and_plots(self, t_decr, b_pe, m_pe, b_pa, m_pa):
        basename = self.fname_inp.split('/')[-1].replace('.h5','.pdf')
        # .pdf adopts the name of fname_inp's last inner subdir
        fname_out_pdf = self.dir_fig+'/'+basename
        pdf_pages = PdfPages(fname_out_pdf)
        #--- 1st page
        self.fit_perp(t_decr, b_pe, m_pe)
        fig, ax = self.plot_perp()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 2nd page
        self.fit_parall(t_decr, b_pa, m_pa)
        fig, ax = self.plot_parall()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 3rd page
        fig, ax = self.plot_TauHist_ii()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)
        """
        #--- 4th page
        fig, ax = self.plot_TauHist(scale='lmin')
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)
        
        #--- 5th page
        fig, ax = self.plot_TauHist(scale='lmax')
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)
        """
        #--- 6th page
        fig, ax = self.plot_ThetaHist()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- 7th page
        fig, ax = self.plot_RuntimeHist()
        pdf_pages.savefig(fig, bbox_inches='tight')
        close(fig)

        #--- Write the PDF document to the disk
        pdf_pages.close()
        print " ---> we generated: " + fname_out_pdf 

    def plot_TauHist(self, scale='omega'):
        hx, hc = self.f['hist_tau'][...] # global version
        fig = figure(1, figsize=(6, 4))
        ax  = fig.add_subplot(111)
        if scale=='omega':
            hx_ = np.power(10.,hx-log10(2.*M_PI))
            ax.set_xlabel('$log_{10}(\Omega \\tau_{back}/2\pi)$')
        elif scale=='lmin':
            hx_ = np.power(10.,hx-log10(self.f['psim/lmin_s'].value))
            ax.set_xlabel('$log_{10}(v \\tau_{back}/\lambda_{min})$')
        elif scale=='lmax':
            hx_ = np.power(10.,hx-log10(self.f['psim/lmax_s'].value))
            ax.set_xlabel('$log_{10}(v \\tau_{back}/\lambda_{max})$')

        ax.plot(hx_, hc, 'k-o', ms=4)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('#')
        ax.grid(True)
        return fig, ax

    def _return_twiny(self, ax, func12, func21, xlabel, offset=-0.2):
        """
        Make a twiny() x-axis, so that we can put it BELOW the
        host x-axis `ax`.
        For the conversion of units, we use `func12` and `func21`.
        The position offset of the twin x-axis is `offset`.
        This is for self.plot_TauHist_ii().
        """
        ax2 = ax.twiny()
        #---------------------
        # Move twinned axis ticks and label from top to bottom
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")

        # Offset the twin axis below the host
        ax2.spines["bottom"].set_position(("axes", offset))

        # Turn on the frame for the twin axis, but then hide all 
        # but the bottom spine
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        for sp in ax2.spines.itervalues():
            sp.set_visible(False)
        ax2.spines["bottom"].set_visible(True)
        #---------------------
        ax2.set_xlabel(xlabel)
        #---------------------
        xmin, xmax = ax.get_xlim()
        xticks = ax.get_xticks()

        x2min, x2max = func21(xmin), func21(xmax)

        new_tick_labels = pow(10., np.arange(
            start = np.floor(log10(x2min)), 
            stop  = np.ceil(log10(x2max)),
            step  = 1.,
        ))
        ax2.set_xscale('log')
        new_tick_locations = func12(new_tick_labels)
        ax2.set_xticks(new_tick_locations[1:])
        flog = lambda x: '$10^{%d}$'% log10(x)
        ax2.set_xticklabels([flog(t) for t in new_tick_labels[1:]])
        # remove minor x-ticks [inherited from the host]
        ax2.xaxis.set_minor_locator(ticker.NullLocator())
        #ax2.xaxis.set_minor_locator(ticker.LogLocator())
        #-- we'll manually set the minor ticks
        mt = []
        for ti, te in zip(new_tick_labels[:-1], new_tick_labels[1:]):
            mt += list(np.linspace(ti, te, 10)) 
        ax2.xaxis.set_ticks(
        ticks = [func12(_mt) for _mt in mt], 
        minor = True,
        )
        ax2.set_xlim(xmin, xmax)
        return ax2

    def plot_TauHist_ii(self,):
        hx, hc = self.f['hist_tau'][...] # global version
        fig = figure(1, figsize=(6, 4))
        ax  = fig.add_subplot(111)

        hx_ = np.power(10.,hx-log10(2.*M_PI))
        ax.set_xlabel('$\Omega \\tau_{back}/2\pi$')

        ax.plot(hx_, hc, 'k-o', ms=4)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('#')
        ax.set_ylim(1.,)
        ax.grid(True)

        #-- units of \lambda_min
        func21 = lambda x: (2.*M_PI/self.f['psim/lmin_s'].value)*x
        func12 = lambda x: (self.f['psim/lmin_s'].value/(2.*M_PI))*x
        ax2 = self._return_twiny(
            ax, func12, func21, 
            xlabel = '$v \\tau_{back}/\lambda_{min}$',
            offset = -0.2
        )

        #--- units of \lambda_max
        func21 = lambda x: (2.*M_PI/self.f['psim/lmax_s'].value)*x
        func12 = lambda x: (self.f['psim/lmax_s'].value/(2.*M_PI))*x
        ax3 = self._return_twiny(
            ax, func12, func21, 
            xlabel = '$v \\tau_{back}/\lambda_{max}$',
            offset = -0.4,
        )

        #--- units of parallel mean-free-path
        RloLc = 1./self.psim['Lc_slab']
        func21 = lambda x: RloLc*(2.*M_PI*x)/(self.lparall)
        func12 = lambda x: self.lparall*x/(RloLc*2.*M_PI)
        ax4 = self._return_twiny(
            ax, func12, func21,
            xlabel = '$v \\tau_{back}/\lambda_{\parallel}$',
            offset = -0.6,
        )


        return fig, ax

    def plot_ThetaHist(self):
        hx, hc = self.f['hist_theta'][...] # global version
        fig = figure(1,figsize=(6,4))
        ax  = fig.add_subplot(111)
        ax.plot(hx, hc, 'k-o', ms=4)
        ax.set_xlabel('$\\theta_{back}$ [deg]')
        ax.set_ylabel('#')
        ax.grid(True)
        return fig, ax

    def plot_RuntimeHist(self):
        """
        histogram of all the runtimes of all particles
        associated to all the B-realizations
        """
        hx, hc = self.f['hist_runtime'][...] 
        total_trun = (hx*hc).sum()/(60.*24.) # [days]
        label='total/days:%2.1f'%total_trun
        dx = hx[1]-hx[0]
        fig = figure(1,figsize=(6,4))
        ax  = fig.add_subplot(111)
        ax.bar(hx, hc, width=dx, align='center', label=label)
        ax.set_xlabel('runtime [min]')
        ax.set_ylabel('#')
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend(loc='best')
        return fig, ax

    def fit_perp(self, t_decr, seed_b, seed_m):
        t   = self.t
        lxx = self.f['lxx'][...]
        lyy = self.f['lyy'][...]
        # select data for fitting
        cc = t>t_decr
        
        seed_xo  = 0.0
        px  = make_fit_lxx(
                [t[cc], lxx[cc]], 
                [seed_b, seed_m, seed_xo]
              )
        py  = make_fit_lxx(
                [t[cc], lyy[cc]], 
                [seed_b, seed_m, seed_xo]
              )
        self.l = {
        'lxx_fit'    : px[0],
        'lxx_fit_p1' : px[1],
        'lxx_fit_p2' : px[2],
        'lyy_fit'    : py[0],
        'lyy_fit_p1' : py[1],
        'lyy_fit_p2' : py[2],
        }
        self.t_decr_kperp = t_decr
        self.lperp = 0.5*(self.l['lxx_fit'] + self.l['lyy_fit'])
        print self.lperp

    def plot_perp(self):
        assert hasattr(self, 't_decr_kperp'), \
            "\n ----> antes de plotear, hay q fitear kperp!!\n"

        # colors for kxx & kyy
        ckxx, ckyy = 'red', 'blue'

        t    = self.t
        cc   = t>self.t_decr_kperp

        lxx  = self.f['lxx'][...]
        lyy  = self.f['lyy'][...]
        lxx_std = self.f['lxx_std']
        lyy_std = self.f['lyy_std']
        l    = self.l
        nB   = self.nB
        # build fitted curves
        fitted_lxx = fun_hyperbola(l['lxx_fit'], l['lxx_fit_p1'], l['lxx_fit_p2'], t[cc])
        fitted_lyy = fun_hyperbola(l['lyy_fit'], l['lyy_fit_p1'], l['lyy_fit_p2'], t[cc])

        fig1 = figure(1, figsize=(6, 4))
        ax1  = fig1.add_subplot(111)
        # plot *all* simulation data
        opt = {'c': 'none', 'alpha':.3, 's':9}
        ax1.scatter(t[~cc], lxx[~cc], edgecolor='red', **opt)
        ax1.scatter(t[~cc], lyy[~cc], edgecolor='blue', **opt)
        # plot *only* fitted data
        opt = {'edgecolor': 'none', 'alpha': 0.4, 's': 9}
        ax1.scatter(t[cc], lxx[cc], c=ckxx, label='lxx', **opt)
        ax1.scatter(t[cc], lyy[cc], c=ckyy, label='lyy', **opt)
        # plot fitted curves
        ax1.plot(t[cc], fitted_lxx, c=ckxx)
        ax1.plot(t[cc], fitted_lyy, c=ckyy)

        # errores kxx
        inf = lxx - lxx_std/np.sqrt(nB)
        sup = lxx + lxx_std/np.sqrt(nB)
        ax1.fill_between(t[cc], inf[cc], sup[cc], facecolor=ckxx, alpha=0.5)

        # errores kyy
        inf = lyy - lyy_std/np.sqrt(nB)
        sup = lyy + lyy_std/np.sqrt(nB)
        ax1.fill_between(t[cc], inf[cc], sup[cc], facecolor=ckyy, alpha=0.5)

        # legend
        ax1.legend()

        FIT_RESULTS = 'fit: y=b+m/(x-xo)' +\
        '\nlxx: %1.1e   lyy:%1.1e'%(l['lxx_fit'],l['lyy_fit']) +\
        '\nlperp: %1.2e' % self.lperp
        RloLc = 1.0/self.psim['Lc_slab']  # [1]
        SIMULAT_PARAMS  = '$R_L/Lc$:{RloLc:1.1e}    perc_slab:{perc_slab:1.2f}'.format(RloLc=RloLc, **self.psim) +\
        '\nNmS, Nm2d: {Nm_s}/{Nm_2d}    $(\sigma/Bo)^2$:{sig:1.1e}'.format(**self.psim) +\
        '\n$L_c^{{2D}}/L_c^{{Slab}}$:{xi:1.1e}'.format(**self.psim)

        TITLE   = '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)

        ax1.set_title(TITLE)
        ax1.set_xlabel('$\Omega t$')
        ax1.set_ylabel('$\lambda_{{perp}}/L_C^{slab}$')
        ax1.grid()

        return fig1, ax1


    def fit_parall(self, t_decr, seed_b, seed_m):
        t   = self.t
        lzz = self.f['lzz'][...]
        # select data for fitting
        cc = t>t_decr
        #----------------------------------- fit k_parall
        seed_xo  = 0.0
        pz  = make_fit_lxx([t[cc], lzz[cc]], [seed_b, seed_m, seed_xo])

        self.l['lzz_fit']    = pz[0]
        self.l['lzz_fit_p1'] = pz[1]
        self.l['lzz_fit_p2'] = pz[2]

        self.t_decr_kparall = t_decr
        self.lparall = self.l['lzz_fit']

    def plot_parall(self):
        assert hasattr(self, 't_decr_kparall'), \
            "\n ----> antes de plotear, hay q fitear kperp!!\n"

        t    = self.t
        cc   = t>self.t_decr_kparall

        lzz  = self.f['lzz'][...]
        lzz_std = self.f['lzz_std'][...]
        l    = self.l
        nB   = self.nB

        # build fitted curves
        fitted_lzz = fun_hyperbola(l['lzz_fit'], l['lzz_fit_p1'], l['lzz_fit_p2'], t[cc])

        fig1 = figure(1, figsize=(6, 4))
        ax1  = fig1.add_subplot(111)
        # plot *all* simulation data
        ax1.scatter(t, lzz, edgecolor='none', c='black', alpha=.3)
        # plot *only* fitted data
        ax1.scatter(t[cc], lzz[cc], edgecolor='none', c='black', label='kzz')
        # plot fitted curves
        ax1.plot(t[cc], fitted_lzz, c='green', lw=3, alpha=.6)

        # errores lzz
        inf = lzz - lzz_std/np.sqrt(nB)
        sup = lzz + lzz_std/np.sqrt(nB)
        ax1.fill_between(t[cc], inf[cc], sup[cc], facecolor='gray', alpha=0.5)
        # legend
        ax1.legend()

        FIT_RESULTS = 'fit: y=b+m/(x-xo)' +\
        '\nlzz: %1.1e'%(l['lzz_fit']) +\
        '\nlparall: %1.2e' % self.lparall
        RloLc = 1.0/self.psim['Lc_slab']  # [1]
        SIMULAT_PARAMS  = '$R_L/Lc$:{RloLc:1.1e}    perc_slab:{perc_slab:1.2f}'.format(RloLc=RloLc, **self.psim) +\
        '\nNmS, Nm2d: {Nm_s}/{Nm_2d}    $(\sigma/Bo)^2$:{sig:1.1e}'.format(**self.psim) +\
        '\n$L_c^{{2D}}/L_c^{{Slab}}$:{xi:1.1e}'.format(**self.psim)

        TITLE   = '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)

        ax1.set_title(TITLE)
        ax1.set_xlabel('$\Omega t$')
        ax1.set_ylabel('$\lambda_{{\parallel}}/L_C^{slab}$')
        ax1.grid()

        return fig1, ax1



def value(fname, value_name):
    cc = read_contents(fname)
    n = len(cc)
    assert n>0, " > NO CONTENTS!!"
    for i in range(n):
        if value_name==cc[i][1]:
            return float(cc[i][0])

    print " --> WRONG NAME (%s) IN"%value_name+\
        " INPUT FILE (%s)!!"%fname
    raise SystemExit


def read_contents(fname):
    f = open(fname, 'r')
    content = []
    for line in f:
        ll = line.split()
        content += [ll]

    return content

def count_lines_in_file(fname):
    file    = open(fname, 'r')
    n   = 0
    for line in file:
        n += 1

    return n

def make_fit_lxx(data, sems):
    x       = data[0]
    y       = data[1]
    # create a set of Parameters
    params = Parameters()
    params.add('b')
    params.add('m')
    params.add('xo')

    SEM_b   = sems[0]
    SEM_m   = sems[1]
    SEM_xo  = sems[2]

    params['b'].value       = SEM_b
    params['b'].vary        = True
    """params['b'].min              = 
    params['b'].max         = """

    params['m'].value       = SEM_m
    params['m'].vary        = True
    """params['m'].min              = 
    params['m'].max         = """

    params['xo'].value       = SEM_xo
    params['xo'].vary        = True
    """params['xo'].min              = 
    params['xo'].max         = """

    METHOD  = "leastsq"#"leastsq"#"lbfgsb"
    result = minimize(residuals_kxx, params, args=(x, y), method=METHOD)
    #print " ----> now dir "
    #print dir(result)

    # write error report
    print " --------> METODO_FITEO: %s" % METHOD
    print " --------> funcion: %s" % 'EXPRESION GIACALONE-99'
    #report_errors(params)

    par     = zeros(3)
    #par[0]  = result.values['b']
    #par[1]  = result.values['m']
    #par[2]  = result.values['xo']
    par[0]  = result.params['b'].value
    par[1]  = result.params['m'].value
    par[2]  = result.params['xo'].value
    return par

def model(N, K, L, t):
    sum = 0.
    for n in range(1, N+1):
        expon = -t * K/(L**2) * ((2.*n-1)*pi/2.)**2
        sum   += (-1.)**(n+1) / (2.*n-1) *exp(expon)

    sum *= 4./pi
    return sum

def make_fit(data, sems):
    x       = data[0]
    y       = data[1]
    # create a set of Parameters
    params = Parameters()
    params.add('N')
    params.add('A')
    params.add('yo')
    params.add('K')
    params.add('L')

    SEM_N   = sems[0]
    SEM_A   = sems[1]
    SEM_yo  = sems[2]
    SEM_K   = sems[3]
    SEM_L   = sems[4]

    params['N'].value       = SEM_N
    params['N'].vary        = False
    """params['A'].min              = 
    params['A'].max         = """

    params['A'].value       = SEM_A
    params['A'].vary        = True
    """params['A'].min         = -.2
    params['A'].max         = +.2"""

    params['yo'].value       = SEM_yo
    params['yo'].vary        = True
    """params['yo'].min         = -.2
    params['yo'].max         = +.2"""

    params['K'].value      = SEM_K
    params['K'].vary       = True
    """params['mu'].min     = 
    params['mu'].max        = """

    params['L'].value     = SEM_L
    params['L'].vary      = False
    """params['sig'].min    = 
    params['sig'].max       ="""

    METHOD  = "leastsq"#"leastsq"#"lbfgsb"
    result = minimize(residuals, params, args=(x, y), method=METHOD)

    # write error report
    print " --------> METODO_FITEO: %s" % METHOD
    print " --------> funcion: %s" % 'EXPRESION GIACALONE-99'
    #report_errors(params)

    par     = zeros(5)
    par[0]  = result.values['N']
    par[1]  = result.values['A']
    par[2]  = result.values['yo']
    par[3]  = result.values['K']
    par[4]  = result.values['L']
    return par

def model_shifted(N, A, yo, K, L, x):
    return yo + A*model(N, K, L, x)

def residuals(params, x, y_data):
    N       = params['N'].value
    A       = params['A'].value
    yo      = params['yo'].value
    K       = params['K'].value
    L       = params['L'].value
    diff    = (model_shifted(N, A, yo, K, L, x)  - y_data)**2.
    #print " diff---> %f" % mean(diff)

    return diff


def fun_hyperbola(b, m, xo, x):
    fun = b + m/(x-xo)
    return fun


def residuals_kxx(params, x, y_data):
    b   = params['b'].value
    m   = params['m'].value
    xo  = params['xo'].value
    diff    = (fun_hyperbola(b, m, xo, x)  - y_data)**2.
    #print " diff---> %f" % mean(diff)

    return diff


def nans(n):
    return nan*ones(n)


def load_trajectories_from_h5(NBs, NPLAS, dir_data):
    data        = []
    nexist      = 0
    nwdata      = 0
    n_noexist   = 0

    fname_inp = '%s/out.h5' % dir_data
    f =h5(fname_inp, 'r')

    for j in range(NBs):
        for i in range(NPLAS):
            path = 'B%02d/pla%03d' % (j, i)
            try:
                xyz = f[path+'/xyz'][...].T
                t   = f[path+'/t'][...]
                data_aux = nans((t.size, 4))
                data_aux[:, 1:] = xyz
                data_aux[:, 0]  = t
                nexist     += 1
                if len(data_aux)>0:
                    nwdata += 1
                    data    += [data_aux]
                    print " ---> SI EXISTE: %s" % path

            except KeyboardInterrupt:
                print " ----> KEYBOARD.... Exiting..."
                raise SystemExit

            except KeyError:
                print " ---> NO EXISTE: %s" % path
                n_noexist += 1

    nt   = t.size
    DATA = zeros((nwdata*nt, 4))
    j    = 0
    for i in range(nwdata):
        j           = i*nt
        DATA[j:j+nt, :]   = data[i] # eje x de 'DATA' es tiempo

    time    = data[0][:,0]        # [seg] con tomar una muestra es suficiente!

    # realz == realizaciones
    # nwdata    : nro de realz q existe Y tienen data
    # n_noexist : nro de realz q solicite y NO existen
    #return DATA, time, nwdata, n_noexist
    return {'DATA': DATA, 
            'time': time, 
            'nwdata': nwdata, 
            'n_noexist': n_noexist,
            'data': data
            }


def load_trajectories(NBs, NPLAS, nfil, ncol, dir_data):
    data        = []
    nexist      = 0
    nwdata      = 0
    n_noexist   = 0
    for j in range(NBs):
        for i in range(NPLAS):
            fname   = '%s/B%02d_pla%03d_traj.dat' % (dir_data, j, i)
            try:
                data_aux = loadtxt(fname, unpack=True)
                nexist  += 1
                if len(data_aux)>0:
                    nwdata += 1
                    data    += [data_aux]
                    print " ---> SI EXISTE: %s" % fname

            except KeyboardInterrupt:
                print " ----> KEYBOARD.... Exiting..."
                raise SystemExit

            except:
                print " ---> NO EXISTE: %s" % fname
                n_noexist += 1
                #print " no existe: %s" % fname
                #pause(300)
        
    DATA    = zeros(nwdata*nfil*ncol).reshape(nwdata*nfil, ncol)
    j   = 0
    for i in range(nwdata):
        j           = i*nfil
        DATA[j:j+nfil][:]   = data[i].T # eje x de 'DATA' es tiempo

    time    = data[0][0,:]        # [seg] con tomar una muestra es suficiente!

    # nwdata    : nro de files q existe Y tienen data
    # n_noexist : nro de files q solicite y NO existen
    return DATA, time, nwdata, n_noexist

def sqr_deviations(DATA, time, every):
    nfil    = len(time)
    n   = nfil/every + (nfil%every!=0)
    x2  = zeros(n)
    y2  = zeros(n)
    z2  = zeros(n)
    tt  = zeros(n)
    j   = 0
    print " calculando <x2>(t), <y2(t)>..."
    for i in range(0, nfil, every):     # yendo de 'every' en 'every'
        cond    = DATA.T[0]==time[i]
        x       = DATA.T[1][cond]
        y       = DATA.T[2][cond]
        z       = DATA.T[3][cond]
        x2[j]   = mean(x*x)      # [AU^2]
        y2[j]   = mean(y*y)          # [AU^2]
        z2[j]   = mean(z*z)          # [AU^2]
        tt[j]   = time[i]        # [seg]
        j   += 1

    return tt, x2, y2, z2


#def sqr_deviations_ii(DATA, time, every):
def sqr_deviations_ii(fname_inp, moreinfo=False):
    """
    Del archivo de entrada, leeremos trayectorias corresponientes
    a PARTICULAS (de las cuales leemos sus trayectorias) y 
    a REALIZACIONES DE CAMPO. 
    Aqui haremos promedios sobre las particulas en cada realizacion.
    'x2[iB,:]', 'y2[iB,:]', 'z2[iB,:]' son promedios (a lo largo del
    tiempo (axis=1)) sobre particulas, para c/realizacion 'iB'; 
    mientras que 'x2_avr', 'y2_avr', 'z2_avr' son promedios sobre 
    los promedios 'x2[]', 'y2[]', 'z2[]' sobre las realizaciones 
    (axis=0); y 'x2_std', 'y2_std', 'z2_std' son las desviaciones 
    estandard correspondientes.
    """
    f        = h5(fname_inp, 'r')
    nB, npla = f['ntot_B'].value, f['ntot_pla'].value
    nt       = f['ntimes'].value

    x2_avr, y2_avr, z2_avr  = zeros((3,nt))
    x2_std, y2_std, z2_std  = zeros((3,nt))

    x2, y2, z2 = nans((3,nB,nt))
    x2std, y2std, z2std = nans((3,nB,nt))

    print " calculando <x2>(t), <y2(t)>... y sus errores"
    for iB in range(nB):
        path = 'B%02d'%iB
        x = f[path+'/x'][...] # [1] (nt, npla)
        y = f[path+'/y'][...] # [1] (nt, npla)
        z = f[path+'/z'][...] # [1] (nt, npla)
        x2[iB,:] = (x*x).mean(axis=1) # promedio sobre particulas
        y2[iB,:] = (y*y).mean(axis=1) # promedio sobre particulas
        z2[iB,:] = (z*z).mean(axis=1) # promedio sobre particulas
        x2std[iB,:] = (x*x).std(axis=1) # std sobre particulas
        y2std[iB,:] = (y*y).std(axis=1) # std sobre particulas
        z2std[iB,:] = (z*z).std(axis=1) # std sobre particulas
        
    # promedios y errores standard sobre las B-realizaciones
    x2_avr  = x2.mean(axis=0)
    y2_avr  = y2.mean(axis=0)
    z2_avr  = z2.mean(axis=0)
    # error: std sobre los valores medios de c/ B-realizacion
    x2_std  = x2.std(axis=0)
    y2_std  = y2.std(axis=0)
    z2_std  = z2.std(axis=0)
    # otra version de error
    x2_std2 = x2std.mean(axis=0)
    y2_std2 = y2std.mean(axis=0)
    z2_std2 = z2std.mean(axis=0)

    time = f['time'][...] # [1] times of trajectory points
    o = {
    'time'      : time,     # [1]
    'x2_mean'   : x2_avr,   # [1]
    'y2_mean'   : y2_avr,   # [1]
    'z2_mean'   : z2_avr,   # ..
    'x2_std'    : x2_std,   # ..
    'y2_std'    : y2_std,   # ..
    'z2_std'    : z2_std,   # ..
    #'x2_std2'   : x2_std2,  # ..
    #'y2_std2'   : y2_std2,  # ..
    #'z2_std2'   : z2_std2,  # ..
    'nB'        : nB,
    'npla'      : npla,
    }
    if moreinfo:
        o.update({
        'x2'   : x2,    'y2'   : y2,   'z2'   : z2,
        'x2std': x2std, 'y2std':y2std, 'z2std':z2std,
        })
    return o


# omega ciclotron para proton:
# Ek: [eV] energ cinetica
# Bo: [G] campo guia
def calc_omega(Bo, Ek):
    q = 0.00000000048032
    mo = 1.6726e-24
    c = 3e10
    Ereposo = 938272013.0

    rigidity = sqrt(Ek*Ek + 2.*Ek*Ereposo);
    gamma = pow(pow(rigidity/Ereposo,2) + 1. , 0.5)
    beta = pow(1. - 1/(gamma*gamma) , 0.5);
    OMEGA   = q * Bo / (gamma * mo * c)     # [s^-1]
    return OMEGA


def fit_Kperp(data, fparams):
    p = make_fit_kxx(data, fparams)
    fit_b   = p[0]
    fit_m   = p[1]
    fit_xo  = p[2]
    k_fit = fit_b
    return k_fit


# param relativistas de un proton
def calc_beta_gamma(Ek):
    Ereposo = 938272013.0          # [eV]
    rigidity = sqrt(Ek*Ek + 2.*Ek*Ereposo)
    gamma = pow(pow(rigidity/Ereposo,2) + 1. , 0.5)
    beta = pow(1. - 1/(gamma*gamma) , 0.5)
    return (beta, gamma)


# inputs:
# k     : [cm^2/s] difussion coef
# Ek    : [eV] kinetic energy of proton
def calc_mfp(kdiff, Ek):
    c = 3e10            # [cm/s]
    beta, g = calc_beta_gamma(Ek)   # parametros relativ del proton
    mfp = 3.*kdiff / (beta*c)   # [cm]
    return mfp


def larmor(Bo, Ek):
    c  = 3e10                       # [cm/s]
    beta, g = calc_beta_gamma(Ek) 
    v = beta*c          # [cm/s]
    wc = calc_omega(Bo, Ek)         # [s-1]
    rl = v/wc                       # [cm]
    return rl

#
##
