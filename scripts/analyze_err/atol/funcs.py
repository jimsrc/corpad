#
from pylab import *
from numpy import *
from lmfit import minimize, Parameters, Parameter, report_errors
from h5py import File as h5
from numpy import (
    zeros, empty, ones, nan, 
    savetxt, loadtxt, mean, median, std,
    nanmean, nanmedian
)
#
def value(fname, value_name):
    cc = read_contents(fname)
    n = len(cc)
    for i in range(n):
        if value_name==cc[i][1]:
            return float(cc[i][0])

    return nan

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

def make_fit_kxx(data, sems):
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
    #print dict(result)
    print " ----> now dir "
    print dir(result)

    # write error report
    print " --------> METODO_FITEO: %s" % METHOD
    print " --------> funcion: %s" % 'EXPRESION GIACALONE-99'
    #report_errors(params)

    par     = zeros(3)
    #par[0]  = result.values['b']
    #par[1]  = result.values['m']
    #par[2]  = result.values['xo']
    par[0]  = result.params['b']
    par[1]  = result.params['m']
    par[2]  = result.params['xo']
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
    print " diff---> %f" % mean(diff)

    return diff


def fun_hyperbola(b, m, xo, x):
    fun = b + m/(x-xo)
    return fun


def residuals_kxx(params, x, y_data):
    b   = params['b'].value
    m   = params['m'].value
    xo  = params['xo'].value
    diff    = (fun_hyperbola(b, m, xo, x)  - y_data)**2.
    print " diff---> %f" % mean(diff)

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

    x2_avr, y2_avr, z2_avr  = zeros(nt), zeros(nt), zeros(nt)
    x2_std, y2_std, z2_std  = zeros(nt), zeros(nt), zeros(nt)

    x2, y2, z2 = nans((nB,nt)), nans((nB,nt)), nans((nB,nt))
    x2std, y2std, z2std = nans((nB,nt)), nans((nB,nt)), nans((nB,nt))

    print " calculando <x2>(t), <y2(t)>... y sus errores"
    for iB in range(nB):
        path = 'B%02d'%iB
        x = f[path+'/x'][...] # [AU] (nt, npla)
        y = f[path+'/y'][...] # [AU] (nt, npla)
        z = f[path+'/z'][...] # [AU] (nt, npla)
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

    tsec = f['time'][...] # [seg] times of trajectory points
    o = {
    't_dim'     : tsec,     # [seg]
    'x2_mean'   : x2_avr,   # [AU]
    'y2_mean'   : y2_avr,   # [AU]
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
            'x2': x2, 'y2': y2, 'z2': z2,
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
