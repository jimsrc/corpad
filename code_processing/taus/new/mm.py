from pylab import *
import numpy as np
from numpy import loadtxt, zeros, ones, array

# freq ciclotron para Bo=5nT
def calc_wc(Ek):
    Bo = (5.0) * 1e-5   # [nT]
    q = 0.00000000048032
    mo = 1.6726e-24
    c = 3e10
    E_reposo=938272013.0 # [eV] proton

    rigidity = np.sqrt(Ek*Ek + 2.0*Ek*E_reposo)

    gamma = np.power(np.power(rigidity/E_reposo,2) + 1. , 0.5)
    OMEGA   = q * Bo / (gamma * mo * c) # [1/s]
    return OMEGA

#----------------------------------------------------------------
Ek		        = 1e6                   # [eV]
Nm		        = 128                   # [1]
slab            = 0.2                   # [1]
sig             = 1.0                   # [1]
Lc_2d, Lc_slab  = 1.0e-2, 1.0e-2        # [AU]
dir_in		    = '../../../output'
dir_in          += '/Ek.%1.1eeV' % Ek
dir_in          += '/Nm%03d' % Nm
dir_in          += '/slab%1.2f' % slab
dir_in          += '/sig.%1.1e' % sig
dir_in          += '/Lc2D.%1.1e_LcSlab.%1.1e' % (Lc_2d, Lc_slab)
#'/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.1.0e-02_LcSlab.1.0e-02' % (int(log10(Ek)), Nm)
print "\n ---> leyendo de: " + dir_in + '\n'
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB, Nplas   = 50, 122
taus        = zeros((NB,Nplas))

for nB in range(NB):
    for npla in range(Nplas):
        fname_inp   = dir_in + '/B%02d_pla%03d_misc.dat' % (nB, npla)
        f           = open(fname_inp, 'r')
        lines       = f.readlines()
        tau_dummy   = lines[1].split()[1] #el <tau> esta en la segunda linea
        taus[nB][npla] = float(tau_dummy)
        f.close()


tau_avr     = np.mean(taus) # es promedio de los-promedios-calculados-por-c++
wc           = calc_wc(Ek) # [1/s]
#tau_avr_phys = tau_avr/4.735689E-01 # [sec]
tau_avr_phys = tau_avr/wc # [sec]
beta         = 1.448440E-01
c            = 3.0e10 #[cm/s]  #300 000 km/s
v            = beta*c
mfp          = v*tau_avr_phys #[cm]
AUincm       = 1.5e13 #[cm]
mfp_         = mfp / AUincm
print " ----> tau*omega: ", tau_avr
print " ----> tau/sec: ", tau_avr_phys
print " ----> mfp/AU: ", mfp_


