from pylab import *
from numpy import *

def gamma(Ek):
	Eo = 938272013.0
        return Ek/Eo + 1.

def calc_kperp(Lc, sigma, b):
        c  = 3e10               # [cm/s]
        v  = b*c                # [cm/s]
        return 5.*v*Lc/12. * sin(3.*pi/5.)*sigma**2

def calc_beta(Ek):
	Eo= 938272013.0
        g = gamma(Ek)
        b = (1. - g**-2)**.5
        return b

fname	= '../../../post/diff.pars_Nm128_slab0.20_sig.1.0e+00_Lc2d.1.0e-02_LcSlab.1.0e-02.dat'
data	= loadtxt(fname, unpack=True)

EKs	= data[0]
RLs	= data[1]
Kperp	= data[2]
Kparall = data[3]
MFPperp	= data[4]
MFPparall=data[5]

Eks	= array(pow(10., linspace(5, 10, 30)))	# [eV]
beta	= calc_beta(Eks)
AUincm	= 1.5e13		# [cm]
Lc	= 1e-2*AUincm		# [cm]
sig	= 1.0
kperp	= calc_kperp(Lc, sig, beta)

fig	= figure(1, figsize=(4, 4))
ax	= fig.add_subplot(111)

ax.plot(EKs, Kperp/(1e18), '-o', c='red', label='our simulations')
ax.plot(Eks, kperp/(1e18)/2., '-', c='black', label='Quasi Linear Theory')
ax.grid()
ax.set_xlabel('$E_K$ [eV]', fontsize=15.)
ax.set_ylabel('$K_\perp$ [$cm^2/s$]  ($\\times 10^{18}$)', fontsize=15.)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e5, 1e10)
ax.set_ylim(1e-2, 5e3)
ax.legend(loc='best')

fname_fig = '../Kperp_vs_Ek_giacalone99.png'
savefig(fname_fig, dpi=100, format='png', bbox_inches='tight')
close()
