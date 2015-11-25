from pylab import *
import numpy as np

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def beta(Ek, Eo):
    g = gamma(Ek, Eo)
    b = (1. - g**-2)**.5
    return b

def gamma(Ek, Eo):
    return Ek/Eo + 1.

def kperp(Lc, sigma, b):
    c  = 3e10       # [cm/s]
    v  = b*c        # [cm/s]
    return 5.*v*Lc/12. * sin(3.*pi/5.)*sigma**2

# omega ciclotron para proton:
# Ek: [eV] energ cinetica
# Bo: [G] campo guia
def calc_omega(Bo, Ek):
        q = 0.00000000048032
        mo = 1.6726e-24
        c = 3e10
        E_reposo = 938272013.0

        rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo);
        gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
        beta = pow(1. - 1/(gamma*gamma) , 0.5);
        OMEGA   = q * Bo / (gamma * mo * c)             # [s^-1]
        return OMEGA

def larmor(Ek):
    Eo = 938272013.0    # [eV]
    b = beta(Ek, Eo)
    c  = 3e10           # [cm/s]
    v = b*c
    wc = calc_omega(5e-5, Ek)   # [s-1]
    rl = v/wc           # [cm]
    return rl

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



fname_fig = './nu.o.omega_vs_Ek.png'
fname_inp = './tau_vs_Ek.txt'
Ek, tau, dum1, dum2 = np.loadtxt(fname_inp, skiprows=True).T

#print Ek, nu

fig     = figure(1, figsize=(6,4))
ax      = fig.add_subplot(111)


ax.plot(Ek, 1./(tau), '-o', color='black', markersize=10)
ax.grid()

ax.set_xlabel('$E_k$ [eV]', fontsize=23)
ax.set_ylabel('$\\nu_\mu / \Omega$', fontsize=26)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.7*1e6, 1.5*1e9)
ax.tick_params(labelsize=19)


AUinkm = 1.5e8      # [km]
AUincm = 1e5*AUinkm # [cm]
Lc = 1e-2*AUincm    # [cm]
rl = larmor(Ek)/Lc
print rl


#------ el otro eje x
ax2 = ax.twiny()
ax2.plot(Ek, 1./(tau), '-o', c='none', markeredgewidth=0)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.7*1e6, 1.5*1e9)
ax2.set_xlabel('$R_L / L_C$', fontsize=25)
ax2.set_xticks(Ek)
X2LABELS = []# array(['0','1','2', '3'])
for i in range(4):
	X2LABELS += [ "%2.2f"%(rl[i]) ]
X2LABELS = array(X2LABELS)
ax2.set_xticklabels(X2LABELS)

ax2.tick_params(labelsize=19)



savefig(fname_fig, dpi=135, bbox_inches='tight')

print " ----> generamos: " + fname_fig
close()
