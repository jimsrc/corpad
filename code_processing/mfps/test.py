from pylab import *
from numpy import *
import os

def calc_beta_gamma(Ek):
	E_reposo = 938272013.0		# [eV]
	rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo)
	gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
	beta = pow(1. - 1/(gamma*gamma) , 0.5)
	return (beta, gamma)

dir_plots = './'
AUincm	= 1.5e13	# [cm]
c 	= 3e10		# [cm/s]
Ek	= 1e7		# [eV]
fnameNm = []
dataNm	= []
kx=[]; ky=[]; kz=[]; kperp=[]
for i in range(0, 5):
	fnameNm  += [ '../../post/Nm%d/k_vs_t_Ek.1e%deV.dat' % (i, int(log10(Ek))) ]
	SIZE_FILE= os.stat(fnameNm[i]).st_size
	if SIZE_FILE==0:
		dataNm+=[0]; kx+=[0]; ky+=[0]
		kperp+=[0]; kz+=[0]
	else:
		dataNm 	 += [ array(loadtxt(fnameNm[i])).T ]
		t	 = dataNm[i][0]
		kx	 += [ dataNm[i][1] ]
		ky	 += [ dataNm[i][2] ] 
		kperp	 += [ .5*(kx[i]+ky[i]) ]
		kz	 += [ dataNm[i][3] ]

#----------------------------------------------
Nm = []
Nm += [[ 460, 'black']]
Nm += [[ 128, 'green']]
Nm += [[ 64, 'red']]
Nm += [[ 16, 'blue']]
Nm += [[ 4, 'magenta']]
fig	= figure(1, figsize=(10, 6))
ax	= fig.add_subplot(111)
for i in range(5):
	try:
		LABEL='Nm: %d' % Nm[i][0]
		ax.plot(t, kperp[i], 'o-', c=Nm[i][1], alpha=.3, label=LABEL)
	except:
		continue

ax.legend(loc='upper left')
ax.grid()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time [1/$\Omega$]')
ax.set_ylabel('$K_{\perp}$(t) [cm2/s]')

fname_fig='%s/kperp_allNms_Ek.1e%deV.png' % (dir_plots, int(log10(Ek)))
savefig(fname_fig, dpi=135, format='png', bbox_inches='tight')
close()

#----------------------------------------------
fig	= figure(1, figsize=(6, 4))
ax	= fig.add_subplot(111)
for i in range(5):
	try:
		LABEL='Nm: %d' % Nm[i][0]
		beta, g = calc_beta_gamma(Ek)
		mfp_perp = 3.*kperp[i]/(beta*c)
		ax.plot(t, mfp_perp/AUincm, 'o-', c=Nm[i][1], alpha=.3, label=LABEL)
	except:
		continue

ax.legend(loc='lower center')
ax.grid()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(15e-5, 1e-2)
ax.set_xlabel('time [1/$\Omega$]', fontsize=17)
ax.set_ylabel('$\lambda_{\perp}$(t) [AU]', fontsize=17)

fname_fig='%s/mfp.perp_allNms_Ek.1e%deV.png' % (dir_plots, int(log10(Ek)))
savefig(fname_fig, dpi=135, format='png', bbox_inches='tight')
close()
