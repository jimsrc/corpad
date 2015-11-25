from pylab import *
from numpy import *

Ek		= 1e7
Nm		= 2
dir_in		= '../../output/output_Ek.1e%deV/Nm%d' % (int(log10(Ek)), Nm)
dir_plots	= './'
NBs		= 50 	# nro de B-realiz
Nplas		= 242	# nro de plas por B-realiz
#----------------
Taus = []
ntau = 0
for npla in range(Nplas):
	print " npla: %d/%d,	ntau:%d"% (npla, Nplas-1, ntau)
	for nB in range(NBs):
		fname_in	= '%s/B%02d_pla%03d_misc.dat' % (dir_in, nB, npla)
		taus, dum1, dum2= loadtxt(fname_in, unpack=True, skiprows=2)
		if(min(taus)==0.0):
			print "nB: %d, npla:%d" % (nB, npla)
		Taus		+= [taus]
		ntau		+= size(taus)

TAUS = zeros(ntau)
i = 0
for j in range(NBs*Nplas):
	n	= size(Taus[j])
	TAUS[i:i+n] = Taus[j]
	i	+= n


#----------- histogs
avrTau	= mean(TAUS)
medTau	= median(TAUS)
TITLE	= 'Nm%d, Ek[eV]: %g' % (Nm, Ek)
LABEL 	= 'N: %2.2g\nmean: %g\nmedian: %g' % (ntau, avrTau, medTau)
NBINS 	= 100
h = hist(TAUS, range=(0., 150.), bins=NBINS, normed=True, label=LABEL)
close()
#---------- plots
fig 	= figure(1, figsize=(6, 4))
ax 	= fig.add_subplot(111)

dx 	= h[1][1] - h[1][0]
hcnts   = h[0]*dx*100.
hbins   = .5*(h[1][0:NBINS] + h[1][1:NBINS+1])

ax.plot(hbins, hcnts, '-o', label=LABEL, alpha=.6)
ax.set_yscale('log')
ax.legend(loc='upper right')
ax.grid()
ax.set_xlabel('scattering time [1/$\Omega}$]')
ax.set_ylabel('# [%]')
ax.set_title(TITLE)

fname_fig = '%s/hist_tau_1e%deV_Nm%d.png' % (dir_plots, int(log10(Ek)), Nm)
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
#show()
close()
