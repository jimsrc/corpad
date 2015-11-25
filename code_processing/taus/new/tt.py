from pylab import *
from numpy import *
"""
def acc_i(v, i):
	sum = 0.0
	for j in range(i+1):
		sum += v[j]
	return sum

def acc(v):
	n = size(v)
	acc_v = zeros(n)
	for i in range(n):
		acc_v[i] = acc_i(v, i)
	return acc_v
"""
#----------------------------------------------------------------
Ek		        = 1e7                   # [eV]
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
exit(1)
pause(10)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dir_plots	= './'
NBs		    = 5 	        # nro de B-realiz
Nplas		= 121 #242	    # nro de plas por B-realiz
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
		cc 	= acc(taus)>2e4
		taus	= taus[cc] 
		Taus	+= [taus]
		ntau	+= size(taus)

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

fname_fig = '%s/hist_tau_1e%deV_Nm%d_half.time.png' % (dir_plots, int(log10(Ek)), Nm)
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
#show()
close()
