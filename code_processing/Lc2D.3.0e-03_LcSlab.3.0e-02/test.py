from pylab import *
import calc_k_vs_t as ck

Ek          = 6e5   # [eV]
Nm          = 128
perc_slab   = 0.01 #0.02 #0.05 #0.10 # 0.00, 0.20, 0.40, 0.60, 1.00
sig         = 1e0
Lc_2d       = 3e-3
Lc_slab     = 3e-2
dir_data    = '../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)
dir_plots   = '../../plots'
dir_out     = '../../post'
dm          = ck.diff_mgr(Ek, dir_data, dir_plots, dir_out)
print "---> leamos trayectorias..."
dm.k_vs_t()

nt      = dm.time.size
t       = dm.DATA[:,0]
x       = dm.DATA[:,1]
y       = dm.DATA[:,2]
z       = dm.DATA[:,3]
dx2     = ck.nans(nt)
dy2     = ck.nans(nt)
dz2     = ck.nans(nt)
x2      = ck.nans(nt)
y2      = ck.nans(nt)
z2      = ck.nans(nt)
dt      = ck.nans(nt)
#cond = dm.DATA[:,0]==dm.time[-1]
tt      = dm.time

for i in range(1,nt):
    print " ---> ", i
    cc0     = t==dm.time[i-1]      # tiempo t-dt
    cc1     = t==dm.time[i]        # tiempo t
    # delta de tiempo
    dt[i]   = dm.time[i] - dm.time[i-1]

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
    x2[i], y2[i], z2[i] = mean(x1**2), mean(y1**2), mean(z1**2)


"""
for i in range(dm.time.size):
    print " ---> ", i
    cond    = t==dm.time[i]
    zz      = z[cond]
    xx      = x[cond]
    yy      = y[cond]
    #hist(xx, bins=30, range=(-1, +1))
    hist(yy, bins=30, range=(-0.05, 0.05), label='t:%g'%dm.time[i])
    legend(loc='best')
    figname = './test/tt_%04d' % i
    savefig(figname, dpi=50)
    close()
"""


