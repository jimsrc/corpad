#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import h5py, os, argparse, glob
from pylab import figure, show, close

#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-fpattern', '--fpattern',
type=str,
default='./test.h5',
help='HDF5 input for R-function profiles',
)
pa = parser.parse_args()

da = {}
flist = glob.glob(pa.fpattern)
assert len(flist)>0, '\n [-] no input files found!\n'
for fnm in flist:
    _Lc = float(fnm.split('_')[-1].replace('.h5','').replace('Lc.',''))
    os.system('ls -lh %s'%fnm); assert os.path.isfile(fnm)
    fi = h5py.File(fnm, 'r')
    psim = {}
    #for nm, var in fi['psim'].iteritems():
    #    psim[nm] = var
    psim = dict([(nm,var.value) for nm,var in fi['psim'].iteritems()])
    Nm = psim['Nm_slab']
    lc = fi['fit/lcorr'].value # fitted value
    print Nm, lc
    da[_Lc] = lc

fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

x = np.log10([da[nm] for nm in da.keys()]) # correlation lengths (fitted)
y = np.log10(da.keys()) # correlation scales
# sort according to 'x'
cc = x.argsort()        # indexes that sort this array
x, y = x[cc], y[cc]

# plot data
ax.plot(x, y, 'o')

# plot fit
m, b = np.polyfit(x, y,deg=1,cov=False)
ax.plot(x, m*x+b, '--', c='red', label='$y = %2.4g x + %2.4g$'%(m,b))
#ax.plot(x, 10.**(m*np.log10(x)+b), '--', c='red', label='fit')

ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel('$log_{10}(\lambda_c)$')      # correlation length (fitted)
ax.set_ylabel('$log_{10}(L_c)$')            # correlation scale
ax.set_title('$Lc$: correlation scale\n$\lambda_c$: correlation length')

fnmfig = './linfit.png'
fig.savefig(fnmfig, dpi=135, bbox_inches='tight')
close(fig)

"""
Using:
Nm   = 129
slab = 0.2
etal,
we find:
y = 1.035*x - 0.4266
where:
y = log10(lambda_c) # log of the REAL correlation length
x = log10(Lc)  # correlation scale of the model

NOTE: we find this in the range of Lc: [0.01, 3.0]
NOTE: similarly, we find that the relation: x = 0.97*y + 0.41
"""

#EOF
