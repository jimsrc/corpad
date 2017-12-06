#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from datetime import datetime, timedelta
from h5py import File as h5
import os, sys
#---- graphics
from matplotlib import cm
from matplotlib.colors import LogNorm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
if 'DISPLAY' in os.environ:
    from pylab import figure, savefig, close


def LinearFit_Hist2D(x, y, H2D):
    """
    Linear fit to 2d-histogram.
    NOTE: It's assumed an square 2d-histogram """
    assert x.size==y.size, 'ERROR: 2d-histogram must be square!'
    bins    = x.size
    xvec    = np.empty(bins*bins)
    yvec    = np.empty(bins*bins)
    # construyo grilla
    for i in range(bins):
        ini = i*bins
        end = (i+1)*bins
        yvec[ini:end] = y
        xvec[ini:end] = x[i]*np.ones(bins)

    H = np.reshape(H2D, (bins*bins))
    wei = np.sqrt(H)+1e-15      # factores de peso
    p = np.polyfit(xvec, yvec, 1, w=wei, cov=True)
    m = p[0][0]
    b = p[0][1]
    return {'m':m, 'b':b}

def contour_2d(fig, ax, x, y, mat, hscale='log', **kargs):
    cb_label = kargs.get('cb_label', 'points per bin square')
    cb_fontsize = kargs.get('cb_fontsize', 15)

    # max/min values of the colorbar
    cbmin, cbmax = kargs.get('vmin',1), kargs.get('vmax',1e3)
    if cbmin is None: cbmin = np.nanmin(mat)
    if cbmax is None: cbmax = np.nanmax(mat)

    opt = {
    'linewidth': 0,
    'cmap': cm.gray_r,                # gray-scale
    'vmin': cbmin, #kargs.get('cbmin',1),
    'vmax': cbmax, #kargs.get('cbmax',1000),
    'alpha': kargs.get('alpha',0.9),
    'levels': np.arange(cbmin, cbmax, (cbmax-cbmin)/kargs.get('nlevels',10))
    }
    if hscale=='log':
        opt.update({'norm': LogNorm(),})
    #--- 2d contour
    surf = ax.contourf(x, y, mat, facecolors=cm.jet(mat), **opt)
    sm = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)
    sm.set_array(mat)
    #--- colorbar
    axcb = fig.colorbar(sm)
    axcb.set_label(cb_label, fontsize=cb_fontsize)
    sm.set_clim(vmin=cbmin, vmax=cbmax)
    return fig, ax

#EOF
