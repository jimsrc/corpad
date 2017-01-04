#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import plot, show, close, figure, grid
from pylab import xlabel, ylabel, xscale, yscale, legend
import funcs as ff


k1, dV1, sgiac = ff.spectra2D_GiacaloneAndJokipii(5e-3, 1e2, 128)
sgiac_norm = sgiac*dV1 / (sgiac*dV1).sum()
#plot(k1, sgiac_norm, label='G&J')
plot(k1, sgiac, label='G&J')

k2, dV2, sshal = ff.spectra2D_Shalchi(5e-3, 1e2, 128)
sshal_norm = sshal*dV2 / (sshal*dV2).sum()
#plot(k2, sshal_norm, label='Shalchi')
plot(k2, sshal, label='Shalchi')


grid(True)
xscale('log'); yscale('log')
legend(loc='best')
show(); close()


#EOF
