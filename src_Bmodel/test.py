#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import Bmodel
#--- parameters
from params import (
    nB, pd, psim, pother, mu, ph, AU_in_cm
)

xyz = np.zeros(3, dtype=np.float32)

m = Bmodel.Bmodel()

pd['xi'] = 1.0
pd['ratio_slab'] = 0.2
pd['sigma_Bo_ratio'] = 0.3
m.set_Bmodel(pdict=pd, nB=nB)
print m.Bxyz(xyz)


"""
Rs = np.arange(0.2, 5., 0.05)
x[1] = np.pi/2. # eliptica

Bm = []
for x[0] in Rs:
    B = bp.return_B(x)
    Bm += [ np.square(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]) ]

Bm = np.array(Bm)
"""
#EOF
