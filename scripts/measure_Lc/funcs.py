#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia


class LcMeasure(object):
    def __init__(self, pd):
        self.m = Bmodel.Bmodel()
        self.m._build_pturb(pd=pd)

    def BiBj_profile(self, ro, dr, ij=(0,0)):
        """
        dada una realizacion de turbulencia, medir los 
        productos Bi(ro)*Bj(ro+dr) para diferentes valores `dr`.
        Output:
        - array of such product.
        """
        n = dr.size
        Bi = np.zeros(n, dtype=np.float32) # Bi(ro)
        Bj = np.zeros(n, dtype=np.float32) # Bj(ro+dr)
        
        Bi[:] = self.m.Bxyz(ro)[ij[0]]
        for i in range(n):
            Bj[i]   = self.m.Bxyz(ro+dr[i])[ij[1]]

        return np.prod([Bi,Bj], axis=0)

    def one_R_realiz(self, Nro=20, dr=None, ij=(0,0), pd=None, nB=0):
        """
        Nro : nro of test-values for `ro`
        dr  : resolution of profile
        ij  : component-coordinates
        pd  : parameters for modeled turbulence
        nB  _ nro of turbulence-realization
        """
        #if 'm' in self.__dict__:
        #    del self.m

        self.m._build_par(nB=nB)

        Lc = pd['Lc_slab']
        BiBj = np.zeros(dr.size, dtype=np.float32)
        for xo, i in zip(np.arange(0., Nro*Lc, 0.5*Lc), range(Nro)):
            ro      = [xo,0.,0.]
            #import pdb; pdb.set_trace()
            BiBj[:] += self.BiBj_profile(ro, dr, ij)

        BiBj /= 1.0*Nro # i want the average (over `ro`)
        return BiBj

#EOF
