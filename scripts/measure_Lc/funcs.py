#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
from src_Bmodel import Bmodel
#from Bparker.Bparker import return_B as Bparker_vector
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia


class LcMeasure(object):
    def __init__(self, pd, nB):
        self.m = Bmodel.Bmodel()
        self.m._build_pturb(pd=pd)
        self.m._build_par(nB=nB)

#EOF
