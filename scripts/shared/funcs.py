#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
from h5py import File as h5
from os.path import isfile, isdir
import numpy as np
from pylab import (
    pause, find, figure, close
)
#from numpy import min, max
import os, sys
from glob import glob
from Bparker.Bparker import return_B as Bparker_vector
from Bparker.Bparker import calc_Rlarmor
from numpy.linalg import norm

def Lc_memilia(r=1.0):
    """ 
    formula from Maria Emilia's thesis (Sec 4.4, p.88).
    """
    #Lc = 0.89*(r**(0.43))*1e6/(1e8) # [AU]
    Lc = 0.0059*(r**(0.43)) # [AU]
    return Lc


def Bo_parker(r=1.0, th=np.pi/2., ph=0.0):
    """ input:
    r  [AU]  : heliodistance
    th [rad] : spherical polar angle (co-latitude)
    ph [rad] : azimuth
    default: (r=1.0, th=pi/2, ph=0.0)
    """
    #--- calc B-field at position (r, th)
    xyz = np.array([r, th, ph], dtype=np.float32) #[AU], [rad], [rad]: position
    Bo  = norm(Bparker_vector(xyz)) # [Gauss]
    return Bo

#EOF
