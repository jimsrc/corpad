#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plots sigmaBoB :D
"""
from pylab import *
import numpy as np
from eval_r import r_eb, sigmaBoB, eb, Va

if __name__=='__main__':
    #--- figs
    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    #ax2 = ax.twinx()

    #--- eb
    #ax.plot(r_eb, sigmaBoB, '-*', label='$e_b/V_A^2$')
    ax.plot(r_eb, eb, '-*', label='$e_b$')
    ax.plot(r_eb, Va*Va, '-s', c='red', label='$V_A^2$')

    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.3, 5.)
    ax.grid()

    ax.legend(loc='lower left')
    savefig('./eb_and_Va2.png', dpi=300, bbox_inches='tight')
    close()
#EOF
