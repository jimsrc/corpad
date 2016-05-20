#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import numpy as np

rep, ep = np.loadtxt('./e+.png.dat')[:,:2].T
rem, em = np.loadtxt('./e-.png.dat')[:,:2].T

d   = np.zeros(rep.size) # has less points than 'tem'
r   = np.zeros(rep.size) # has less points than 'tem'
em_ = np.zeros(rep.size) # has less points than 'tem'
ep_ = np.zeros(rep.size) # has less points than 'tem'
#--- sync e+ && e-
for i in range(rep.size):
    diff = np.min(np.abs(rep[i] - rem))
    cc   = diff==np.abs(rep[i] - rem)
    if rem[cc].size>1:
        print rem[cc]
        r[i]   = rep[i]
        em_[i] = em[cc].mean() # agarro el promedio de los "degenerados"
        ep_[i] = ep[i]

    elif rem[cc].size==1:
        r[i]   = rep[i]
        em_[i] = em[cc]
        ep_[i] = ep[i]

    else:
        print "wierd..."
        break

#--- sync e+ && rA
r_  = np.nan*np.ones(rep.size)
rA_ = np.nan*np.ones(rep.size)
rd, rA = np.loadtxt('rA.png.dat').T
for i in range(rep.size):
    diff = np.min(np.abs(rep[i] - rd))
    cc   = diff==np.abs(rep[i] - rd)
    if rd[cc].size>1:
        print rd[cc]
        r_[i]   = rep[i]
        rA_[i] = rA[cc].mean() # agarro el promedio de los "degenerados"
        ep_[i] = ep[i]

    elif rd[cc].size==1:
        r_[i]   = rep[i]
        rA_[i] = rA[cc]
        ep_[i] = ep[i]

    else:
        print "wierd..."
        break

rA_[:2] = np.nan # se confunde el algoritmo. Para estos valores de r no hay data (Helios)

#--- deducimos el e^b   :D :D
eb = (ep_+em_)/(1.0+rA_)

if __name__=='__main__':
    #--- figs
    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    ax2 = ax.twinx()

    #--- rA
    ax2.plot(rd, rA, '-v', color='grey', ms=1, lw=0.5)
    ax2.plot(r_, rA_, '-v', color='grey', ms=1, lw=2., label='$r_A$')
    ax2.legend(loc='best')

    #--- e+ && e-
    ax.plot(r, ep_, '-o', alpha=0.5, label='$e+$')
    ax.plot(r, em_, '-s', alpha=0.5, label='$e-$')
    ax.grid()
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.3, 5.0); ax.set_ylim(1e2, 2.5e4)

    #--- eb
    ax.plot(r, eb, '-*', label='$e_b$')

    ax.legend(loc='lower left')
    savefig('./ep_em_vs_helioradius.png', dpi=300, bbox_inches='tight')
    close()

#EOF
