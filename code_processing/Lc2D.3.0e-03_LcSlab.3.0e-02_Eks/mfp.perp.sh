AUincm=1.5e13
Lc=0.03*AUincm

file='/home/jim/simulacion/pla_stochastic/composite_turbulence/theory/WNLT/Duu.2D/mfp.profiles/mfp.perp.compos_slab0.20_sig1.00_Lc2D.3.0e-03_LcSlab.3.0e-02.dat'

set ylabel "lambda_perp [AU]"
set xlabel "Rl / Lc   [1]"
set log x
set log y
set grid

pl '../../post/diff.pars_Nm128_slab.0.20_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' u ($2/Lc):($5/AUincm) w lp pt 7
repl file u 1:2 w lp pt 7 ps .3 lc 3
