AUincm=1.5e13
Lc=0.01*AUincm

set grid
set yrange[0:0.25]
pl 'diff.pars_Nm128_Ek.6.0e+05_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' u 1:($1>0?$5/AUincm:1/0) w lp pt 7
