set grid
set key box right top
set xlabel "RIGIDITY [GV]"
set ylabel "K_perp / K_parall"
set log x
set log y

FILE='sigma_1.0_coeficientes_difusion_vs_energia.dat'
pl FILE u 1:($2/$3) pt 7 ps 1.4 ti "K_perp / K_parall\n sigma=1.0"
