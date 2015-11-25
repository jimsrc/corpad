reset; clear
set grid
set key box right top
set xlabel "RIGIDITY [GV]"
set ylabel "K_perp / K_parall"
set log x
set log y

FILE_1='sigma_0.3_coeficientes_difusion_vs_energia.dat'
FILE_2='sigma_1.0_coeficientes_difusion_vs_energia.dat'
pl FILE_1 u ($2-1):($3/$4) pt 7 ps 1.4 ti "sigma=0.3"
repl FILE_2 u ($2-1):($3/$4) pt 7 ps 1.4 lt 3 ti "sigma=1.0"

