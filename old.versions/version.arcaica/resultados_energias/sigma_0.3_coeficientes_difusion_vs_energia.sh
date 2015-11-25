set grid
set key box left top

set xlabel "Ek [mo*c^2]"
set ylabel "[cm^2 / s]"
set log x
set log y
FILE='sigma_0.3_coeficientes_difusion_vs_energia.dat'
pl FILE u ($2-1):3 pt 7 ps 1.3 ti "K_perp"
repl FILE u ($2-1):4 pt 9 ps 1.3 lt 3 ti "K_parall\nsigma=0.3"
