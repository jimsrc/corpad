set grid
set key box right bottom

set xlabel "RIDIGITY [GV]"
set ylabel "[cm^2 / s]"
set log x
set log y
FILE='coeficientes_difusion_vs_energia.dat'
pl FILE u 1:2 pt 7 ps 1.3 ti "K_perp"
repl FILE u 1:3 pt 9 ps 1.3 lt 3 ti "K_parall"
