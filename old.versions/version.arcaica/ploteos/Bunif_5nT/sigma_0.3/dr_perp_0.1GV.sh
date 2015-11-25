set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_0.1GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.1 GV"'
#----------------------------------------------
usin='u 1:($2*1e10)'
sty='pt 7'

#-----------------------ajuste lineal para los datos
f(x) = m*x + b
fit f(x) FILE @usin via m, b
omega=0.478589
pendiente_fisica=m*omega
pendiente(m) = sprintf("pendiente= %3.1e cm^2/seg", pendiente_fisica)

#-----------------------ploteo
pl FILE @usin @sty @title
repl f(x) lt -1 ti pendiente(m)
