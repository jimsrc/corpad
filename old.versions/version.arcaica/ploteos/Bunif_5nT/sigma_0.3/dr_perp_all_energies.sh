set grid
set key right bottom
FILE_1='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
FILE_2='../plas_0.1GV/avrg_displacements/sqr_displacements.dat'
FILE_3='../plas_1.0GV/avrg_displacements/sqr_displacements.dat'
FILE_4='../plas_10.0GV/avrg_displacements/sqr_displacements.dat'

set xlabel "time [1/omega]"
set ylabel "< dr_perp^2 > [cm^2]"
set log y
set macros
ti_1='ti "0.01 GV"'
ti_2='ti "0.1 GV"'
ti_3='ti "1 GV"'
ti_4='ti "10 GV"'
usin='u 1:($2*1e10)'

pl FILE_1 @usin @ti_1
repl FILE_2 @usin @ti_2
repl FILE_3 @usin @ti_3
repl FILE_4 @usin @ti_4
