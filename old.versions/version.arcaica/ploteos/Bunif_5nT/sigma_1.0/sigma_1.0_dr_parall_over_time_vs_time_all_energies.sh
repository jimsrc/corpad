reset; clear
set grid
set key box right bottom
set xlabel "t [1/omega]"
set ylabel "<dr_parall^2> / t"

#---------------------------------para c/sigma y valores de normalizacion en [cm²/s]
SIGMA='"sigma=1.0"'
FILE_1='../../../plas_0.01GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_2='../../../plas_0.1GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_3='../../../plas_1.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_4='../../../plas_10.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
usin_1='u 1:($3*1e10/$1 / 0.9e20)'
usin_2='u 1:($3*1e10/$1 / 2.6e21)'
usin_3='u 1:($3*1e10/$1 / 1.1e23)'
usin_4='u 1:($3*1e10/$1 / 7.4e24)'
title_1='ti "0.01GV (x0.9e19)"'
title_2='ti "0.1GV (x2.6e21)"'
title_3='ti "1GV (x1.1e23)"'
title_4='ti "10GV (x7.4e24)"'

set label @SIGMA at 200, 0.1 font "arial, 13"
pl FILE_1 @usin_1 w lp @title_1
repl FILE_2 @usin_2 w lp @title_2
repl FILE_3 @usin_3 w lp @title_3
repl FILE_4 @usin_4 w lp @title_4
