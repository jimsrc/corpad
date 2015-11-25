reset; clear
set grid
set key box right bottom
set xlabel "t [1/omega]"
set ylabel "<dr_perp^2> / t"

#---------------------------------para c/sigma y valores de normalizacion en [cm²/s]
SIGMA='"sigma=1.0"'
FILE_1='../../../plas_0.01GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_2='../../../plas_0.1GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_3='../../../plas_1.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_4='../../../plas_10.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
usin_1='u 1:($2*1e10/$1 / 3.0e19)'
usin_2='u 1:($2*1e10/$1 / 4.1e20)'
usin_3='u 1:($2*1e10/$1 / 4.5e21)'
usin_4='u 1:($2*1e10/$1 / 5.9e22)'
title_1='ti "0.01GV (x)"'
title_2='ti "0.1GV (x)"'
title_3='ti "1GV (x)"'
title_4='ti "10GV (x)"'

set label @SIGMA at 200, 0.1 font "arial, 13"
pl FILE_1 @usin_1 w lp @title_1
repl FILE_2 @usin_2 w lp @title_2
repl FILE_3 @usin_3 w lp @title_3
repl FILE_4 @usin_4 w lp @title_4
