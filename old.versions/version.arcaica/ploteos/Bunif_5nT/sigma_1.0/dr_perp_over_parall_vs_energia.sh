reset; clear
set grid
set macros
set key box right top
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > / < dr_parall^2 >"
set xrange [0:]
set yrange [0.001:]
#----------------------para cada energia de pla
FILE_1='../../../plas_0.01GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_2='../../../plas_0.1GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_3='../../../plas_1.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
FILE_4='../../../plas_10.0GV/sigma_1.0/1e3_GYROPERIODS/avrg_displacements/sqr_displacements.dat'
#----------------------------------------------
usin='u 1:($2 / $3)'
sty='w lp ps 0.8'
set log y

#-----------------------ajuste lineal para los datos
#f(x) = m*x + b
#fit f(x) FILE @usin via m, b
omega=0.478589				# VALOR NO-RELATIVISTA!!
#pendiente_fisica=m*omega
#pendiente(m) = sprintf("pendiente= %3.1e cm^2/seg", pendiente_fisica)

#-----------------------ploteo
pl FILE_1 @usin @sty lt 1 ti "0.01 GV"
repl FILE_2 @usin @sty lt 3 ti "0.1 GV"
repl FILE_3 @usin @sty lt 4 ti "1 GV"
repl FILE_4 @usin @sty lt 5 ti "10 GV"
#repl f(x) lt -1 ti pendiente(m)
