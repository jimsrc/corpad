set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_parall^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.01 GV"'
#----------------------------------------------
usin='u 1:($3*1e10)'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_parall^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_0.1GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.1 GV"'
#----------------------------------------------
usin='u 1:($3*1e10)'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_parall^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_10.0GV/avrg_displacements/sqr_displacements.dat'
title='ti "10 GV"'
#----------------------------------------------
usin='u 1:($3*1e10)'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_parall^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_1.0GV/avrg_displacements/sqr_displacements.dat'
title='ti "1 GV"'
#----------------------------------------------
usin='u 1:($3*1e10)'
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
set grid
set key right bottom
FILE_1='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
FILE_2='../plas_0.1GV/avrg_displacements/sqr_displacements.dat'
FILE_3='../plas_1.0GV/avrg_displacements/sqr_displacements.dat'
FILE_4='../plas_10.0GV/avrg_displacements/sqr_displacements.dat'

set xlabel "time [1/omega]"
set ylabel "< dr_parall^2 > [cm^2]"
set log y
set macros
ti_1='ti "0.01 GV"'
ti_2='ti "0.1 GV"'
ti_3='ti "1 GV"'
ti_4='ti "10 GV"'
usin='u 1:($3*1e10)'

pl FILE_1 @usin @ti_1
repl FILE_2 @usin @ti_2
repl FILE_3 @usin @ti_3
repl FILE_4 @usin @ti_4
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.01 GV"'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_10.0GV/avrg_displacements/sqr_displacements.dat'
title='ti "10 GV"'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > [cm^2]"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_1.0GV/avrg_displacements/sqr_displacements.dat'
title='ti "1 GV"'
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
set grid
set macros
set key right bottom
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > / < dr_parall^2 >"
set xrange [0:]
set yrange [0:]
#----------------------para cada energia de pla
FILE='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.01 GV"'
#----------------------------------------------
usin='u 1:($2 / $3)'
sty='pt 7'

#-----------------------ajuste lineal para los datos
#f(x) = m*x + b
#fit f(x) FILE @usin via m, b
omega=0.478589
#pendiente_fisica=m*omega
#pendiente(m) = sprintf("pendiente= %3.1e cm^2/seg", pendiente_fisica)

#-----------------------ploteo
pl FILE @usin @sty @title
#repl f(x) lt -1 ti pendiente(m)
set grid
set macros
set key box right top
set xlabel "time [1/omega] \n omega = 0.478589 s^-1"
set ylabel "< dr_perp^2 > / < dr_parall^2 >"
set xrange [0:]
set yrange [0.001:]
#----------------------para cada energia de pla
FILE_1='../plas_0.01GV/avrg_displacements/sqr_displacements.dat'
FILE_2='../plas_0.1GV/avrg_displacements/sqr_displacements.dat'
FILE_3='../plas_1.0GV/avrg_displacements/sqr_displacements.dat'
FILE_4='../plas_10.0GV/avrg_displacements/sqr_displacements.dat'
title='ti "0.1 GV"'
#----------------------------------------------
usin='u 1:($2 / $3)'
sty='w lp ps 0.8'
set log y

#-----------------------ajuste lineal para los datos
#f(x) = m*x + b
#fit f(x) FILE @usin via m, b
omega=0.478589
#pendiente_fisica=m*omega
#pendiente(m) = sprintf("pendiente= %3.1e cm^2/seg", pendiente_fisica)

#-----------------------ploteo
pl FILE_1 @usin @sty lt 1 ti "0.01 GV"
repl FILE_2 @usin @sty lt 3 ti "0.1 GV"
repl FILE_3 @usin @sty lt 4 ti "1 GV"
repl FILE_4 @usin @sty lt 5 ti "10 GV"
#repl f(x) lt -1 ti pendiente(m)
