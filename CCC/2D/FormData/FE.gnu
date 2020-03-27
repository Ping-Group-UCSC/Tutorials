#!/usr/bin/env gnuplot

# general settings
set term postscript eps color "Helvetica,20"

GAP = 4.7055

set ylabel '{/Helvetica-Bold FE (eV)}' offset -2,0,0
set xlabel '{/Helvetica-Bold {/Symbol e}_F {/Symbol - }VBM (eV)}'


set yrange [-1217.4:-1213.9]
# set ytics 1
# set mxtics 2
unset ytics

set xrange [0:GAP]
set xtics 1
set mytics 2

set parametric
set trange [-10:10]

# set key horizontal center outside top
set key at 0.3,-1217.4+0.2 bottom left
set border lw 3
# set grid lw 1.5 dt 2 lc rgb '#666666'
set size 0.7,1


####################################################################################################
set output "FE.eps"
p \
'edit_fe_Ti.dat' w lp ls 7 lw 4 ps 2 t 'Ti'
