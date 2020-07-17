#!/usr/bin/env gnuplot

# general settings
set term postscript eps color "Helvetica,20"

cbm = 2.22
polaron = 1.872886331

set label 1 'free polaron' at polaron+0.1,1 rotate by 270 tc rgb '#222222'

set ylabel '{/Helvetica-Bold FE (eV)}' offset 0,0,0
set xlabel '{/Helvetica-Bold {/Symbol e}_F {/Symbol - }VBM (eV)}'

set xrange [0:cbm]
set yrange [-0.8:4.8]

set ytics 1
set mxtics 2
set xtics 1
set mytics 2

set parametric
set trange [-5:5]

# first plot
set output "ef1.eps"
p \
polaron,t w l lw 3 dt 2 lc rgb '#666666' t '', \
'vo.dat' w lp ls 7 lw 4 ps 2 lc rgb 'red' t 'V_{O}'

# second plot
set term postscript eps color "Helvetica,40"
set tics font ',36'
set border lw 3
set output "ef2.eps"
set label 2 at 0.8,2.5 "+2" tc rgb '#000000'
set label 3 at 1.3,3.3 "+1" tc rgb '#000000'
set label 4 at 2.0,3.6  "0" tc rgb '#000000'
set size  1,2
set key at 0.3,4.6 left
replot
