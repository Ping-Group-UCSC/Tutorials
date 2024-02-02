#!/usr/bin/gnuplot -persist
set xtics ( "Gamma" 0,  "X" 71,  "W" 107,  "L" 143,  "Gamma" 230,  "K" 322 )
unset key
set terminal png
set output "phbs.png"
nRows = real(system("awk '$1==\"kpoint\" {nRows++} END {print nRows}' bandstruct.kpoints"))
nCols = 6
#set arrow from 71,0. to 71,70 nohead lw 0.5 lt 2 lc rgb 'black'
set ylabel "Phonon energy (meV)" font ",15"
#set yrange [0:]
plot for [i=1:nCols] "dft_phfrq.dat" u 0:(column(i)*27211.386) w l lw 2 lc rgb "red"
