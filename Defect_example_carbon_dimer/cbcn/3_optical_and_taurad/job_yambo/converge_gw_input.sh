#!/bin/bash 
bands='600 750 900 1050 1200 1350 1500 1750 1900 2000'
blocks='1 2 3 4 5 6 7 8 9 10'
for i in ${bands}
 do
  for j in ${blocks}
   do
    sed  -e "s/1 | 2000/1 | $i/g"  gw_ff.in  > tmp$i
    sed  -e "s/NGsBlkXp= 1/ NGsBlkXp= $j/g"  tmp$i > gw_conv_$i'b_'$j'Ry.in'
    rm tmp*
   done
 done
