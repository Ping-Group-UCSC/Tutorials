#!/bin/bash 
bands='170 180 190 200'
blocks='1 2 3 4 5 6'
for i in ${bands}
do
  for j in ${blocks}
  do
    sed  -e "s/100 | 180/100 | $i/g"  rpa.in  > tmp$i
    sed  -e "s/NGsBlkXd= 1/NGsBlkXd= $j/g"  tmp$i > rpa_conv_100-$i'b_'$j'Ry.in'
    rm tmp*
  done
done
