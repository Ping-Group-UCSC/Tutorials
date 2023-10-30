#!/bin/bash
RL=`grep "Max WF components" ./out/*.save/l_stderr  | awk '{print $7}'`
RL=`echo "scale=2;$RL/(sqrt(45./15.))^3" | bc -l`
RL=`echo "$RL/1+1" | bc`
echo "RL = $RL"
sed "s/XXX/$RL/g" init.in0 > init.in
