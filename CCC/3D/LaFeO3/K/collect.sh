#!/usr/bin/env bash

conv=27.2113961318

fil0=Q0/main.out
fil1=Q-1/main.out
fil2=Q-2/main.out
filp=../Prist/main.out

e0=`grep Etot $fil0 | tail -1 | awk -F Etot: '{print $2}' | awk '{print $1}'`
e1=`grep Etot $fil1 | tail -1 | awk -F Etot: '{print $2}' | awk '{print $1}'`
e2=`grep Etot $fil2 | tail -1 | awk -F Etot: '{print $2}' | awk '{print $1}'`
ep=`grep Etot $filp | tail -1 | awk -F Etot: '{print $2}' | awk '{print $1}'`

c1=`grep Net $fil1 | awk '{print $3}'`
c2=`grep Net $fil2 | awk '{print $3}'`

vbm=`grep HOMO: $filp | awk '{print $2}'`

t1=`echo $conv $e1 $e0 $c1 $vbm | awk '{print $1*($2-$3+$4-$5)}'`
t2=`echo $conv $e2 $e1 $c2 $c1 $vbm | awk '{print $1*($2-$3+$4-$5-$6)}'`
t3=`echo $conv $e2 $e0 $c2 $vbm | awk '{print $1*($2-$3+$4-2*$5)}'`

#echo $e0 $e1 $e2 $ep $c1 $c2 $vbm

echo "ctl(0/-1) ctl(-1/-2) | ctl(0/-2)"
echo "$t1    $t2    |  $t3"
