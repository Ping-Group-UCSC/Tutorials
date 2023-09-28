#!/bin/bash
DIRY="/sdcc/u/jxuucsc/yambo_codes/yambo-4.1.4/bin"
$DIRY/ypp -F amp.in -J bse
$DIRY/ypp_feng -e s -J bse
