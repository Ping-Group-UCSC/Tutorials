#!/bin/bash
grep "S:" lindbladLinear.out > ftmp
awk '{print $5"  "$12"  "$13"  "$14}' ftmp > st.dat
sed -i '1d' st.dat
rm ftmp
