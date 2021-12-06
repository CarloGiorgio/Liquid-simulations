#!/bin/bash
g++ -std=c++17 LJ_project.cpp -lm -D ense=1 -Wall -O3 -o NPT.x

T=1.4
rho0=0.1
for i in $(seq 0.01 0.01 0.24)
do
    ./NPT.x "$rho0" "$T" "$i" data/NPT/NPT_"$i"_"$T"_.txt
    echo "Done pressure $i "
done

#while read rho T p;do
#./NPT.x "$rho" "$T" "$p" data/NPT/NPT_"$p"_"$T"_.txt
#echo "Done pressure $p "
#done<"eos_NTV_1.4_.txt"

python3 eos_NPT.py "$T"_ NPT True