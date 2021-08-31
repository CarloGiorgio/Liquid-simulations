#!/bin/bash
g++ -std=c++11 prova.cpp -o NPT_brute -O3 -Wall
T=1.4
rho0=0.3
for i in $(seq 0.1 0.1 1)
do
    ./NPT_brute "$rho0" "$T" "$i" data/NPT/NPT_"$i"_"$T"_brute.txt
    echo "Done pressure $i "
done
python3 eos_NPT.py "$T"_brute True
