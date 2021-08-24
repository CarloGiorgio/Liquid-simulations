#!/bin/bash
g++ -std=c++11 prova.cpp -o prova -O3 -Wall
T=1.5
rho0=0.3
for i in $(seq 0.03 0.2 4)
do
    ./prova "$rho0" "$T" "$i" NTP_"$i"_"$T".txt
    echo "Done pressure $i "
done
