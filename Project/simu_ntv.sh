#!/bin/bash
g++ -std=c++11 prova.cpp -o NTV_LL -O3 -Wall
T=1.4
for i in $(seq 0.05 0.025 0.6)
do
    ./NTV_LL "$i" "$T" data/NTV/NTV_"$i"_"$T"_LL.txt
    echo "Done density $i "
done
python3 eos_NVT.py "$T"_LL True
