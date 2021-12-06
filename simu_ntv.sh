#!/bin/bash
g++ -std=c++11 LJ_project.cpp -lm -D ense=0 -Wall -O3 -o NTV.x
T=1.4
while IFS="" read -r i || [ -n "$i" ]
do
  echo "Density $i "
    ./NTV.x "$i" "$T" data/NTV/NTV_"$i"_"$T"_.txt
done < dens.txt
python3 eos_NVT.py "$T"_ True

#for i in $(seq 0.05 0.01 0.3)
#do
#    ./NTV.x "$i" "$T" data/NTV/NTV_"$i"_"$T"_.txt
#    echo "Done density $i "
#done
#python3 eos_NVT.py "$T"_ True
