#!/bin/bash
g++ -std=c++11 LJ_project.cpp -lm -D ense=0 -D gr -Wall -O3 -o NTV_gr.x
T=1.4
bins=80

while IFS="" read -r i || [ -n "$i" ]
do
  echo "Density $i "
  ./NTV_gr.x "$i" "$T" "$bins" data/NTV/NTV_"$i"_"$T"_g.txt data/gr/gr_"$i"_"$T"_.txt
done < dens.txt
python3 gr.py



#for i in $(seq 0.05 0.01 0.3)
#do
#    ./NTV_gr.x "$i" "$T" "$bins" data/NTV/NTV_"$i"_"$T"_g.txt data/gr/gr_"$i"_"$T"_.txt
#    echo "Done density $i "
    #python3 gr.py data/gr/gr_"$i"_"$T"_.txt
#done
#python3 eos_NVT.py "$T"_g True
