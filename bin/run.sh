#! /bin/bash

./build.py -j5

./macsim --l3_repl==PDP <a.txt | tee -i PDP.out

mkdir PDP_old_16
mv *.out* ./PDP_old_16
mv *.csv* ./PDP_old_16

#./macsim | tee -i LRU.out

#mkdir LRU2
#mv *.out* ./LRU2
