#!/bin/bash
# Ni sd reference, then the uniform 4p scale scan (as for Fe), results appended to all.txt
cd $(dirname $(readlink -f $0))
./run.sh sd >> all.txt
for l in 0.6 0.5 0.7 0.4 0.55 0.65 0.8 1.0 1.2 1.5 2.0 3.0; do ./run.sh $l >> all.txt; done
echo SCANDONE >> all.txt
