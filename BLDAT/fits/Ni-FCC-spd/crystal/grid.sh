#!/bin/bash
cd $(dirname $(readlink -f $0))
for a in 0.45 0.6 0.75; do for c in 0.45 0.6 0.75; do
  [ $a == 0.6 ] && [ $c == 0.6 ] && continue; ./run.sh $a $c >> grid.txt; done; done
echo GRIDDONE >> grid.txt
