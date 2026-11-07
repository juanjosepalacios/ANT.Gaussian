#!/bin/bash
d=t$1_$2; mkdir -p $d; cd $d; python3 ../mk2.py $1 $2 INPUT
OMP_NUM_THREADS=8 /home/palacios/CRYSTAL23OMP/bin/Linux-ifx_i64_omp/v1.0.1/crystalOMP < INPUT > out 2>&1
echo "$1 $2 $(grep 'TOTAL ENERGY(DFT)' out | tail -1 | awk '{print $4}') $(grep -A1 'TOTAL ATOMIC SPINS' out | tail -1)"
