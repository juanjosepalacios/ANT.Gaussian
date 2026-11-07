#!/bin/bash
lam=$1; d=s$lam; mkdir -p $d; cd $d; python3 ../mk.py $lam INPUT
OMP_NUM_THREADS=8 /home/palacios/CRYSTAL23OMP/bin/Linux-ifx_i64_omp/v1.0.1/crystalOMP < INPUT > out 2>&1
e=$(grep "TOTAL ENERGY(DFT)" out | tail -1 | awk '{print $4}')
sp=$(grep -A1 "TOTAL ATOMIC SPINS" out | tail -1 | awk "{print \$1}")
echo "$lam E=$e spin=$sp cyc=$(grep -c 'CYC  ' out) $(grep -m1 -iE 'ERROR|DEPENDEN' out)"
