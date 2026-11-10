#!/bin/bash
# run.sh sd | <scale> | <s1> <s2>
cd $(dirname $(readlink -f $0))
if [ "$1" == sd ]; then d=sd; mkdir -p $d; cp base.d12 $d/INPUT
elif [ -n "$2" ]; then d=s$1_$2; mkdir -p $d; python3 mk2.py $1 $2 $d/INPUT
else d=s$1; mkdir -p $d; python3 mk.py $1 $d/INPUT; fi
cd $d
OMP_NUM_THREADS=${NT:-4} /home/palacios/CRYSTAL23OMP/bin/Linux-ifx_i64_omp/v1.0.1/crystalOMP < INPUT > out 2>&1
e=$(grep "TOTAL ENERGY(DFT)" out | tail -1 | awk '{print $4}')
sp=$(grep -A1 "TOTAL ATOMIC SPINS" out | tail -1 | awk "{print \$1}")
echo "$d E=$e spin=$sp cyc=$(grep -c 'CYC  ' out) $(grep -m1 -iE 'ERROR|DEPENDEN' out)"
