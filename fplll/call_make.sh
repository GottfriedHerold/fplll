#!/bin/sh

cd fplll/fplll/
SIM_HASH_LEN=64
SIM_HASH_NUM=2
THRESHOLD_LVLS_2SIEVE_LB={{32-7, 64-12}}

export SIM_HASH_LEN
export SIM_HASH_NUM
touch newsieve/Typedefs.h


k=3
d=64
b=30

cd ../../
filename_out=sieve_k${k}_d${d}_b${b}
make -f Makefile >> $filename_out &
time fplll/fplll/newlatsieve -k${k} -d${d} -b${b} >$filename_out 2>&1 
