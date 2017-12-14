#!/bin/bash
SIM_HASH_LEN=64
SIM_HASH_NUM=2
export SIM_HASH_LEN
export SIM_HASH_NUM
touch newsieve/Typedefs.h
make

k=3
d=60
b=30

filename_out=sieve_k${k}_d${d}_b${b}
time nohup ./newlatsieve -k${k} -d${d} -b${b} >$filename_out 2>&1 &
