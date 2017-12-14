#!/bin/bash
SIM_HASH_LEN=64
SIM_HASH_NUM=2
export SIM_HASH_LEN
export SIM_HASH_NUM
touch newsieve/Typedefs.h
make
