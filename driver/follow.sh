#!/bin/bash

# This script lets you follow the current output of a running matlab session
# exp:    
#     sh follow.sh main_driver.m

sin=$1 

if [[ "$1" == *.m ]]
then
   s=${sin:0:-2}
else
   s=$sin
fi
tail -f ./out/out_$s 
