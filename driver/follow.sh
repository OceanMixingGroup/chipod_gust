#!/bin/bash

# This script lets you follow the current output of a running matlab session
# exp:    
#     sh follow.sh main_driver.m

s=$1
if [[ "$s" == *.m ]]
then
   s=${s%?}
   s=${s%?}
fi

tail -f ./out/out_$s 
