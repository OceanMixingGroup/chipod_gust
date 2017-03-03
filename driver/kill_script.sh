#!/bin/bash

# This script kills a currently running matlab session
# exp:    
#     sh kill_script.sh main_driver.m

sin=$1 

if [[ "$1" == *.m ]]
then
   s=${sin:0:-2}
else
   s=$sin
fi
kill  $(cat ./out/pid_$s)
