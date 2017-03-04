#!/bin/bash

# This script kills a currently running matlab session
# exp:    
#     sh kill_script.sh main_driver.m

s=$1
if [[ "$s" == *.m ]]
then
   s=${s%?}
   s=${s%?}
fi

kill  $(cat ./out/pid_$s)
