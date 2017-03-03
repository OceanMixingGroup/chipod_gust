#!/bin/bash

# This script executes any matlab script in a nohup status
# exp:    
#     sh run.sh main_driver.m

sin=$1 

if [[ "$1" == *.m ]]
then
   s=${sin:0:-2}
else
   s=$sin
fi

nohup matlab.2014b -nodisplay -nosplash -r $s > ./out/out_$s & echo $! > ./out/pid_$s
