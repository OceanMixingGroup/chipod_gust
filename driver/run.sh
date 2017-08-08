#!/bin/bash

# This script executes any matlab script in a nohup status
# exp:    
#     sh run.sh main_driver.m

s=$1
if [[ "$s" == *.m ]]
then
   s=${s%?}
   s=${s%?}
fi

nohup matlab.2017a -nodisplay -nosplash -r $s > ./out/out_$s & echo $! > ./out/pid_$s
