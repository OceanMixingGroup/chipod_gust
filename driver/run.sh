#!/bin/bash

nohup matlab.2014b -nodisplay -nosplash -r $1 > ./out/out_$1 & echo $! > ./out/pid_$1
