#!/bin/bash
#   created by: 
#        Johannes Becherer
#        Tue Feb  7 11:35:24 PST 2017
#
#        This shell script initalizes the mfile folder  of a chipod/gust directory by
#           1) downloading the latest software versoin from github
#           2) creating an out put directory
#           3) copying the standard driver function in this directory 


if [[ $(hostname) == matlab* ]]
then
      echo You are logged on to the matlab server! 
      echo The github connection will be done via 128.193.69.189
      ssh   mixing@128.193.69.189 'cd ~/ganges/work/chipod_gust; sh update.sh'
      scp -r mixing@128.193.69.189:~/ganges/work/chipod_gust ./
                  
else
      echo The chipod_gust software package is loaded from github
      
      git clone https://github.com/OceanMixingGroup/chipod_gust
fi
   

mkdir out
cp ./chipod_gust/driver/* ./
