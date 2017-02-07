#!/bin/bash
#   created by: 
#        Johannes Becherer
#        Tue Feb  7 11:35:24 PST 2017
#
#        This shell script initalizes the mfile folder  of a chipod/gust directory by
#           1) downloading the latest software versoin from github
#           2) creating an out put directory
#           3) copying the standard driver function in this directory 

git clone https://github.com/OceanMixingGroup/chipod_gust

mkdir out

cp ./chipod_gust/driver/* ./
