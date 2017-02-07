#!/bin/bash
#   created by: 
#        Johannes Becherer
#        Tue Feb  7 11:35:24 PST 2017
#
#   This function will reset the mfile directory to its pre initialization stag


echo Are you sure that you want to reset the mfile directory?
echo this action will delete all driver routines from this directory
echo type y/n

read test

if [[ $test == y* ]]
then
   echo !!!!!!!!!!!! GOOOOOO !!!!!!!!!
      rm -rf ./chipod_gust
      rm -r ./out

      rm ./*.m
      rm run.sh
      rm kill_script.sh
      rm monitorMatlab.sh
else
   echo action stoped
fi



