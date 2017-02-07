#
#/bin/bash
#   created by: 
#        Johannes Becherer
# Tue Feb  7 11:31:29 PST 2017
# 
# This shell script gets the lates version of chipod_gust software from github
#

if [[ $(hostname) == matlab* ]]
then
      echo You are logged on to the matlab server! 
      echo The github connection will be done via 128.193.69.189
      ssh   mixing@128.193.69.189 'cd ~/ganges/work/chipod_gust; bash update.sh'

      scp -r mixing@128.193.69.189:~/ganges/work/chipod_gust/* ./

else
      echo The chipod_gust software package is updated from github
      git pull origin master
fi
