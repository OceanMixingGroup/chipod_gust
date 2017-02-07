#!/bin/bash

if [[ $(hostname) == matlab* ]]
then
watch -n 5 "echo 'matlab1';echo; ssh matlab1 'ps aux |  grep matlab | grep -v grep';echo;echo 'matlab2';echo;ssh matlab2 'ps aux |  grep matlab | grep -v grep';echo; echo 'matlab3';echo;ssh matlab3 'ps aux |  grep matlab | grep -v grep';echo;echo 'matlab4';echo;ssh matlab4 'ps aux | grep matlab | grep -v grep';echo; echo 'matlab5';echo;ssh matlab5 'ps aux |  grep matlab | grep -v grep ';"
else
watch -n 5 "ps aux |  grep matlab | grep -v grep"  
fi
