#!/bin/bash
#echo "INSTALLING SMOGV2"
#echo
#echo "NOTE: MAKE SURE YOU HAVE INSTALLED XML::Simple, Data::Dumper, Exporter and perl PDL"
#echo
#echo "SETTING UP PROPER ENVIRONMENTAL VARIABLES"
curPath=$1
#echo "Setting Vars to $curPath"
export PATH=$PATH:$curPath
export PERLLIB=$PERLLIB:$curPath
export SMOG_PATH=$curPath
#echo "DONE INSTALLING"

