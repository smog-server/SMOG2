#!/bin/bash
#echo "Configuring SMOGV2"
#echo
#echo "NOTE: MAKE SURE YOU HAVE INSTALLED String::Util, XML::Simple, Data::Dumper, Exporter and perl PDL"
#echo
#echo "SETTING UP PROPER ENVIRONMENTAL VARIABLES"
curPath=$1
#echo "Setting Vars to $curPath"
export PATH=$PATH:$curPath
export PERLLIB=$PERLLIB:$curPath/modules
export SMOG_PATH=$curPath
#echo "DONE INSTALLING"

