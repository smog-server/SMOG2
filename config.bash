#!/bin/bash
#echo "Configuring SMOG V 2.Beta"
#echo
#echo "NOTE: MAKE SURE YOU HAVE INSTALLED String::Util, XML::Simple, Data::Dumper, Exporter and perl PDL"
#echo
#echo "SETTING UP PROPER ENVIRONMENTAL VARIABLES"
smog2Path=$1
#echo "Setting Vars to $curPath"
export PATH=$PATH:$smog2Path
export PERLLIB=$PERLLIB:$smog2Path/modules
export SMOG_PATH=$smog2Path
#echo "DONE INSTALLING"

