#!/bin/bash
#########################################
# This is just the quick_check script
# that only makes the tops instead of 
# checking them
#########################################

export smog_exec=`which smog2`
if [ -z "$smog_exec" ]
then
	echo "

COULD NOT FIND USABLE VERSION OF SMOG 2!  QUITTING

"
	exit 1
fi

export perl4smog=`tail -n 1 $smog_exec| awk '{print $1}'`

export SMOG_PATH=`tail -n 1 $smog_exec| awk '{print $2}' | sed "s/src\/smogv2$//g"`

echo executable to be tested: $smog_exec
echo perl to be used: $perl4smog
echo smog path: $SMOG_PATH

dir=$(pwd)
count=0
if [ -z "$SMOG_PATH" ]; then
	echo ERROR: quick-check needs to know where SMOG2 is
	echo Use configure script in SMOG2 distribution to set up SMOG_PATH
	exit 1
fi
echo "Using the executables found in $SMOG_PATH"

for pdb in 1AKEapo_v2 1ypa.2chains 2GIS_noSAM_v2
do
for graining in AA CA AAgauss
do

if [ -e $pdb.$graining.pdb ]
then

echo Making $pdb with a $graining graining

if [ $graining == "AA" ]
then
	grainingText="AA"
elif [ $graining == "CA" ]
then
	grainingText="CA"
elif [ $graining == "AAgauss" ]
then
	grainingText="t $SMOG_PATH/share/templates/SBM_AA+gaussian"
fi

$perl4smog $smog_exec -i $pdb.$graining.pdb -$grainingText -dname temporary  &> output
mv temporary.top $pdb.$graining.top
rm temporary*

fi

done
done

rm output
