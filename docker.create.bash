#!/bin/bash
# this script will build a new docker container from the current smog2 repo
if command -v gcp &> /dev/null;
then
        CP=gcp
else
        CP=cp
fi

mkdir docker.tmp

for i in `cat MANIFEST`
do
	$CP --parent $i docker.tmp
done

if [ ! -n "$1" ]
then
	git show --format="%H" --no-patch > docker.tmp/.gitcommit
elif [ "$1" != "release" ]
then
	echo ERROR: Arg $1 not recognized.  Only supported value is "release"
	echo BUILD NOT COMPLETED
	exit 
fi

docker build --no-cache -t smog2:gitversion .
rm -r docker.tmp
