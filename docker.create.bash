#!/bin/bash
# this script will build a new docker container from the smog2 repo.
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
docker build --no-cache -t smog2:latest .
rm -r docker.tmp
