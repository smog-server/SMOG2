#!/bin/bash
# this script will build a new docker container from the current smog2 repo
containername=smog2:gitversion
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

docker build --no-cache -t $containername .
rm -r docker.tmp

echo "
Done building container $containername

It is strongly recommended that you run smog-check
to ensure that the container is properly configured/

To do this, run the container. Once inside, issue the commands:
cd /opt/smog2/SMOG-CHECK
./smog-check
./smog-tool-check

Since containers are static, you will only need to run
smog-check the first time.

"
