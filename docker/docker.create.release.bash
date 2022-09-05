#!/bin/bash
# this script will build a new docker container from the current smog2 repo
tag=v2.4.4
containername=smogserver/smog2:$tag
docker build --no-cache -t $containername -f Dockerfile.$tag .

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
