#!/bin/bash
# this script will build a new docker container from the current smog2 repo
# multiple tags can be added by adding more -t flags at build time
tag=$1
option=$2
docker_username="smogserver"
containername=$docker_username/smog2:$tag
if [ "$option" == "cross" ]
then
	# this build a cross-platform container. It also pushes the container to your DockerHub repo. For a simpler build, see the commented command, below
	echo Will attempt to build and push a cross-platform version of the SMOG 2 container for DockerHub account $docker_username and SMOG version $tag
	docker buildx build --push --platform linux/arm/v7,linux/arm64/v8,linux/amd64 --no-cache -t $containername -f Dockerfile.$tag .
elif [ "$option" == "" ]
then
	echo Will attempt to build a local native-platform version of the SMOG 2 container: container called $containername
	docker build --no-cache -t $containername -f Dockerfile.$tag .
else

	echo \"$option\" is not a supported option 
	exit
fi

if [ $? == 0 ]
then

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

else

echo "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BUILD OF DOCKER CONTAINER DID NOT COMPLETE PROPERLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"

fi
