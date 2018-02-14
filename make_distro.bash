#!/bin/bash

if [ -d distro ]
then
	echo directory distro exists. quitting
else

	mkdir distro
	for name in `cat MANIFEST`
	do
		cp --parent $name distro
		sed "s/optim//g" configure.smog2 > distro/configure.smog2
	done
	echo done making distro

	echo DONT FORGET TO UPDATE VERSION NUMBER IN smogv2
fi
