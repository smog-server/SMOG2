#!/bin/bash

if [ -d distro ]
then
	echo directory distro exists. quitting
else

	mkdir distro
	for name in `cat MANIFEST`
	do
		cp --parent $name distro
	done
	echo done making distro
fi
