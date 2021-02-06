#!/bin/bash

if [ $# -lt 1 ]
then

if [ -d distro ]
then
	echo directory distro exists. quitting
else
	mkdir distro
	for name in `cat MANIFEST`
	do
		if command -v gcp &> /dev/null;
		then
			gcp --parent $name distro
		else
			cp --parent $name distro
		fi
	done
	sed "s/optim//g" configure.smog2 > distro/configure.smog2
	echo done making distro

	echo DONT FORGET TO UPDATE VERSION NUMBER IN smogv2
fi

else

if [ $1 = manifest ]
then

# Steps:
# 1. makes a don't keep file with names of files
#     to remove
# 2. build grep remove line
# 3. add all files and remove the files that we don't want 
# 4. go back and check MANIFEST yourself and make sure 
#     nothing unwanted made it in

if [ -f MANIFEST.bak ]
then
echo "ERROR: Backup MANIFEST already exists"
echo "Quitting so that I don't overwrite it"
exit
fi

cp MANIFEST MANIFEST.bak
#^ means match the beginning in grep
cat > dontkeep << EOF
^./SBM_AA
^./SBM_AA+gaussian
^./SBM_calpha
^./SBM_calpha+gaussian
^./bin
^.git
^./dontkeep
^./MANIFEST.bak
EOF

remove="dummyXYZ"
for token in `cat dontkeep`
do
remove=$remove"\\|$token"
done
rm dontkeep

find . -type f | grep -v "$remove" | sed "s/\.\///g" > MANIFEST

else
echo Error: parameter $1 not recognized
echo "bash make_distro.bash manifest <-- builds a new MANIFEST file"
fi
fi
