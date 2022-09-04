echo Detecting configured version of SMOG2
# grab all of the environment settings used when running smog2
export smog_exec=`which smog2`
if [ -z "$smog_exec" ]
then
        echo "

COULD NOT FIND USABLE VERSION OF SMOG 2!  QUITTING

"
        exit 1
fi
SMOG_PATH=`tail -n 1 $smog_exec| awk '{print $2}' | sed "s/src\/smogv2$//g"`
echo Found SMOG2 at $SMOG_PATH
SCM=$SMOG_PATH/src/tools/SCM.jar
TEMPLATE_DIR=`pwd`/share/scm/templates

problem=0

for name in CI2
do

	dir=$TEMPLATE_DIR/$name
	
	for coarse in AA AACA CA
	do
		if [ ! -f $dir/$name.$coarse.contacts ]
		then
			echo template files are missing from $TEMPLATE_DIR
			exit 1
		fi
		
		java -jar $SCM -g $dir/$name.gro -t $dir/$name.top --default -o $name.$coarse.contacts --coarse $coarse &> scm.output
		
		scmfinished=$(grep Finished scm.output | wc | awk '{print $1}')
		if [ $scmfinished -lt 1 ]
		then
			echo SCM.jar didnt finish running. Check that it exists in $SMOG_PATH/src/tool or look at scm.output for run output.
			exit 1
		fi
		
		if [ ! -f $name.$coarse.contacts ]
		then
			echo SCM.jar did not print an output contact file. Look at scm.output for run output.
			exit 1
		fi
		
		
		numlines=$(diff $dir/$name.$coarse.contacts $name.$coarse.contacts | wc | awk '{print $1}')
		if [ $numlines -gt 0 ]
		then
			echo option --coarse $coarse differed from template
			problem=1
		fi
		
		rm $name.$coarse.contacts
	done
done

echo
if [ $problem -gt 0 ]
then
	echo There were issues!
	exit 1
else
	echo "PASSED ALL"
	rm scm.output
fi


