#!/bin/bash

#set this to your directory
SCM=../SCM.jar
problem=0

for name in CI2
do

dir=templates/$name

for coarse in AA AACA CA
do

java -jar $SCM -g $dir/$name.gro -t $dir/$name.top --default -o $name.$coarse.contacts --coarse $coarse &> /dev/null

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
else
echo Simple test passed!
fi
echo
