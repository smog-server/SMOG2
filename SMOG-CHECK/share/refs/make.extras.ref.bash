#!/bin/bash
for i in atomtypes bondtypes angletypes dihedraltypes nonbond_params
do
	grep "^$i" ../templates/SBM_AA/extras > extras.ref.$i 
done
