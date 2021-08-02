#########################################################################################
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#            Jeffrey Noel, Mariana Levi, Antonio Oliveira, Vinícius Contessoto,
#             Mohit Raghunathan, Joyce Yang, Prasad Bandarkar, Udayan Mohanty,
#                          Ailun Wang, Heiko Lammert, Ryan Hayes,
#                               Jose Onuchic & Paul Whitford
#
#            Copyright (c) 2015,2016,2018,2021, The SMOG development team at
#                        Rice University and Northeastern University
#
#          SMOG 2, Shadow and openSMOG are available at http://smog-server.org
#
#          You can direct questions to info@smog-server.org, or the smog-users forum,
#          which you can find at https://mailman.rice.edu/mailman/listinfo/smog-users
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#########################################################################################

#####################################################################
# setRatios: Adjust ratios of strengths as obtained from .sif
#####################################################################
package setRatios;
use strict;
use warnings FATAL => 'all';
use Exporter;
use PDL;
use templateParser;
use smog_common;

our @ISA = 'Exporter';
our @EXPORT = qw(setRatios getSetDiheCounts %fTypes);

sub DiheCountsHelper
{
	my %uniqueBonds;
	my($diheArr) = @_;
	my @countsIndex;
	my @counts;
	my $size = $diheArr->dim(1);
	## Count number of dihedrals passing through a bond ##
	my $tindex=0;
	$counts[0]=1;
	my @rescale;
	for(my $i=0;$i<$size;$i++)
	{
		my @A = $diheArr->slice(":,$i:$i")->list;
		if($#A == 0){
			$countsIndex[$i]=-1;
			$rescale[$i]=0;
			next;
		}elsif($#A >= 6){ 
			$rescale[$i]=1;
		}else{
			$rescale[$i]=0;
		}
		my $b=$A[1];
		my $c=$A[2];
		my $eG=$A[6];
		## only count the dihedral if it is not an improper
		if($eG >= 0 ){	
	
	       		if($b>$c){my $tt=$c;$c=$b;$b=$tt};
	
	       		if(exists $uniqueBonds{"$b-$c--$eG"}){  
	       			my $tt=$uniqueBonds{"$b-$c--$eG"};
	       			$countsIndex[$i]=$tt;
	       			$counts[$tt]++;
	       		}else{
	       			$tindex++;
	       			$uniqueBonds{"$b-$c--$eG"}=$tindex;
	       			$countsIndex[$i]=$tindex;
	       			$counts[$tindex]++;
	       		}
		}else{
			$countsIndex[$i]=0;		
		}	
	}
	
	for(my $i=0;$i<$size;$i++)
	{
		if($rescale[$i]==1){
			set($diheArr,5,$i,1/$counts[$countsIndex[$i]]);
		}
	}
	undef %uniqueBonds;
}

# For each chain, count the total number of bonds through a dihedral,
# utilizes getDiheCounts and setDiheCounts
sub getSetDiheCounts
{
	my($diheFunctHandle) = @_;
	foreach my $chain(keys %{$diheFunctHandle})
	{
		DiheCountsHelper($diheFunctHandle->{$chain});
	}	
}

# Set the dihedral strength through normalization.
sub setRatios	
{
	my($diheFunctHandle,$inputPDL,$atomNum,$atomTypes) = @_;
	my $sum; my $diheStrengthTotal;
	my %rescalePDL;
	foreach my $chain(keys %{$diheFunctHandle})
	{
		adjustFactorsHelper1($diheFunctHandle->{$chain},$inputPDL->{$chain},$atomTypes,\$sum,\%rescalePDL,$chain);
	}
	
	foreach my $chain(keys %{$diheFunctHandle})
	{
		adjustFactorsHelper2($diheFunctHandle->{$chain},$inputPDL->{$chain},$atomNum,$atomTypes,$diheStrengthTotal,$sum,\%rescalePDL,$chain);
	}
}

sub adjustFactorsHelper1
{
	my($diheArr,$inputPDL,$atomTypes,$sum,$rescalePDL,$chain) = @_;
	my $normalize;
	my $relativeRatio;
 	my $size = $diheArr->dim(1);
	my @tempArr;
	my $tmpsum=${$sum};
 	for(my $i=0;$i<$size;$i++)
 	{

		my @buffer=$diheArr->slice(":,$i:$i")->list;
        	if($#buffer <6){
			next;
		}
        	my($a,$b,$c,$d,$func,$count,$eG) = @buffer;

		## Convert from relative index to absolute indexing ##
		$b = sclr(slice($inputPDL,"3:3,$b,:"));
		$c = sclr(slice($inputPDL,"3:3,$c,:"));
		my $resTypeb = $atomTypes->{$b}->[1];
    		my $resTypeUse;

 		## $eG == -1 is IMPROPER SKIP ##
		if($eG < 0) 
		{			
			push(@tempArr,pdl($count,0));
			next;
		}
		$eG = $eGTable{$eG}; 
   		if(defined $termRatios->{$resTypeb}->{"energyGroup"}->{$eG})
   		{
       			$normalize = $termRatios->{$resTypeb}->{"energyGroup"}->{$eG}->{"normalize"};
       			$resTypeUse = $resTypeb;
		}else{
       			## CASE WHERE THERE IS A DIHEDRAL BETWEEN TWO DIFFERENT RESTYPES ##
			my $resTypec = $atomTypes->{$c}->[1];
       			if(! defined $termRatios->{$resTypec}->{"energyGroup"}->{$eG})
       			{
    				$a = sclr(slice($inputPDL,"3:3,$a,:"));
				$c = sclr(slice($inputPDL,"3:3,$c,:"));
				$d = sclr(slice($inputPDL,"3:3,$d,:"));
        	    		smog_quit("energyGroup $eG is not defined for $resTypeb-$resTypec ($a-$b-$c-$d)\n");
       			}
       			$normalize = $termRatios->{$resTypec}->{"energyGroup"}->{$eG}->{"normalize"};

	 		if(!defined $normalize)
			{
	    			$a = sclr(slice($inputPDL,"3:3,$a,:"));
	    			$c = sclr(slice($inputPDL,"3:3,$c,:"));
				$d = sclr(slice($inputPDL,"3:3,$d,:"));
	    			my $resTypea = $atomTypes->{$a}->[1];
	    			my $resTypec = $atomTypes->{$c}->[1];
	    			my $resTyped = $atomTypes->{$d}->[1];
				smog_quit("Normalize option not set for $resTypea-$resTypeb-$resTypec-$resTyped of energyGroup $eG with atom indices $a-$b-$c-$d");
			}
       			$resTypeUse = $resTypec;
   		}

    		## Normalize option is set ##	
		if($normalize eq 1) 
		{
			$relativeRatio=$termRatios->{$resTypeUse}->{"energyGroup"}->{$eG}->{"intraRelativeStrength"};
			$count*=$relativeRatio;
			if($residues{$atomTypes->{$b}->[5]}->{'atomCount'} !=0 and $residues{$atomTypes->{$c}->[5]}->{'atomCount'}!=0){
        			$tmpsum+=$count;
			}
        		set($diheArr,5,$i,$count);
			push(@tempArr,pdl($count,1));
		}else{
			push(@tempArr,pdl($count,0));
		}
	
 	}
	if(@tempArr){
		$rescalePDL->{$chain}=cat(@tempArr);
 	}
        ${$sum}=$tmpsum;
}

sub adjustFactorsHelper2
{
	my($diheArr,$inputPDL,$totalAtoms,$atomTypes,$diheStrengthTotal,$sum,$rescalePDL,$chain) = @_;
	my $totalStrength;
	my $contactTotal;
 	my $size = $diheArr->dim(1);
	my $diheLeftOver = 0;
        ## epsilonC+epsilonD ##
	$totalStrength = $termRatios->{"interRelativeTotal"};
	## epsilonC ##			
	$contactTotal = $termRatios->{"contactRelative"};
	## leftOver = totalAtoms*(1-epsilonC/(epsilonC+epsilonD)) ##			
	$diheLeftOver = $totalAtoms*(1 - ($contactTotal/$totalStrength));


	if(!defined $rescalePDL->{$chain}){
		return;
	}
	for(my $i=0;$i<$size;$i++)
	{
		my $normalize=sclr(slice($rescalePDL->{$chain},"1:1,$i:$i"));

    		## Normalize option is set ##	
		if($normalize eq 1) 
		{
			my $count=sclr(slice($rescalePDL->{$chain},"0:0,$i:$i"));
			$count = ($count/$sum)*($diheLeftOver); 
			set($diheArr,5,$i,$count);
		}
	}
}

1;
