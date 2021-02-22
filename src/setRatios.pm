#########################################################################################
#
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#                     Jeffrey Noel, Mariana Levi, Mohit Ranghunathan,
#                 Heiko Lammert, Ryan Hayes, Jose Onuchic & Paul Whitford
#
#                     Copyright (c) 2015, The SMOG development team at
#                        Rice University and Northeastern University
#
#              SMOG 2 & Shadow are available at http://smog-server.org
#
#                        Direct questions to: info@smog-server.org
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
my %uniqueBonds;

sub DiheCountsHelper
{
	my($diheArr,$inputPDL) = @_;
	my @countsIndex;
	my @counts;
	my $size = $diheArr->dim(1);
	## Count number of dihedrals passing through a bond ##
	my $tindex=0;
	$counts[0]=1;
	for(my $i=0;$i<$size;$i++)
	{
		my @A = $diheArr->slice(":,$i:$i")->list;
		if($#A == 0){
			$countsIndex[$i]=-1;
			next;
		} 
		my $b=$A[1];
		my $c=$A[2];
		my $eG=$A[6];
		$b = sclr(slice($inputPDL,"3:3,$b,:"));
		$c = sclr(slice($inputPDL,"3:3,$c,:"));
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
		if($countsIndex[$i]==-1){
			next;
		}
		my @A = $diheArr->slice(":,$i:$i")->list;
		unless($#A < 6){
			set($diheArr,5,$i,1/$counts[$countsIndex[$i]]);
		}
	}
}

# For each chain, count the total number of bonds through a dihedral,
# utilizes getDiheCounts and setDiheCounts
sub getSetDiheCounts
{
	my($diheFunctHandle,$whichPDL) = @_;
	foreach my $chain(keys %{$diheFunctHandle})
	{
		DiheCountsHelper($diheFunctHandle->{$chain},$whichPDL->{$chain});
	}	
}

# Set the dihedral strength through normalization.
sub setRatios	
{
	my($diheFunctHandle,$inputPDL,$atomNum,$atomTypes) = @_;
	my %residueRatio;
	my $sum; my $diheStrengthTotal;
	undef %uniqueBonds;
	my %rescalePDL;
	foreach my $chain(keys %{$diheFunctHandle})
	{
		adjustFactorsHelper1($diheFunctHandle->{$chain},$inputPDL->{$chain},$atomNum,$atomTypes,\%residueRatio,\$sum,\%rescalePDL,$chain);
	}
	
	foreach my $chain(keys %{$diheFunctHandle})
	{
		adjustFactorsHelper2($diheFunctHandle->{$chain},$inputPDL->{$chain},$atomNum,$atomTypes,$diheStrengthTotal,\$sum,\%rescalePDL,$chain);
	}
}

sub adjustFactorsHelper1
{
	my($diheArr,$inputPDL,$totalAtoms,$atomTypes,$residueRatio,$sum,$rescalePDL,$chain) = @_;
	my $normalize;
	my $relativeRatio;
 	my $size = $diheArr->dim(1);
	my @tempArr;
 	for(my $i=0;$i<$size;$i++)
 	{

		my @buffer=$diheArr->slice(":,$i:$i")->list;
        	if($#buffer <6){next;}
        	my($a,$b,$c,$d,$func,$count,$eG) = @buffer;

		## Convert from relative index to absolute indexing ##
    		$a = sclr(slice($inputPDL,"3:3,$a,:"));
		$b = sclr(slice($inputPDL,"3:3,$b,:"));
		$c = sclr(slice($inputPDL,"3:3,$c,:"));
		$d = sclr(slice($inputPDL,"3:3,$d,:"));
		my $resnamea = $atomTypes->{$a}->[5];
		my $resnameb = $atomTypes->{$b}->[5];
		my $resnamec = $atomTypes->{$c}->[5];
		my $resnamed = $atomTypes->{$d}->[5];
		my $resTypeb = $atomTypes->{$b}->[1];
		my $resTypec = $atomTypes->{$c}->[1];
    		my $resTypeUse;

 		## $eG == -1 is IMPROPER SKIP ##
		if($eG < 0) 
		{			
		push(@tempArr,pdl($count,0));
		 next;
		}
		$eG = $eGTable{$eG}; ## Obtain user defined residue name ##
   		if(!defined $termRatios->{$resTypeb}->{"energyGroup"}->{$eG})
   		{
       			## CASE WHERE THERE IS A DIHEDRAL BETWEEN TWO DIFFERENT RESTYPES ##
       			if(! defined $termRatios->{$resTypec}->{"energyGroup"}->{$eG})
       			{
        	    		smog_quit("energyGroup $eG is not defined for $resTypeb-$resTypec ($a-$b-$c-$d)\n");
       			}
       			$normalize = $termRatios->{$resTypec}->{"energyGroup"}->{$eG}->{"normalize"};
       			$resTypeUse = $resTypec;
   		}
   		else{
       			$normalize = $termRatios->{$resTypeb}->{"energyGroup"}->{$eG}->{"normalize"};
       			$resTypeUse = $resTypeb;
		}
 		if(!defined $normalize)
		{
    			my $resTypea = $atomTypes->{$a}->[1];
    			my $resTyped = $atomTypes->{$d}->[1];
			smog_quit("Normalize option not set for $resTypea-$resTypeb-$resTypec-$resTyped of energyGroup $eG with atom indices $a-$b-$c-$d");
		}
    		## Normalize option is set ##	
		if($normalize eq 1) 
		{
			$relativeRatio=$termRatios->{$resTypeUse}->{"energyGroup"}->{$eG}->{"intraRelativeStrength"};
			if(!defined $relativeRatio){
				smog_quit("normalize=0, but intraRelativeStrength not defined for energyGroup $eG. Check .sif file");
			}

			$count*=$relativeRatio;
			unless($residues{$resnamea}->{'atomCount'}==0 || $residues{$resnameb}->{'atomCount'}==0
				|| $residues{$resnamec}->{'atomCount'}==0 || $residues{$resnamed}->{'atomCount'}==0){
        			${$sum}+=$count;
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


	for(my $i=0;$i<$size;$i++)
	{
		if(!defined $rescalePDL->{$chain}){next;}
		my $normalize=sclr(slice($rescalePDL->{$chain},"1:1,$i:$i"));

    		## Normalize option is set ##	
		if($normalize eq 1) 
		{
			my $count=sclr(slice($rescalePDL->{$chain},"0:0,$i:$i"));
			$count = ($count/${$sum})*($diheLeftOver); 
			set($diheArr,5,$i,$count);
		}
	}
}

1;
