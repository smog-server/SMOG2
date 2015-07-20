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
#              SMOG v2 & Shadow are available at http://smog-server.org
#
#                        Direct questions to: info@smog-server.org
#
#########################################################################################

#####################################################################
# setRatios: Adjust ratios of strengths as obtained from .sif
#####################################################################
package setRatios;
#####################
## COMPILE HEADERS ##
#####################
use strict;
use warnings;
####################
## MODULE HEADERS ##
####################
use Exporter;
use PDL;
use templateParser;

our @ISA = 'Exporter';
our @EXPORT = qw(setRatios getSetDiheCounts %fTypes);
my %uniqueBonds;

sub getDiheCountsHelper
{
 my($diheArr,$inputPDL) = @_;
 my $size = $diheArr->dim(1);
 my @tempArr;
 ## Count number of dihedrals passing through a bond ##
 for(my $i=0;$i<$size;$i++)
 {
 	my @A = $diheArr->slice(":,$i:$i")->list;
 	if($#A == 0){next;} 
	my $a=$A[0];
	my $b=$A[1];
	my $c=$A[2];
	my $d=$A[3];
	my $func=$A[4];
	my $cD=$A[5];
	my $eG=$A[6];
	$a = sclr(slice($inputPDL,"3:3,$a,:"));
	$b = sclr(slice($inputPDL,"3:3,$b,:"));
	$c = sclr(slice($inputPDL,"3:3,$c,:"));
	$d = sclr(slice($inputPDL,"3:3,$d,:"));

 ## only count the dihedral if it is not an improper
        if($eG >= 0 ){	
		if(exists $uniqueBonds{"$b-$c--$eG"}){$uniqueBonds{"$b-$c--$eG"}++;}
		elsif (exists $uniqueBonds{"$c-$b--$eG"}) {$uniqueBonds{"$c-$b--$eG"}++;}
		else {$uniqueBonds{"$b-$c--$eG"}=1;}
	}	
 }
}

sub setDiheCountsHelper
{
 my($diheArr,$inputPDL) = @_;
 my $size = $diheArr->dim(1);
 my $count=0;
 
	for(my $i=0;$i<$size;$i++)
	{
	#	my ($a,$b,$c,$d,$func,$cD,$eG) = $diheArr->slice(":,$i:$i")->list;
		my @A = $diheArr->slice(":,$i:$i")->list; 
 		if($#A == 0){next;} 
#		if($#A >0){ 
			my $a=$A[0];
			my $b=$A[1];
			my $c=$A[2];
			my $d=$A[3];
			my $func=$A[4];
			my $cD=$A[5];
			my $eG=$A[6];

			$a = sclr(slice($inputPDL,"3:3,$a,:"));
			$b = sclr(slice($inputPDL,"3:3,$b,:"));
			$c = sclr(slice($inputPDL,"3:3,$c,:"));
			$d = sclr(slice($inputPDL,"3:3,$d,:"));
			
			$count = (exists $uniqueBonds{"$b-$c--$eG"}?
					$uniqueBonds{"$b-$c--$eG"}
					:$uniqueBonds{"$c-$b--$eG"});

			if($eG >=0){
				set($diheArr,5,$i,1/$count);
			}else{
				set($diheArr,5,$i,1);
			}
#	 	}
	}
}

##
# Count the total number of bonds through a dihedral,
# utilizes getDiheCounts and setDiheCounts
sub getSetDiheCounts
{
 my($diheFunctHandle,$whichPDL) = @_;
 foreach my $res(keys %{$diheFunctHandle})
 {
	getDiheCountsHelper($diheFunctHandle->{$res},$whichPDL->{$res});
 }
 foreach my $res(keys %{$diheFunctHandle})
 {
	setDiheCountsHelper($diheFunctHandle->{$res},$whichPDL->{$res});
 
 }	
}

##
# Set the dihedral strength through normalization.
sub setRatios	
{
 my($diheFunctHandle,$inputPDL,$atomNum,$atomTypes) = @_;
 my $energyGroupSum; my $scaleFactor;my %residueRatio;
 my $sum; my $diheStrengthTotal;
 undef %uniqueBonds;

 foreach my $res(keys %{$diheFunctHandle})
 {
		adjustFactorsHelper1($diheFunctHandle->{$res},$inputPDL->{$res},$atomNum,$atomTypes,\%residueRatio,\$sum);
 }

 foreach my $res(keys %{$diheFunctHandle})
 {
		adjustFactorsHelper2($diheFunctHandle->{$res},$inputPDL->{$res},$atomNum,$atomTypes,$diheStrengthTotal,\$sum);
 }
}


sub adjustFactorsHelper1
{
	my($diheArr,$inputPDL,$totalAtoms,$atomTypes,$residueRatio,$sum) = @_;
	my $totalStrength;my $normalize;
	my $contactTotal;my $relativeRatio;
 	my $size = $diheArr->dim(1);
 	for(my $i=0;$i<$size;$i++)
 	{
	if($size <6){last;}
	my($a,$b,$c,$d,$func,$count,$eG) = $diheArr->slice("0:6,$i:$i")->list;
	## Convert from relative index to absolute indexing ##
    	$a = sclr(slice($inputPDL,"3:3,$a,:"));
	$b = sclr(slice($inputPDL,"3:3,$b,:"));
	$c = sclr(slice($inputPDL,"3:3,$c,:"));
	$d = sclr(slice($inputPDL,"3:3,$d,:"));
    	my $resTypea = $atomTypes->{$a}->[1];
	my $resTypeb = $atomTypes->{$b}->[1];
	my $resTypec = $atomTypes->{$c}->[1];
    	my $resTyped = $atomTypes->{$d}->[1];
    	my $resTypeUse;

 	if(!$resTypeb || !$resTypec) {smog_quit ("Atom $b or $c doesn't have a residue type");}

 	## $eG == -1 is IMPROPER SKIP ##
	if($eG < 0) 
	{			
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
     smog_quit("Normalize option not set for $resTypea-$resTypeb-$resTypec-$resTyped of energyGroup $eG with atom indices $a-$b-$c-$d");
    }
    ## Normalize option is set ##	
	if($normalize eq 1) 
	{
		$relativeRatio=$termRatios->{$resTypeUse}->{"energyGroup"}->{$eG}->{"intraRelativeStrength"};
		if(!defined $relativeRatio){
			smog_quit("normalize=0, but intraRelativeStrength not defined for energyGroup $eG. Check .sif file");
		}

        ${$sum}+=($count*$relativeRatio);
		$count*=$relativeRatio;
        set($diheArr,5,$i,$count);	
	}

	
 	}
	
}

sub adjustFactorsHelper2
{
	my($diheArr,$inputPDL,$totalAtoms,$atomTypes,$diheStrengthTotal,$sum) = @_;
	my $totalStrength;my $normalize;
	my $contactTotal;my $diheRelative;
 	my $size = $diheArr->dim(1);
	my $totalDihedral;

for(my $i=0;$i<$size;$i++)
 {
	if($size <6){last;}
	my($a,$b,$c,$d,$func,$count,$eG) = $diheArr->slice("0:6,$i:$i")->list;
	## Convert from relative index to absolute indexing ##
	$a = sclr(slice($inputPDL,"3:3,$a,:"));
	$d = sclr(slice($inputPDL,"3:3,$d,:"));
    $b = sclr(slice($inputPDL,"3:3,$b,:"));
	$c = sclr(slice($inputPDL,"3:3,$c,:"));


    my $resTypea = $atomTypes->{$a}->[1];
	my $resTypeb = $atomTypes->{$b}->[1];
	my $resTypec = $atomTypes->{$c}->[1];
    my $resTyped = $atomTypes->{$d}->[1];
    my $resTypeUse;

 	if($eG < 0) {next;} ## IMPROPER NO NORMALIZATION --> Handled Internally as $eG==-1
	$eG = $eGTable{$eG}; ## Obtain user defined residue name ##
    if(!defined $termRatios->{$resTypeb}->{"energyGroup"}->{$eG})
    {
       ## CASE WHERE THERE IS A DIHEDRAL BETWEEN TWO DIFFERENT RESTYPES ##
       if(! defined $termRatios->{$resTypec}->{"energyGroup"}->{$eG})
       { 
           smog_quit("energyGroup $eG is not defined for $resTypeb-$resTypec\n");
       }
       $normalize = $termRatios->{$resTypec}->{"energyGroup"}->{$eG}->{"normalize"};
       $resTypeUse = $resTypec;
    }
    else{
        $normalize = $termRatios->{$resTypeb}->{"energyGroup"}->{$eG}->{"normalize"};
        $resTypeUse = $resTypeb;
    }
	
    
    ## Normalize option is set ##	
	if($normalize eq 1) 
	{
		my $diheLeftOver = 0;
                ## epsilonC+epsilonD ##
		$totalStrength = $termRatios->{"interRelativeTotal"};
		## epsilonC ##			
		$contactTotal = $termRatios->{"contactRelative"};
		## leftOver = totalAtoms*(1-epsilonC/(epsilonC+epsilonD)) ##			
		$diheLeftOver = $totalAtoms - $totalAtoms*($contactTotal/$totalStrength);
		$count = ($count/${$sum})*($diheLeftOver); 
		set($diheArr,5,$i,$count);
	}
    }
}

1;
