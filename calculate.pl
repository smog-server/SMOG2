#!/usr/bin/perl -w
#####################################################################
# calculate.pl: Do the calculation necessary for the values
# Author: Mohit Raghunathan													
# Date: May 2012															 
#####################################################################

#####################
## COMPILE HEADERS ##
#####################
use strict;
use warnings;
####################
## MODULE HEADERS ##
####################
use XML::Simple;
use Data::Dumper;
use pdbParser;
use PDL;
use mathFunctions;
use setRatios;
use Getopt::Long;
use Carp;
use XML::SAX::ParserFactory;
use XML::Validator::Schema; 
use IPC::Open3;


######################
## GLOBAL VARIABLES ##
######################

my $angToNano = 0.1;
my @diheArr;
my $OUTPUT;
my $help;my $inputPDB; my $inputFolder;
my $ignH;my $ignChainContact=0;my $ignHContacts;
my $topFile = "temp.top";
my $groFile = "temp.gro";
my $DCAFile = "";
my $shadowFile = "shadow.contacts";
my $ndxFile = "temp.ndx";
my $mapOpt = "shadow";my $contactOpt = "6 1";
my $numContacts = 0;
my $coarseFolder = ""; ## Coarse Grain Folder ##
my $contactRes = ""; ## Duplicate ##
my $softTerm=0; ## Boolean to for soft termination

#####################
## CALCULATE BONDS ##
#####################


##
# Combined bond calculation,sorting, and printing
sub printOrderedBonds
{

  my($bondFunctHandle1,$whichPDL1,$bondHandle2,$whichPDL2,$coarseGrained) = @_;
  my @bondCache; ## CACHE BONDS TO ORDER

  ## print directive headers ##
  print $OUTPUT "\n[ bonds ]\n";
  print $OUTPUT ";ai     aj      func    r0(nm)  Kb\n";
  ## Cache Bonds ##
  if($coarseGrained) {
  calculateBonds($bondFunctHandle1,$whichPDL1,\@bondCache);
  coarseCalculateBonds($bondHandle2,$whichPDL2,\@bondCache);
  }
  else
  {
	calculateBonds($bondFunctHandle1,$whichPDL1,\@bondCache);
  	connCalculateBonds($bondHandle2,$whichPDL2,\@bondCache);
  }

  ## Sort bonds by i then j ## 
  @bondCache = sort {($a->{"i"} <=> $b->{"i"}) || ($a->{"j"} <=> $b->{"j"}) } @bondCache;
  
  ## Print our sorted bonds ##
  foreach my $p(@bondCache)
  {print $OUTPUT $p->{"v"};}



}


##
# calculate bond angle of a specific residue  PDL
sub calculateBondsHelper
{
 my($bondHash,$inputPDL,$bondCache) = @_;
 my $bondArr = $bondHash->{"bonds"};
 my $bondFunc =$bondHash->{"functions"};
 my $bondIndex = 0;
 foreach my $bonds(@{$bondArr})
 {
	my ($atomOne,$atomTwo) = split("-",$bonds);
	my $outer = slice($inputPDL,"0:2,$atomOne,:")-slice($inputPDL,"0:2,$atomTwo,:");
	my @bondList = (sqrt(inner($outer,$outer))*$angToNano)->flat()->list();
	my @indexList1 = slice($inputPDL,"3:3,$atomOne,:")->flat()->list();
	my @indexList2 = slice($inputPDL,"3:3,$atomTwo,:")->flat()->list();
	for(my $i=0;$i<scalar(@bondList);$i++)
	{
		#print $OUTPUT bondOutput($bondFunc->[$bondIndex],$indexList1[$i],$indexList2[$i],$bondList[$i]);} 
  push @{$bondCache}, {'i' => $indexList1[$i],'j' => $indexList2[$i], 'v' => bondOutput($bondFunc->[$bondIndex],$indexList1[$i],$indexList2[$i],$bondList[$i])};
 }
		$bondIndex++;

 }
 
}

##
# Setup the calculation of bond angles of all the residues,
#	then call calculateBondsHelper() with a specific residue
sub calculateBonds
{
 my($bondFunctHandle,$whichPDL,$bondCache) = @_;
 #print $OUTPUT "\n[ bonds ]\n";
 #print $OUTPUT ";ai     aj      func    r0(nm)  Kb\n";
 foreach my $res(keys %{$bondFunctHandle})
 {
	if(!exists $whichPDL->{$res}){next;}
	calculateBondsHelper($bondFunctHandle->{$res},$whichPDL->{$res},$bondCache);
 }
}

##
#Setup the calculate of bond angles of all combined residues,
#then call connCalculateBondsHelper() with a specific residue
sub connCalculateBonds
{
 my($bondHandle,$whichPDL,$bondCache) = @_;
 foreach my $counter(keys %{$bondHandle})
 {
	connCalculateBondsHelper($bondHandle->{$counter},$whichPDL->{$counter},$bondCache);
 }
}

##
#Setup the calculate of bond angles of all combined residues,
#then call connCalculateBondsHelper() with a specific residue
sub coarseCalculateBonds
{
 my($bondHandle,$whichPDL,$bondCache) = @_;
 foreach my $counter(keys %{$bondHandle})
 {
	coarseCalculateBondsHelper($bondHandle->{$counter},$whichPDL->{$counter},$bondCache);

 }
}
##
#calculate bond angle of a specific combined residue PDL
sub coarseCalculateBondsHelper
{
 my($bondInfo,$inputPDL,$bondCache) = @_;

 my $size = $bondInfo->dim(1);

 for(my $i=0;$i<$size;$i++)
 {

 my $atomOne = sclr(slice($bondInfo,"0:0,$i:$i"));
 my $atomTwo = sclr(slice($bondInfo,"1:1,$i:$i"));
 my $func = sclr(slice($bondInfo,"2:2,$i:$i"));
 my $outer = slice($inputPDL,"0:2,$atomOne,:")-slice($inputPDL,"0:2,$atomTwo,:");
 my $output = sclr(sqrt(inner($outer,$outer))*$angToNano);
 $atomOne = sclr(slice($inputPDL,"3:3,$atomOne,:"));
 $atomTwo = sclr(slice($inputPDL,"3:3,$atomTwo,:"));
 
 #print $OUTPUT bondOutput(intToFunc("bonds",$func,""),$atomOne,$atomOne,$output);
 push @{$bondCache}, {'i' => $atomOne,'j' => $atomTwo, 'v' => bondOutput(intToFunc("bonds",$func,""),$atomOne, $atomTwo,$output)};
 }
 
}

##
#calculate bond angle of a specific combined residue PDL
sub connCalculateBondsHelper
{
 my($bondInfo,$inputPDL,$bondCache) = @_;
 my $atomOne = sclr(slice($bondInfo,"0:0"));
 my $atomTwo = sclr(slice($bondInfo,"1:1"));
 my $func = sclr(slice($bondInfo,"2:2"));
 my $outer = slice($inputPDL,"0:2,$atomOne,:")-slice($inputPDL,"0:2,$atomTwo,:");
 my $output = sclr(sqrt(inner($outer,$outer))*$angToNano);
 
 $atomOne = sclr(slice($inputPDL,"3:3,$atomOne,:"));
 $atomTwo = sclr(slice($inputPDL,"3:3,$atomTwo,:"));
 
 #print $OUTPUT bondOutput(intToFunc("bonds",$func,""),$atomOne,$atomOne,$output);
 push @{$bondCache}, {'i' => $atomOne,'j' => $atomTwo, 'v' => bondOutput(intToFunc("bonds",$func,""),$atomOne, $atomTwo,$output)};
 
}

##
# Convert internal bond function format to gromacs function format,
# internal format is f()+g()
sub bondOutput
{
 my($inputFunc,$i,$j,$values) = @_;
 if(!$inputFunc){confess "ERROR AT BONDOUTPUT::NO FUNCTION DEFINED\n";}
 ## Remove white spaces from string ##
 $inputFunc =~ s/^\s+|\s+$//g;
 ## If $input contains a + ##
 my @funcs = split(/\+/,$inputFunc); ## split combination of functions
 my $outputString="";
 my @paramArr;
 my $fType=""; my $formattedString;
	## Parse each part of the function ##
 foreach my $fun(@funcs)
 {
	my $paramFormat="";
	($fType,@paramArr) = bondFuncParser($fun);
     ## $paramArr[0] is bond angle
	if($paramArr[0] =~ /\?/){$paramArr[0]=$values;} ## non-native option used
	 
	 ## Check bonds threshold ##
	 if($paramArr[0] < $interactionThreshold->{"bonds"}->{"shortBond"}
	 || $paramArr[0] > $interactionThreshold->{"bonds"}->{"longBond"})
	 {confess("ERROR: BOND betwwen atom $i $j exceeded bonds threshold with value $paramArr[0]\n");}
	
	
 ## Format output <i j function p1 p2 ... pn>
 foreach(@paramArr){$paramFormat="$paramFormat %12.9e";}
	$formattedString = sprintf("%9d%9d%6d$paramFormat",$i,$j,$fType,@paramArr);
	$outputString = "$outputString$formattedString\n";
 }
 return $outputString;
 
 }


##
# Given a user defined function parse values
sub bondFuncParser
{
	my($bondFunc) = @_;
	my $fType;my @paramArr;
	my $funcName = ""; my $params = "";
 ($bondFunc =~ m/(.*)\((.*)\)/);
	$funcName = $1;
	$params = $2;
	$fType = returnFunction($funcName); ## Obtain Gromacs function type
 @paramArr = split(",",$params); ## Parse all parameters
 foreach $params(@paramArr){$params =~ s/^\s+|\s+$//g;}
	return ($fType,@paramArr);
}


######################
## CALCULATE ANGLES ##
######################


##
# Combined bond calculation,sorting, and printing
sub printOrderedAngles
{
		my($connAngleFunctionals,$connPDL) = @_;
  my @angleCache; ## CACHE BONDS TO ORDER

  ## print directive headers ##
  print $OUTPUT "\n[ angles ]\n";
  print $OUTPUT ";ai  aj   ak  func  th0(deg)   Ka\n";

  ## Cache Bonds ##
  connCalculateAngles($connAngleFunctionals,$connPDL,\@angleCache);

  ## Sort angles by i,j then k ## 
  @angleCache = sort {($a->{"i"} <=> $b->{"i"}) || ($a->{"j"} <=> $b->{"j"}) || ($a->{"k"} <=> $b->{"k"})} @angleCache;
  
  ## Print our sorted bonds ##
  foreach my $p(@angleCache)
  {print $OUTPUT $p->{"v"};}



}


##
# calculate 3 angle value of a specific residue PDL
sub calculateAnglesHelper
{
 my($angleArr,$inputPDL) = @_;
 
 foreach my $angles(@{$angleArr})
 {
	my ($atom1,$atom2,$atom3) = split("-",$angles);
 my $left = norm(slice($inputPDL,"0:2,$atom1,:")-slice($inputPDL,"0:2,$atom2,:"));
	my $right = norm(slice($inputPDL,"0:2,$atom3,:")-slice($inputPDL,"0:2,$atom2,:"));
	my $acos = inner($left,$right);
 
 }

}

##
# Setup the calculation of 3 angle value of all the residues,
# then call calculateAnglesHelper() with a specific residue
sub calculateAngles
{
 my($angleFunctHandle,$whichPDL) = @_;
 foreach my $res(keys %{$angleFunctHandle})
 {
	calculateAnglesHelper($angleFunctHandle->{$res}->{"angles"},$whichPDL->{$res});
 }
}

##
# Setup the calculate of bond angles of all combined residues,
# then call connCalculateAnglesHelper() with a specific residue
sub connCalculateAngles
{
 my($angleFunctHandle,$whichPDL,$angleCache) = @_;
 #print $OUTPUT "\n[ angles ]\n";
 #print $OUTPUT ";ai  aj   ak  func  th0(deg)   Ka\n";

 foreach my $res(keys %{$angleFunctHandle})
 {
	connCalculateAnglesHelper($angleFunctHandle->{$res},$whichPDL->{$res},$angleCache);
 }
}

##
# calculate bond angle of a specific combined residue PDL
sub connCalculateAnglesHelper
{
 my($angleArr,$inputPDL,$angleCache) = @_;
 my $size = $angleArr->dim(1);
 
 for(my $i=0;$i<$size;$i++)
 {
		my ($atom1,$atom2,$atom3,$func) = $angleArr->slice(":,$i:$i")->list;
  my $left = norm(slice($inputPDL,"0:2,$atom1,:")-slice($inputPDL,"0:2,$atom2,:"));
		my $right = norm(slice($inputPDL,"0:2,$atom3,:")-slice($inputPDL,"0:2,$atom2,:"));
		my $acos = mathFunctions::rad_to_deg(mathFunctions::acos(sclr(inner($left,$right))));
	
		my $atomi1 = sclr(slice($inputPDL,"3:3,$atom1,:"));
		my $atomi2 = sclr(slice($inputPDL,"3:3,$atom2,:"));
		my $atomi3 = sclr(slice($inputPDL,"3:3,$atom3,:"));

	
		#print $OUTPUT angleOutput(intToFunc("angles",$func,""),$atomi1,$atomi2,$atomi3,$acos);
  push @{$angleCache}, {'i' => $atomi1,'j' => $atomi2,'k' =>$atomi3, 'v' => angleOutput(intToFunc("angles",$func,""),$atomi1,$atomi2,$atomi3,$acos)};

 }

}

## 
# Convert internal angle function format to gromacs function format
sub angleOutput
{
 my($inputFunc,$i,$j,$k,$values) = @_;
	if(!$inputFunc){confess "ERROR AT ANGLEOUTPUT::NO FUNCTION DEFINED\n";}
 ## Remove white spaces from string ##
 $inputFunc =~ s/^\s+|\s+$//g;
 ## If $input contains a + ##
	my @funcs = split(/\+/,$inputFunc);
	my $outputString="";
	my $value="";my $kd="";
	my $fType="";my @paramArr;
	my $formattedString="";
	foreach my $fun(@funcs)
	{
	my $paramFormat="";
	($fType,@paramArr) = angleFuncParser($fun);
	if($paramArr[0] =~ /\?/){$paramArr[0]=$values;} ## non-native option used
	
	## Check angles threshold ##
	 if($paramArr[0] < $interactionThreshold->{"angles"}->{"smallAngles"})
	 	 {confess("ERROR: ANGLES between atom $i $j $k exceeded angles threshold with value $paramArr[0]\n");}
	
 ## Format output <i j k function p1 p2 ... pn>
 foreach(@paramArr){$paramFormat="$paramFormat %12.9e";}
	$formattedString = sprintf("%9d%9d%9d%6d$paramFormat",$i,$j,$k,$fType,@paramArr);
	$outputString = "$outputString$formattedString\n";
	}
	return $outputString; 
}

##
# Given a user defined angle function parse values
sub angleFuncParser
{
	my($angFunc) = @_;
	my $fType;my @paramArr;
	my $funcName = ""; my $params = "";
 ($angFunc =~ m/(.*)\((.*)\)/);
	$funcName = $1;
	$params = $2;
	$fType = returnFunction($funcName); ## Obtain Gromacs function type
 @paramArr = split(",",$params); ## Parse all parameters
 foreach $params(@paramArr){$params =~ s/^\s+|\s+$//g;}
	return ($fType,@paramArr);
}



#########################
## CALCULATE DIHEDRALS ##
#########################


##
# Combined dihedral calculation,sorting, and printing
sub printOrderedDihedrals
{

		my($connDiheFunctionals,$connPDL) = @_;
  my @diheCache; ## CACHE BONDS TO ORDER

  ## print directive headers ##
		print $OUTPUT "\n[ dihedrals ]\n";
 	print $OUTPUT ";ai  aj  ak  al  func  phi0(deg) kd mult\n";

  ## Cache Bonds ##
  connCalculateDihedrals($connDiheFunctionals,$connPDL,\@diheCache);

  ## Sort dihedrals by i,j then k ## 
  @diheCache = sort {($a->{"j"} <=> $b->{"j"}) || ($a->{"k"} <=> $b->{"k"}) || ($a->{"i"} <=> $b->{"i"})} @diheCache;
  
  ## Print our sorted bonds ##
  foreach my $p(@diheCache)
  {print $OUTPUT $p->{"v"};}

  print $OUTPUT "\n";

}


##
# calculate dihedral angle of a specific residue PDL
sub calculateDihedralsHelper
{
 my($diheArr,$inputPDL) = @_;
 foreach my $dihedrals(@{$diheArr})
 {
	my ($atom1,$atom2,$atom3,$atom4) = split("-",$dihedrals);
	
 my $b1 = slice($inputPDL,"0:2,$atom2,:")-slice($inputPDL,"0:2,$atom1,:");
	my $b2 = slice($inputPDL,"0:2,$atom3,:")-slice($inputPDL,"0:2,$atom2,:");
	my $b3 = slice($inputPDL,"0:2,$atom4,:")-slice($inputPDL,"0:2,$atom3,:");
	
	my $n1 = norm(crossp($b1,$b2));
	my $n2 = norm(crossp($b2,$b3));
	my $m1 = crossp($n1,norm($b2));
	
	my $cos = atan2(inner($m1,$n2),inner($n1,$n2));
	my $negcos = $cos; #$cos->where($cos < 0);

	$negcos*=-1;$negcos+=mathFunctions::pi();
	my $acos = sclr(mathFunctions::rad_to_deg($cos));
	
	my $atomi1 = sclr(slice($inputPDL,"3:3,$atom1,:"));
	my $atomi2 = sclr(slice($inputPDL,"3:3,$atom2,:"));
	my $atomi3 = sclr(slice($inputPDL,"3:3,$atom3,:"));
	my $atomi4 = sclr(slice($inputPDL,"3:3,$atom4,:"));
 }
}

##
# Setup the calculate of dihedral angles of all the residues,
# then call calculateDihedralsHelper() with a specific residue
sub calculateDihedrals
{
 my($diheFunctHandle,$whichPDL) = @_;
 foreach my $res(keys %{$whichPDL})
 {
	calculateDihedralsHelper($diheFunctHandle->{$res}->{"dihedrals"},$whichPDL->{$res});
 }
}



##
# Setup the calculate of bond angles of all combined residues,
# then call connCalculateBondsHelper() with a specific residue
sub connCalculateDihedrals
{
 my($diheFunctHandle,$whichPDL,$diheCache) = @_;
 #print $OUTPUT "\n[ dihedrals ]\n";
 #print $OUTPUT ";ai  aj  ak  al  func  phi0(deg) kd mult\n";
 foreach my $res(keys %{$diheFunctHandle})
 {
		connCalculateDihedralsHelper($diheFunctHandle->{$res},$whichPDL->{$res},$diheCache);
 }
}

sub connCalculateDihedralsHelper
{
 my($diheArr,$inputPDL,$diheCache) = @_;
 my $size = $diheArr->dim(1);
 for(my $i=0;$i<$size;$i++)
 {
		my ($atom1,$atom2,$atom3,$atom4,$func,$cD,$eG) = $diheArr->slice(":,$i:$i")->list;
	
	   
        my $b1 = slice($inputPDL,"0:2,$atom2,:")-slice($inputPDL,"0:2,$atom1,:");
		my $b2 = slice($inputPDL,"0:2,$atom3,:")-slice($inputPDL,"0:2,$atom2,:");
		my $b3 = slice($inputPDL,"0:2,$atom4,:")-slice($inputPDL,"0:2,$atom3,:");
	
		my $n1 = norm(crossp($b1,$b2));
		my $n2 = norm(crossp($b2,$b3));
		my $m1 = crossp($n1,norm($b2));
	
		my $cos = atan2(inner($m1,$n2),inner($n1,$n2));
		my $negcos = $cos;#$cos->where($cos < 0);

		$negcos*=-1;$negcos+=mathFunctions::pi();
		my $acos = sclr(mathFunctions::rad_to_deg($cos));
	
		my $atomi1 = sclr(slice($inputPDL,"3:3,$atom1,:"));
		my $atomi2 = sclr(slice($inputPDL,"3:3,$atom2,:"));
		my $atomi3 = sclr(slice($inputPDL,"3:3,$atom3,:"));
		my $atomi4 = sclr(slice($inputPDL,"3:3,$atom4,:"));

		## $eG >=0 dihedral is proper dihedral of any functional format
        if($eG >= 0)
		{
	 	    #print $OUTPUT dihedralOutput(intToFunc("dihedrals",$func,$eG),$atomi1,$atomi2,$atomi3,$atomi4,$acos,$cD);
			my $val = dihedralOutput(intToFunc("dihedrals",$func,$eG),$atomi1,$atomi2,$atomi3,$atomi4,$acos,$cD,$eG);
			if(!$val) {next;}
			push @{$diheCache}, {'i' => $atomi1,'j' => $atomi2,'k' =>$atomi3, 'l' => $atomi4, 'v' => dihedralOutput(intToFunc("dihedrals",$func,$eG),$atomi1,$atomi2,$atomi3,$atomi4,$acos,$cD,$eG)};
		}
		## Improper dihedral
		else
		{
			#print $OUTPUT dihedralOutput(intToFunc("impropers",$func,$eG),$atomi1,$atomi2,$atomi3,$atomi4,$acos,$cD);}
			push @{$diheCache}, {'i' => $atomi1,'j' => $atomi2,'k' =>$atomi3, 'l' => $atomi4, 'v' => dihedralOutput(intToFunc("impropers",$func,$eG),$atomi1,$atomi2,$atomi3,$atomi4,$acos,$cD,-1)};
 	   }

 }

}

## Convert .stop dihedral function format to gromacs function format
sub dihedralOutput
{
 my($inputFunc,$ai,$aj,$ak,$al,$values,$ratio,$eG) = @_;
 ## Remove white spaces from string ##
 $inputFunc =~ s/^\s+|\s+$//g;
 ## If $input contains a + ##
	my @funcs = split(/\+/,$inputFunc);
	my $outputString="";
	my $fType="";my $exportRatio;       
 	my $formattedString="";
	my $valueSend;my @paramArr;
	foreach my $fun(@funcs)
	{
	my $paramFormat="";
	($fType,@paramArr) = diheFuncParser($fun);
	
	## 0 DIHEDRAL ##
	if($fType == 0){return "";}
	## If dihedral angle is to be native and/or has scaling ##
	if($paramArr[0] =~ /^\?$/){$paramArr[0] = $values;}
 	elsif($paramArr[0] =~ /\?\*(.*)/){$paramArr[0]=$values*eval($1);}

	## If dihedral kd is normalized or has scaling&normalize ##
	if($paramArr[1] =~ /^\?$/){$paramArr[1]=$ratio;} ##? has to be consistent with normalization##
 	elsif($paramArr[1] =~ /\?\*(.*)/){$paramArr[1]=eval($1)*$ratio;}
 	else{$paramArr[1]*=$ratio;} ## Simple scaling ##
 
	### MINOR KLUDGE TO FIX RIGID/IMPROPER DIHEDRALS BEING OFF BY +-180 ###
	if($fType == 2)
	{
		if($paramArr[0]>180){$paramArr[0]-=180;}else{$paramArr[0]+=180;}
	}
 
	## Format output <i j k function p1 p2 ... pn>
 	foreach(@paramArr){$paramFormat="$paramFormat %12.9e";}

 ## MINOR KLUDGE TO FIX MULTIPLICITY FORMATTING FOR TYPE 1 ##
 ## BLAME PAUL##
 if($fType eq 1) {$paramFormat=" %12.9e	%12.9e %u";}

	$formattedString = sprintf("%9d%9d%9d%9d%3d$paramFormat",$ai,$aj,$ak,$al,$fType,@paramArr);
	$outputString = "$outputString$formattedString\n";
	}
	return $outputString;
 
 }
 
##
# Given a user defined dihedral function parse values
sub diheFuncParser
{
	my($diheFunc) = @_;
 my $fType;my @paramArr;
	my $funcName = "";my $params = "";
	($diheFunc =~ m/(.*)\((.*)\)/);
	$funcName = $1; $params = $2;

	$fType = returnFunction($funcName);
	@paramArr = split(",",$params); ## Parse all parameters;
	foreach $params(@paramArr){$params =~ s/^\s+|\s+$//g;}
	return ($fType,@paramArr);
}



#########################
## CALCULATE 1-4 PAIRS ##
#########################
## NOTE: 1-4 ARE ASSUMED TO BE NON NATIVE ##




########################
## CALCULATE CONTACTS ##
########################

sub calculateContacts
{
	my($contactPDL,$inputPDL,$atomTypes,$numCon,$numAtom) = @_;
	my $atoma; my $atomb;my $fType; my $dist;
	my $nbtypea; my $nbtypeb;
	my $resTypea; my $resTypeb;
	my $contactInterScale;
 	my $dihedralInterScale;
	my $multFactor = 1;
	my $c0; my $c1; 
	my $epsilon;my $cG;my $funct;
 	my $totalStrength; my $outputString;
	my $contactIntraScale;
	my $contactIntraTotal;	
	my $count;
	
	print $OUTPUT "\n\n[ pairs ]\n";
	
	## REMOVE ##
	my $totalEpsilon=0;

	## Sum all initial contact strengths ##	
	for(my $i=0;$i<$numCon;$i++)
	{
			my $normalize;
			$atoma = sclr(slice($contactPDL,"1:1,$i:$i"));
	  		$atomb = sclr(slice($contactPDL,"3:3,$i:$i"));
			$nbtypea = $atomTypes->{$atoma}->[0];
			$nbtypeb = $atomTypes->{$atomb}->[0];
			$resTypea = $atomTypes->{$atoma}->[1];
			$resTypeb = $atomTypes->{$atomb}->[1];
			($funct,$cG) = getContactFunctionals($nbtypea,$nbtypeb);
			if(!$funct || !$cG)
			{confess("\n No contact Function defined for nbType contacts $nbtypea-$nbtypeb\n");}
			$normalize = $termRatios->{"contactGroup"}->{$cG}->{"normalize"};
			if(!$normalize){next;}		
			$contactIntraTotal = $termRatios->{"cintraRelativeTotal"};
			$contactIntraScale = $termRatios->{"contactGroup"}->{$cG}->{"intraRelativeStrength"};
			$count+=($contactIntraScale/$contactIntraTotal);
	}

	## Adjust ratio ##
	for(my $i=0;$i<$numCon;$i++)
	{
			my $contactLeftOver;my $diheTotal; my $normalize;
			my $deltaMin=-1;my $scale=-1;my $deltaMax=-1;
			my $resIdxA=0;my $resIdxB=0;
			my $atomA="";my $atomB="";
			
	  		$atoma = sclr(slice($contactPDL,"1:1,$i:$i"));
	  		$atomb = sclr(slice($contactPDL,"3:3,$i:$i"));
			$nbtypea = $atomTypes->{$atoma}->[0];
			$nbtypeb = $atomTypes->{$atomb}->[0];
			$resTypea = $atomTypes->{$atoma}->[1];
			$resTypeb = $atomTypes->{$atomb}->[1];
			($funct,$cG) = getContactFunctionals($nbtypea,$nbtypeb);
			if(!$funct || !$cG)
			{confess("\n No contact Function defined for nbType contacts $nbtypea-$nbtypeb\n");} 
			
			$normalize = $termRatios->{"contactGroup"}->{$cG}->{"normalize"};

            ## Stacking scaling ##
            $epsilon = 1;
            if(exists $contactSettings->{"contactScaling"}->{$resTypea}
            && exists $contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb})
            {
              $deltaMin=$contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb}->{"deltaMin"};
              $deltaMax=$contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb}->{"deltaMax"};

              $scale = $contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb}->{"scale"};
              $scale = eval($scale);
              $resIdxA = $allAtoms{$atoma}->[2];
              $resIdxB = $allAtoms{$atomb}->[2];
              
              $atomA = $allAtoms{$atoma}->[3];
              $atomB = $allAtoms{$atomb}->[3];
              ##Atom to boolean##
              $atomA = exists($contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb}->{"atomList"}->{$atomA});
              $atomB = exists($contactSettings->{"contactScaling"}->{$resTypea}->{$resTypeb}->{"atomList"}->{$atomB});

              
              if(abs($resIdxB-$resIdxA)>=$deltaMin 
              && abs($resIdxB-$resIdxA)<=$deltaMax
              && $atomA && $atomB)
              {
               $epsilon*=$scale;              
              }
              else{$epsilon=1;}
            
            }

			## Non-Native option for contacts ##
			if(!$normalize) {
			$dist = sclr(slice($contactPDL,"4:4,$i:$i"));
			$c0=$dist;
	  		$c1=$dist;
			print $OUTPUT contactOutput($funct,$atoma,$atomb,$c0,$c1,$epsilon),"\n";
			}	
			
			else {	
			## DIHE TO CONTACT SCALING ##
			## DIHEDRAL TO CONTACT RATIO IS GLOBAL IRREGARDLESS OF RESIDE TYPE ##
			$totalStrength = $termRatios->{"interRelativeTotal"};
			$diheTotal = $termRatios->{"energyRelative"};
			$contactLeftOver = $numAtom - $numAtom*($diheTotal/$totalStrength);
			$contactLeftOver = $contactLeftOver/$count;
   			$multFactor = $contactLeftOver;

			## CONTACT TO CONTACT SCALING ##			
			$contactIntraTotal = $termRatios->{"cintraRelativeTotal"};
			$contactIntraScale = $termRatios->{"contactGroup"}->{$cG}->{"intraRelativeStrength"};
			$multFactor = $multFactor*($contactIntraScale/$contactIntraTotal);

			$epsilon = $epsilon*$multFactor;
	  		$dist = sclr(slice($contactPDL,"4:4,$i:$i"));
	  		$c0=$dist;
	  		$c1=$dist;
	  		$totalEpsilon+=$epsilon;
			print $OUTPUT contactOutput($funct,$atoma,$atomb,$c0,$c1,$epsilon),"\n";
			}
	}
	
	#print "TOTAL EPSIOLON IS $totalEpsilon\n";
}

##
# Convert internal contact function format to gromacs function format,
# internal format is f()+g()
sub contactOutput
{
 my($inputFunc,$i,$j,$c0,$c1,$epsilon) = @_;
 if(!$inputFunc){confess "ERROR AT CONTACTOUTPUT::NO FUNCTION DEFINED\n";}
 ## Remove white spaces from string ##
 $inputFunc =~ s/^\s+|\s+$//g;
 ## If $input contains a + ##
 my @funcs = split(/\+/,$inputFunc); ## split combination of functions
 my $outputString="";
 my @paramArr;
 my $fType=""; my $formattedString;
 ## Parse each part of the function ##
 foreach my $fun(@funcs)
 {
	my $paramFormat="";
	($fType,@paramArr) = contactFuncParser($fun);

	## LJ CONTACTS 6-12 ##	
	if($fType == 1){contactParseLJ612(\@paramArr,$epsilon,$c0,$c1);$fType=1;}
        ## LJ CONTACTS 10-12 ##
    elsif($fType == 2){contactParseLJ1012(\@paramArr,$epsilon,$c0,$c1);$fType=1;}
        ## Gaussian ##
    elsif($fType == 3){contactParseGaussian(\@paramArr,$epsilon,$c0);}
    else{confess("Contact function type $fType is not defined\n");}

    ## PAIRS Ftype ALWAYS 1 unless gaussian 6 ##
    if($fType==3){$fType = 6;}
 	## Format output <i j function p1 p2 ... pn>
 	foreach(@paramArr){$paramFormat="$paramFormat	%12.9e";}
	$formattedString = sprintf("%d %d %d$paramFormat",$i,$j,$fType,@paramArr);
	$outputString = "$outputString$formattedString";
 }
 return $outputString;
 
}

##
# Parse LJ 6-12 Params
sub contactParseLJ612
{
  my($paramArr,$epsilon,$c6,$c12) = @_;
  ##i=0 is epsilon##
  ##i=1 is c6##
  ##i=2 is c12##
  
  if(($paramArr->[1] =~ /^\?$/ xor $paramArr->[2] =~ /^\?$/))
  {confess("Contact function ill defined\n");}

  ## Exponentiate c6 c12 ##
  $c6=$c6**6;$c12=$c12**12;

  if($paramArr->[1] =~ /^\?$/){$paramArr->[1]=$c6;}
  if($paramArr->[2] =~ /^\?$/){$paramArr->[2]=$c12;}
   

  ## native option used
  if($paramArr->[0] =~ /^\?$/)
  {
    $paramArr->[1]=2*$epsilon*$c6;
    $paramArr->[2]=$epsilon*$c12;
  }
  ## non-native option used
  else{$paramArr->[1]=2*$paramArr->[0]*$c6*$epsilon;$paramArr->[2]=$paramArr->[0]*$c12*$epsilon;} 
  ## Shift first value off ##
  shift @{$paramArr};

}


##
# Parse LJ 6-12 Params
sub contactParseLJ1012
{
  my($paramArr,$epsilon,$c10,$c12) = @_;
  ##i=0 is epsilon##
  ##i=1 is c6##
  ##i=2 is c12##

  if(($paramArr->[1] =~ /^\?$/ xor $paramArr->[2] =~ /^\?$/))
  {confess("Contact function ill-defined\n");}
  
   ## Exponentiate c10 c12 ##
  $c10=$c10**10;$c12=$c12**12;

  if($paramArr->[1] =~ /^\?$/){$paramArr->[1]=$c10;}
  if($paramArr->[2] =~ /^\?$/){$paramArr->[2]=$c12;}
  

  ## native option used
  if($paramArr->[0] =~ /^\?$/)
  {
    $paramArr->[1]=6*$epsilon*$c10;
    $paramArr->[2]=5*$epsilon*$c12;
  }
  ## non-native option used
  else{$paramArr->[1]=6*$paramArr->[0]*$c10*$epsilon;$paramArr->[2]=5*$paramArr->[0]*$c12*$epsilon;} 
  ## Shift first value off ##
  shift @{$paramArr};

}

##
# Parse Gaussian Params
sub contactParseGaussian
{

 my($paramArr,$epsilon,$r0) = @_;
 my $A=0;my $sigma_G=0;my $a=0;
  ##i=0 is epsilon_c##
  ##i=1 is epsilon_nc#
  ##i=2 is sigma_gaussian##
  ##i=3 is r ##


  if($paramArr->[1] =~ /\?/ && !($paramArr->[2] =~ /^\?$/))
  {confess("ERROR AT CONTACT FUNCTION PARSING:: CONFLICTING INPUTS\n");}

  ## Epsilon_c ##
  if($paramArr->[0] =~ /^\?$/){$paramArr->[0]=$epsilon;}
  else{$paramArr->[0]=$epsilon*($paramArr->[0]);}
  ## Epsilon_nc ##
  if($paramArr->[1] =~ /\?/){confess("ERROR AT CONTACT FUNCTION PARSING:: CONFLICTING INPUTS\n");}

  ## Check if sigma_G has r0 dependency ##
  if($paramArr->[2] =~ /\?/)
  {
     #Replace ? with $r0 
     $paramArr->[2] =~ s/\?/$r0/g;
     #Eval expression
     $paramArr->[2] = eval($paramArr->[2]);
  }
  if($paramArr->[3] =~ /^?$/){$paramArr->[3]=$r0;}
  
  ($A,$r0,$sigma_G,$a) = ($paramArr->[0],$paramArr->[3],$paramArr->[2],$paramArr->[1]);
  @{$paramArr}=($A,$r0,$sigma_G,$a);

}


##
# Given a user defined contact function parse values
sub contactFuncParser
{
	my($contactFunc) = @_;
	my $fType;my @paramArr;
	my $funcName = ""; my $params = "";
 	#($contactFunc =~ m/(.*)\((.*)\)/);
 	($contactFunc =~ /^([^(]+)\((.*)\)$/);
 	$funcName = $1;
	$params = $2;
	$fType = returnFunction($funcName); ## Obtain Gromacs function type
 	@paramArr = split(",",$params); ## Parse all parameters
 	foreach $params(@paramArr){$params =~ s/^\s+|\s+$//g;}
	return ($fType,@paramArr);
}



##########################
## CALCULATE EXCLUSIONS ##
##########################


sub calculateExclusions
{
  my($fileName,$noChainFlag) = @_;
  my $numContacts = 0; my $garbage = 0;
  my $line = "";
  my $type1;my $type2; my $contact1; my $contact2;
  my $dist;
  
  ## OPEN .contact FILE ##
  unless (open(MYFILE, $fileName)) {
    print "Cannot read from '$fileName'.\nProgram closing.\n";
    <STDIN>;
    exit;
  }
 
  print $OUTPUT "\n[ exclusions ]\n";
  print $OUTPUT ";ai	aj\n";
  
  while($line = <MYFILE>)
  {
	($contact1,$type1,$contact2,$type2,$dist) = split(/\s+/,$line);
	## INTER CHAIN CONTACT IGNORE ##
 	if($noChainFlag && ($allAtoms{$type1}->[2] ne $allAtoms{$type2}->[2])){next;}
	print $OUTPUT $type1,"	",$type2,"\n";
  }



}

sub calculateAtomTypes
{
	my($uniqeHandle) = @_;
	my $formattedString = "";	
	print $OUTPUT "\n [ atomtypes ] \n";
	print $OUTPUT ";name  mass     charge   ptype c6       c12\n";
	
	## OBTAIN UNIQUE ATOM TYPES ##
	my %obtained;my $nbtype ="";
	my $c12;my $c6;my $ptype;my $charge;my $mass;
	
	foreach my $atoms(keys %{$uniqeHandle})
	{
		$nbtype = $uniqeHandle->{$atoms}->[0];
  		if(!$uniqeHandle->{$atoms}->[0]){print Dumper $uniqeHandle->{48}; exit;}

		if(exists $obtained{$nbtype}){next;}
		$c6 = $interactions->{"nonbonds"}->{$nbtype}->{"c6"};
		$c12 = $interactions->{"nonbonds"}->{$nbtype}->{"c12"};
		$ptype = $interactions->{"nonbonds"}->{$nbtype}->{"ptype"};
		$mass = $interactions->{"nonbonds"}->{$nbtype}->{"mass"};
		$charge = $interactions->{"nonbonds"}->{$nbtype}->{"charge"};
		if(!$c6 || !$c12 || !$ptype || !$mass || !$charge)
		{confess("Nonbond Params at atom types not set for nbType $nbtype\n");}
		
		$formattedString = sprintf("%6s  $mass  $charge  $ptype %13.5e%13.5e",$nbtype,$c6,$c12);
		print $OUTPUT $formattedString,"\n";
		$obtained{$nbtype}=1;
	}
}

sub calculateAtoms
{
	my($inputPDB) = @_;
	open(PDB,$inputPDB);
	my $line = "";
	print $OUTPUT "\n[ atoms ]\n";	
	print $OUTPUT ";nr  type  resnr residue atom  cgnr\n";
	my $atomName;my $atomNum; my $resName; my $resNum;
	my $charge = 0; my $mass = 1; my $atomType="none";
	my $counter=0; my $resCounter=1;my $resNumPrev=1;
	while($line = <PDB>)
	{
		if($line =~ m/END/) {last;}
		if($line =~ m/ATOM/ || $line =~ m/HETATM/)
		{
			$atomName = substr($line, 12, 4);
			$atomName =~ s/^\s+|\s+$//g;
			#$atomNum = substr($line,6,5);
			
			$resName = substr($line,17,4);
			$resName =~ s/^\s+|\s+$//g;
			$resNum = substr($line,22,4);
			$resNum =~ s/^\s+|\s+$//g;
   			if($resNum ne $resNumPrev)
			{
				$resNumPrev=$resNum;
				$resCounter++;
			}
			if(!$residues{$resName}->{"atoms"}->{$atomName}){next;}
			$counter++;$atomNum=$counter;
			$atomType = $residues{$resName}->{"atoms"}->{$atomName}->{"nbType"};
			printf $OUTPUT "%6d %10s %6d %6s %6s %6d\n", $atomNum, $atomType, $resCounter, $resName, $atomName, $atomNum;
			

		}
	}


}




sub printDefaultParamsHeader
{
	print $OUTPUT "; SBMV2 FOR GROMACS/NAMD\n";
	print $OUTPUT "\n[ defaults ]\n";
	print $OUTPUT ";nbfunc comb-rule gen-pairs\n";
	print $OUTPUT "           1           1 no\n";

}

sub printMoleculeTypes
{

	print $OUTPUT "\n[ moleculetype ]\n";
	print $OUTPUT ";name   nrexcl\n";
	print $OUTPUT "Macromolecule           3\n";

}

sub printDefaultParamsFooter
{
	print $OUTPUT "\n[ system ]\n";
 	print $OUTPUT ";name\n";
 	print $OUTPUT "Macromolecule\n";

 	print $OUTPUT "\n[ molecules ]\n";
 	print $OUTPUT ";name   #molec\n";
 	print $OUTPUT "Macromolecule   1\n";



}



##
# Convert PDB to GRO
# Store chain info to ndx file
sub convertPDBToGroNdx
{
	my ($inputPDB,$outputGRO,$outputNDX) = @_;
	my $counter = 0 ;	
	my $line = "";
	my $resName="";
	my $resNum="";
	my $atomName="";
	my $atomNum="";
 	my $x="";my $y="";my $z="";
 	my $output = "";
 	my $chain = "";
 	my $chainCounter = 1;
 	my %chainHash;
	my $resNumPrev=1;my $resCounter=1;
 	open(PDB,$inputPDB); ## INPUT PDB
	open(GRO,">tempGro");
	while($line=<PDB>)
	{
		
	 if($line =~ m/END/) {last;}
  	 if($line =~ m/TER/) {$chainCounter++;next;}
         if($line =~ m/ATOM|HETATM/)
         {
		$atomName = substr($line, 12, 4);
		$atomName =~ s/^\s+|\s+$//g;
		$resName = substr($line,17,4);
		$resName =~ s/^\s+|\s+$//g;
        if(!exists $residues{$resName}->{"atoms"}->{$atomName})
		{next;}
		
		$x = substr($line, 30, 8);$x*=0.10;
		$y = substr($line, 38, 8);$y*=0.10;
		$z = substr($line, 46, 8);$z*=0.10;
		
	
	 	$chain = substr($line,21,1);
  		if($chain eq " ") {$chain = $chainCounter;}
		$counter++;$atomNum=$counter;
		$resNum = substr($line,22,4);
		$resNum =~ s/^\s+|\s+$//g;
		if($atomNum == 1){$resNumPrev=$resNum;}
		if($resNum ne $resNumPrev)
			{
				$resNumPrev=$resNum;
				$resCounter++;
			}

  		## Save atom index to chain ##
  		push(@{$chainHash{$chain}},$atomNum);
		if($atomNum >= 100000){$atomNum = $atomNum%100000;}
	    if($resCounter >= 100000){$resCounter = $resCounter%100000;}
		$output = sprintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",$resCounter,$resName,$atomName,$atomNum,$x,$y,$z);
		print GRO $output;
		
	 }
	 	


	}

 ## ADJUST GRO FILE ##
	close(PDB);
	close(GRO);open(GROIN,"tempGro");
	open(GRO,">$outputGRO");
	print GRO "SMOGV2\n";
	print GRO $counter,"\n";
	while(<GROIN>) { print GRO $_;}
	close(GRO); close(GROIN);
	unlink("tempGro");
	
 ## CREATE NDX FILE ##
 open(NDX,">$outputNDX");
 my $c;my $v;
 while(($c,$v) =  each %chainHash)
 {
			print NDX "[ $c ]\n";
   print NDX join("\n",@{$v});
   print NDX "\n";
 }
 close(NDX);
 


}

sub validateTemplate
{
 my ($file,$type) = @_;
 my $validator = XML::Validator::Schema->new(file => "$ENV{SMOG_PATH}/schemas/$type.xsd");
 my $parser = XML::SAX::ParserFactory->parser(Handler => $validator);
 eval { $parser->parse_uri($file) };
 die "\nFailed at validating $file: $@" if $@;

}


sub parseInputFolder
{
 my ($folderName) = @_;
 opendir(my $folder,$folderName);
 my $bif;my $sif;my $b;my $nb;
 while(my $file = readdir($folder))
 {
   if($file =~ m/\.bif/)
   {$bif = $file;$bif = "$folderName/$bif";validateTemplate($bif,"bif");next;}
   if($file =~ m/\.sif/)
   {$sif = $file;$sif = "$folderName/$sif";validateTemplate($sif,"sif");next;}
   if($file =~ m/\.b/)
   {$b = $file;$b = "$folderName/$b";validateTemplate($b,"b");next;}
   if($file =~ m/\.nb/)
   {$nb = $file;$nb = "$folderName/$nb";validateTemplate($nb,"nb");next;}
  }

 setInputFileName($bif,$sif,$b,$nb);

}


sub usage
{
  print "\n***************************\n";
  print "********** SMOG2 **********";
  print "\n***************************\n";
  print "Unknown option: @_\n" if ( @_ );
  print "usage: smogv2 [-o file.top] [-i input.pdb] [-g file.gro] [--help|-?] [-t templateFolder/] [-c extra.contacts] [-ignh]\n";
		print "Necessary Inputs\n";		
		print "-i [input.pdb]: Input pdb to generate Hamiltonian\n";
		print "-tAA [templateFolder]: Folder containing templates of molecular and interaction definitions (All Atom)\n";
                print "-tCG [templateFolder]: Folder containing templates of molecular and interaction definitions (Coarse Grain)\n";
                print "-contactRes [templateFolder]: Folder containing templates for contact resolution\n";
		print "Optional Inputs\n";
		print "-g [file.gro]: Output .gro file name\n";
		print "-o [file.top]: Output .top file name\n";
        print "-s [file.shadow]: Output .shadow file name\n";
		print "-n [file.ndx]: Output .ndx file name\n";
		print "-c [extra.contacts]: Input contact map file\n";
		print "-ignH: Ignore hydrogens in PDB\n";
		print "-ignHContacts: Ignore hydrogen contacts\n";
		print "-ignChainContacts: Ignore contacts between chains\n";
		print printCitation();
  exit;
}

sub stripHydrogen
{
		my($inputFile) = @_;
		system("perl $ENV{SMOG_PATH}/tools/stripH.pl $inputFile");
		printf("Stripping Hydrogens off $inputFile. Creating ignH_$inputFile\n");
		$inputPDB = "$inputFile.new";
		printf("Generating Hamiltonian with $inputPDB\n");
		<STDIN>
		#printCitation();
}


sub setContactParams
{
			my $method = $contactSettings->{"method"};
			my $ignoreHydrogens="";
			if($ignHContacts){$ignoreHydrogens="--ignoreH";}
			if($method =~ m/shadow/)
			{
					my $radius = $contactSettings->{"shadowRadius"};
					my $dist = $contactSettings->{"contactDistance"};
				    return "-m shadow -c $dist -s $radius --distance $ignoreHydrogens";
			}
			elsif($method =~ m/cutoff/) 
			{
			  my $dist = $contactSettings->{"contactDistance"};
			  return "-m cutoff -c $dist -rd 0 --distance $ignoreHydrogens";
			
			}
			else {confess "Error at contact map method settings. $method is not a supported contactmap method.";}
}


sub printCitation
{
print "\n************************\n";
print "******* CITATION *******";
print "\n************************\n";
print 'SMOGV2 & Shadow is part of  http://smog-server.org
Direct questions to: info@smog-server.org
Work utilizing SMOGV2 should cite:

Noel JK, Whitford PC, Sanbonmatsu KY & Onuchic JN (2010)
SMOG@ctbp: simplified deployment of structure-based models in GROMACS.
Nucleic Acids Research, 38, W657‚Äì61. DOI: 10.1093/nar/gkq498

Work using the Shadow contact map should cite:

Noel, JK, Whitford PC & Onuchic JN (2012)
The shadow map: a general contact definition for capturing the
dynamics of biomolecular folding and function.
Journal of Physical Chemistry B 116, 8692‚Äì8702.
';
print "***************************************************\n";


}

##
## MEMORY MANAGEMENT ROUTINES
##
sub determine_optimal_memory {

    # We'll set a ceiling for the memory allocation.  On a 32-bit OS this is going
    # to be 1500m (the max it can safely handle), on a 64-bit OS we won't take more
    # than 6GB
    my $max_memory = 1500;

    # We need not only a 64 bit OS but 64 bit java as well. It's easiest to just test
    # java since the OS support must be there if you have a 64 bit JRE.

    my ($in,$out);
    open3(\*IN,\*OUT,\*OUT,"java -version") or print_error("Can't find java");
    close IN;
    while (<OUT>) {
        if (/64-Bit/) {
            $max_memory = 6000;
        }
    }
    close OUT;
    # The way we determine the amount of physical memory is OS dependent.
    my $os = $^O;

    my $physical;
    if ($os =~ /Win/) {
        confess("SMOG2 does not run on windows\n");
    }
    elsif ($os =~/darwin/ or $os =~ /bsd/i) {
        $physical = get_osx_memory($max_memory);
    }
    else {
        $physical = get_linux_memory($max_memory);
    }

    warn "Raw physical free memory is $physical\n";

    # We then set the memory to be the minimum of 2/3 of the physical
    # memory or the ceiling, whichever is lower.
    $physical = int(($physical/3)*2);

    if ($max_memory < $physical) {
        return $max_memory;
    }

    warn "Using $physical MB of RAM to launch seqmonk\n";
    return $physical;

}

sub get_linux_memory {
    # We get the amount of physical memory on linux by parsing the output of free

    open (MEM,"free -m |") or print_error("Can't launch free on linux: $!");

    while (<MEM>) {
        if (/^Mem:\s+(\d+)/) {
            return $1;
        }
    }

    close MEM;

    confess("Couldn't parse physical memory from the output of free");
}

sub get_osx_memory {

    # We get the amount of physical memory on OSX by parsing the output of top

    open (MEM,"top -l 1 -n 0 |") or print_error("Can't get amount of memory on OSX: $!");

    my $total_mem = 0;

    while (<MEM>) {
        if (/^PhysMem:.*,\s+(\d+)M\s+unused/) {
            $total_mem += $1;
        }    
    }

    close MEM;
    unless ($total_mem) {
        confess("Could't parse physical memory from the output of top");
    }

    return $total_mem;

}

sub freeMemoryForShadow
{
undef %resPDL;undef %connPDL;
undef %connAngleFunctionals;undef %connDiheFunctionals;
undef %connBondFunctionals;undef $contactPDL;
}


##
## ALL-ATOM PROGRAM DRIVER ##
## Perform All-Atom Topology Generation
sub allAtom
{
	my ($inputFolder,$inputPDB,$groFile,$ndxFile,$shadowFile,$topFile) = @_;

	print "\n************************\n";
	print "******* ALL ATOM *******";
	print "\n************************\n";
        
        ## TEST IF SMOG PATH EXISTS ##
        if(!$ENV{"SMOG_PATH"}){confess("SMOG_PATH environment variable not set\n");}
	$| = 1;
	## OPEN .top handle ##
	open($OUTPUT, ">$topFile");
        print "Parsing All-Atom Templates.....";
	parseInputFolder($inputFolder); ## PARSE FOLDER --> FILE

	#####################
        ## PARSE TEMPLATES ##
        #####################

	## parse Bif File ##
	parseBif();
	## parse Sif File ##
	parseSif();
	## parse Bonds/Nonbonds File
	parseBonds();parseNonBonds();
        print "Done\n";

	## Find Bonds ##
	createBondFunctionals();
	## Find Angles ##
	createDihedralAngleFunctionals();

	###############
        ## PARSE PDB ##
        ###############
	print "Parsing .pdb.....";
	parseATOM($inputPDB);
	catPDL();
	print "Done\n";

	#####################################
	## DO CALCULATIONS FOR HAMILTONIAN ##
	#####################################
	## COUNT NUMBER OF DIHEDRALS ##
	getSetDiheCounts(\%connDiheFunctionals,\%connPDL);
	## SET DIHEDRAL RATIOS ##
	setRatios(\%connDiheFunctionals,\%connPDL,$totalAtoms,\%allAtoms);
	## PRINT TOPOLOGY HEADERS ##
	printDefaultParamsHeader();
	## CALCULATE ATOM TYPE HEADERS ##
	calculateAtomTypes(\%allAtoms);
	## PRINT MOLECULE TYPE ##
	printMoleculeTypes();
	## PRINT ATOMS HEADER ##
	calculateAtoms($inputPDB);
	## CALCULATE BONDS/CONNECTED BONDS ##
	printOrderedBonds(\%bondFunctionals,\%resPDL,\%connBondFunctionals,\%connPDL,0);
	## CALCULATE ANGLES ##
	printOrderedAngles(\%connAngleFunctionals,\%connPDL);
	## CALCULATE DIHEDRALS ##
	printOrderedDihedrals(\%connDiheFunctionals,\%connPDL);

    

	## CREATE GRO FILE ##
	print "Creating .gro...";
	convertPDBToGroNdx($inputPDB,$groFile,$ndxFile);
	print "Done\n";
	## CALL SHADOW/PARSE CONTACTS ##
	print "Finding contacts via shadow algorithm...";
	
	
	my $SCMparams = setContactParams();
	my $absAtomCount = keys(%allAtoms);
    my $memoryMax = "";
    if($absAtomCount >= 1000000)
    {
     $memoryMax = 4000;
     print "JAVA HEAP SIZE INCREASED to $memoryMax Mb\n";
     $memoryMax = "-Xmx$memoryMax"."m";
    } 
	
	## Delete Preexisting contact file ##
	if(-e $shadowFile){unlink($shadowFile);}
	
	## Free Memory for Shadow ##
	freeMemoryForShadow();
	
	if((!exists $ENV{SMOG_PATH}) || 
		!(-e "$ENV{SMOG_PATH}/tools/SCM.jar")){confess "Can't find Shadow. Make sure SMOG_PATH is set correctly.\n"}
	system("java $memoryMax -jar $ENV{SMOG_PATH}/tools/SCM.jar -g $groFile -t $topFile -o $shadowFile  $SCMparams > shadow.log");
	$numContacts = parseCONTACT($shadowFile,"",$ignChainContact,0);
	print "Done\n";


	if($numContacts == 0){confess "There are 0 contacts,check shadow.error log file for possible errors";}
	print "NOTE: Total atom count used for normalization is $totalAtoms with $absAtomCount atoms in the system\n";
	print "Calculating Hamiltonian...";
	## CALCULATE PAIRS and EXCLUSION ##
 	calculateContacts($contactPDL,\%resPDL,\%allAtoms,$numContacts,$totalAtoms);
	calculateExclusions($shadowFile,$ignChainContact);
	## PRINT TOPOLOGY FOOTERS ##
	printDefaultParamsFooter();
	print "Done\n";
	print "Your Hamiltonian is ready!\nFiles generated $topFile,$groFile,$shadowFile,$ndxFile,shadow.log\n";
	close($OUTPUT);
        print "Press any key to continue...\n";
	#<STDIN>; ## Wait for user return #
	printCitation();
}

##
## COARSE GRAIN PROGRAM DRIVER ##
## Perform Coarse Grain Topology Generation
sub coarseGrain
{
	my ($inputFolder,$coarseFolder,$inputPDB,$groFile,$ndxFile,$shadowFile,$topFile) = @_;

	print "\n**********************************\n";
	print "********* COARSE GRAINED *********";
	print "\n**********************************\n";
	$| = 1;
	## TEST IF SMOG PATH EXISTS ##
        if(!$ENV{"SMOG_PATH"}){confess("SMOG_PATH environment variable not set\n");}

	## OPEN .top handle ##
	open($OUTPUT, ">$topFile");
	

	##############################
        ## PARSE ALL ATOM TEMPLATES ##
        ##############################
	print "Parsing All-Atom Templates.....";
	parseInputFolder($inputFolder);
	## parse Bif File ##
	parseBif();
	## parse Sif File ##
	parseSif();
	## parse Bonds/Nonbonds File
	parseBonds();parseNonBonds();
        print "Done\n";

	## Find Bonds ##
	createBondFunctionals();
	## Find Angles ##
	createDihedralAngleFunctionals();
	###############
        ## PARSE PDB ##
        ###############
	print "Parsing .pdb.....";
	parseATOM($inputPDB);
	catPDL();
	print "Done\n";

	#####################################
	## DO CALCULATIONS FOR HAMILTONIAN ##
	#####################################
	## COUNT NUMBER OF DIHEDRALS ##
	getSetDiheCounts(\%connDiheFunctionals,\%connPDL);
	## SET DIHEDRAL RATIOS ##
	setRatios(\%connDiheFunctionals,\%connPDL,$totalAtoms,\%allAtoms);
	## PRINT TOPOLOGY HEADERS ##
	printDefaultParamsHeader();
	## CALCULATE ATOM TYPE HEADERS ##
	calculateAtomTypes(\%allAtoms);
	## PRINT MOLECULE TYPE ##
	printMoleculeTypes();
	## PRINT ATOMS HEADER ##
	calculateAtoms($inputPDB);
	## CALCULATE BONDS/CONNECTED BONDS ##
	printOrderedBonds(\%bondFunctionals,\%resPDL,\%connBondFunctionals,\%connPDL,0);
	## CALCULATE ANGLES ##
	printOrderedAngles(\%connAngleFunctionals,\%connPDL);
	## CALCULATE DIHEDRALS ##
	printOrderedDihedrals(\%connDiheFunctionals,\%connPDL);

        print "Finding coarse grain contacts.....";
        convertPDBToGroNdx($inputPDB,$groFile,$ndxFile);
	my $SCMparams = setContactParams();
	my $absAtomCount = keys(%allAtoms);
    my $memoryMax = "";
    if($absAtomCount >= 1000000)
    {
     $memoryMax = 4000;
     print "JAVA HEAP SIZE INCREASED to $memoryMax Mb\n";
     $memoryMax = "-Xmx$memoryMax"."m";
    } 
    
	if((!exists $ENV{SMOG_PATH}) || 
		!(-e "$ENV{SMOG_PATH}/tools/SCM.jar")){confess "Can't find Shadow. Make sure SMOG_PATH is set correctly.\n"}
	system("java $memoryMax -jar $ENV{SMOG_PATH}/tools/SCM.1.31.jar --coarse CA -g $groFile -t $topFile -o $shadowFile $SCMparams  > shadow.log");
        print "Done\n";
	
	## CLEAR ALL-ATOM MEMORY ##
        clearBifMemory();clearPDBMemory();
	###################################
        ## PARSE COARSEGRAINED TEMPLATES ##
        ###################################

	## OPEN .top handle ##
	open($OUTPUT, ">$topFile");
	print "Parsing Coarse Grained Templates.....";
	parseInputFolder($coarseFolder);
	## parse Bif File ##
	parseBif();
	## parse Sif File ##
	parseSif();
	## parse Bonds/Nonbonds File
	parseBonds();parseNonBonds();
        print "Done\n";

	## Find Bonds ##
	createBondFunctionals();
	## Find Angles ##
	createDihedralAngleFunctionals();
        
	###############
    ## PARSE PDB ##
    ###############
	print "Coarse Graining system.....";
	parseATOMCoarse($inputPDB);
	catPDL();
	print "Done\n";
	
	## CREATE GRO FILE ##
	print "Creating .gro.....";
	convertPDBToGroNdx($inputPDB,$groFile,$ndxFile);
	print "Done\n";
	$absAtomCount = keys(%allAtoms);
	print "NOTE: Total atom count used for normalization is $totalAtoms with $absAtomCount atoms in the system\n";
	print "Creating Topology.....";
	
	## COUNT NUMBER OF DIHEDRALS ##
	getSetDiheCounts(\%connDiheFunctionals,\%connPDL);
	## SET DIHEDRAL RATIOS ##
	setRatios(\%connDiheFunctionals,\%connPDL,$totalAtoms,\%allAtoms);
	## PRINT TOPOLOGY HEADERS ##
	printDefaultParamsHeader();
	## CALCULATE ATOM TYPE HEADERS ##
	calculateAtomTypes(\%allAtoms);
	## PRINT MOLECULE TYPE ##
	printMoleculeTypes();
	## PRINT ATOMS HEADER ##
	calculateAtoms($inputPDB);
	## CALCULATE BONDS/CONNECTED BONDS ##
	printOrderedBonds(\%bondFunctionals,\%resPDL,\%connBondFunctionals,\%connPDL,1);
	## CALCULATE ANGLES ##
	printOrderedAngles(\%connAngleFunctionals,\%connPDL);
	## CALCULATE DIHEDRALS ##
	printOrderedDihedrals(\%connDiheFunctionals,\%connPDL);

	$numContacts = parseCONTACT($shadowFile,$DCAFile,$ignChainContact,1);
	if($numContacts == 0){confess "There are 0 contacts,check shadow.error log file for possible errors";}
	## CALCULATE PAIRS and EXCLUSION ##
 	calculateContacts($contactPDL,\%resPDL,\%allAtoms,$numContacts,$totalAtoms);
	calculateExclusions($shadowFile,$ignChainContact);
	## PRINT TOPOLOGY FOOTERS ##
	printDefaultParamsFooter();
	print "Done\n";
	print "Your Hamiltonian is ready!\nFiles generated $topFile,$groFile,$shadowFile,$ndxFile,shadow.log\n";
	close($OUTPUT);
	print "Press any key to continue...\n";
	#<STDIN>; ## Wait for user return ##
	printCitation();




}

##################
## PARSE INPUTS ##
##################

## parse command line options ##
# -o output top file name, default is temp.top
# -v verbose option
usage() if ( @ARGV < 1 or
          !GetOptions('help|?' => \$help, 'o=s' => \$topFile, 'tAA=s' => \$inputFolder, 'i=s' => \$inputPDB, 'g=s' => \$groFile,
	'c=s' => \$DCAFile, 'ignH' => \$ignH, 'ignChainContacts' => \$ignChainContact, 's=s' => \$shadowFile, 'n=s' => \$ndxFile,
	 'tCG=s' => \$coarseFolder, 'contactRes=s' => \$contactRes, 'ignHContacts'=>\$ignHContacts,
	 'softTerm' => \$softTerm)
          or defined $help 
          or !defined $inputPDB 
          or (!defined $inputFolder  && !defined $contactRes));

## Coarse grain input check ##
if($contactRes eq "" && $coarseFolder ne "")
{confess("Please specify a contact resolution when using -tCG option with -contactRes option\n");}
if($contactRes ne "" && $coarseFolder ne "")
{$inputFolder=$contactRes;} 

## PreProcessing ##
## Check if -ignh option is set ##
stripHydrogen($inputPDB) if defined $ignH;

####################
## PROGRAM DRIVER ##
####################

## PROCESS All-Atom or Coarse Grain ##
## ALL-ATOM ##
if($coarseFolder eq ""){allAtom($inputFolder,$inputPDB,$groFile,$ndxFile,$shadowFile,$topFile);}
## COARSE-GRAIN ## 
else{coarseGrain($inputFolder,$coarseFolder,$inputPDB,$groFile,$ndxFile,$shadowFile,$topFile);}
## END ##
