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
# templateParser: Parses all XML files into a hash 
#####################################################################
package templateParser;

#####################
## COMPILE HEADERS ##
#####################
use strict;
use warnings;
####################
## MODULE HEADERS ##
####################
use XML::Simple qw(:strict);
use Exporter;
use String::Util 'trim';
use Storable qw(dclone);
use smog_common;

## DECLARATION TO SHARE DATA STRUCTURES ##
our @ISA = 'Exporter';
our @EXPORT = 
qw(checkFunction getEnergyGroup $energyGroups $interactionThreshold $countDihedrals $termRatios %residueBackup %fTypes $functions %eGRevTable %eGTable intToFunc funcToInt %residues %bondFunctionals %angleFunctionals %connections %dihedralAdjList adjListTraversal adjListTraversalHelper $interactions setInputFileName parseBif parseSif parseBonds createBondFunctionals createDihedralAngleFunctionals parseNonBonds getContactFunctionals $contactSettings clearBifMemory @topFileBuffer @linesInDirectives Btypespresent NBtypespresent PAIRtypespresent EGinBif checkenergygroups bondtypesused pairtypesused checkBONDnames checkNONBONDnames checkPAIRnames checkREScharges checkRESimpropers round);

######################
## GLOBAL VARIABLES ##
######################

my %SMOGversions;
my $smi=0;
foreach my $ver("2.0", "2.0.1", "2.0.2", "2.0.3", "2.1", "2.2", "2.3"){
 $SMOGversions{$ver}=$smi;
 $smi++; 
}
#########################
## XML PARSED VARIBLES ##
#########################

## HOLDS INFORMATION ON RESIDUES ##
## residue => 
##			"residueType" => residue type (ie. amino,rna),
##			"atoms" => hash of atoms with nbtype, btype,
##			"impropers" => list of 4 atom impropers
##			"bonds" => hash of bonds info with key as "atomA-atomB"
##
our %residues;
our %residueBackup;
our $functions;
our $contactSettings;
our $interactionThreshold;
our $countDihedrals;

my $settings; our $termRatios;
our $interactions;
our %funcTable;our %funcTableRev;
our %eGTable;our %eGRevTable;
our @topFileBuffer;our @linesInDirectives;
## SORTED FUNCTIONALS ##
our %bondFunctionals;
our %dihedralFunctionals;
our %angleFunctionals;
our %connections;
our %dihedralAdjList;

## Create new XML::Simple object ##
my $xml = new XML::Simple;
my $data;
my $residueHandle;

## INPUT FILE NAMES ##
our $bifxml = "bif.xml";
our $sifxml = "sif.xml";
our $bondxml = "b.xml";
our $nbondxml = "nb.xml";

our %Btypespresent;
our %NBtypespresent;
our %PAIRtypespresent;
our %EGinBif;	
our %EGinSif;	
our %pairtypesused;
our %bondtypesused;
our %fTypes;


my %bondHandle;
my @improperHandle;
my %dihedralHandle;

###########################
## CLEAR VARIABLE MEMORY ##
###########################
sub clearBifMemory {
	## Store bif info for PDB looping
	%residueBackup = %{ dclone (\%residues) };
	undef %residues;
	undef $functions;
	undef $contactSettings;
	undef $termRatios;
	undef $interactions;
	undef %funcTable;
	undef %funcTableRev;
	undef %eGRevTable;
	undef %bondFunctionals;
	undef %dihedralFunctionals;
	undef %angleFunctionals;
	undef %connections;
	undef %dihedralAdjList;
	undef @topFileBuffer;
	undef @linesInDirectives;
        undef %NBtypespresent;
        undef %Btypespresent;
        undef %PAIRtypespresent;
	undef %pairtypesused;
	undef %bondtypesused;
}

########################
## SET INPUTFILE NAME ##
########################
sub setInputFileName {
	my ($a,$b,$c,$d,$e) = @_;
	$bifxml = $a;
	$sifxml = $b;
	$bondxml = $c;
	$nbondxml = $d;
}

sub round
{
	my ($val)=@_;
	my $round;
	if($val<0){
   		$round=int(abs($val)+0.5);
		$round=$round*(-1);
	}else{
		$round=int($val+0.5);
 	}
	return $round;	
}

# checkFunction: verifies that function type is supported 
sub checkFunction
{
	my($funcString) = @_;
	$funcString =~ s/\(.*//g;
	if(!exists $fTypes{"$funcString"}){smog_quit ("\"$funcString\" is not a supported function type in SMOG");}
	if(!exists $functions->{$funcString}){smog_quit ("Function \"$funcString\" is being used, but is not defined in .sif file");}
}

sub checkREScharges
{
	my $string="";
	foreach my $res(keys %residues){
		my $charge=0;
		foreach my $at(keys %{$residues{$res}->{"atoms"}}){
			my $atomType = $residues{$res}->{"atoms"}->{$at}->{"nbType"};
			if(defined $residues{$res}->{"atoms"}->{$at}->{"charge"}){
				# explicit charge on specific atoms takes priority
				$charge+=$residues{$res}->{"atoms"}->{$at}->{"charge"};
			}elsif(exists $interactions->{"nonbonds"}->{$atomType}){
				# charge added based on nbType
				$charge+=$interactions->{"nonbonds"}->{$atomType}->{"charge"};
			}
		}
		if(defined $residues{$res}->{"totalcharge"}){
			my $tc=$residues{$res}->{"totalcharge"};
			# check that total charge is the value expected
			if(abs($charge-$tc)>0.001){
				$charge=int($charge*10**8)/10**8;
				$string .= "Residue $res has total charge of $charge, but \"totalcharge\" defined in .bif is $tc\n";
			}
		}else{
			# check that charge is integer
			if(abs($charge - round($charge)) >0.001){
				my $t=round($charge);
				$string .="Residue defined in .bif, $res, has non-integer total charge of $charge\n";
			}
		}
	}
	return $string;
}

sub checkRESimpropers
{
	my $string="";
	foreach my $res(keys %residues){
		# handle for bonds in residue $res
		my $bondshandle=$residues{$res}->{"bonds"};
		my @improperHandle = @{$residues{$res}->{"impropers"}};
               	foreach my $improper(@improperHandle){
			my @IMPARR;
			for(my $I=0;$I<4;$I++){
				$IMPARR[$I]=$improper->[$I]
			}
			my $passed=0;
			my @ips=@IMPARR;
			for(my $I=0;$I<4;$I++){
				my %blist;
				my $centeratom=$IMPARR[0];
				# shift so that we always check the first three
				$IMPARR[0]=$IMPARR[1];
				$IMPARR[1]=$IMPARR[2];
				$IMPARR[2]=$IMPARR[3];
				$IMPARR[3]=$centeratom;
				# make list of atoms that are bonded to $I
				foreach my $tb(keys %{$bondshandle}){
					my @B=split(/-/,$tb);
					if("$B[0]" eq "$centeratom"){
						$blist{$B[1]}=1;
					}
					if("$B[1]" eq "$centeratom"){
						$blist{$B[0]}=1;
					}
				}
			
				# check bond list and see if the other three atoms are bonded to the possible center atom
				my @found=(0,0,0);
				foreach my $ba(keys %blist){
					for(my $I=0;$I<3;$I++){
						if("$ba" eq "$IMPARR[$I]"){
							$found[$I]=1;
						}
					}
				}
				if($found[0]==1 && $found[1]==1 && $found[2]==1){
					$passed=1;
				}
			}
			if($passed == 0){
				$string  .= "Improper dihedral in residue $res, defined by atoms @ips, is not defined by three atoms bonded to a central atom.  There may be missing bonds, or an incorrectly-defined improper, in the .bif file.\n";
			}
		}
	}
	return $string;
}

sub checkBONDnames
{
	my $string="";
	foreach my $name(keys %bondtypesused){
		unless($name =~ /^[a-zA-Z0-9_]+$|^\*$/){
			$string .="Only letters, numbers and _, or a solitary *, can appear in bond/angle/dihedral definitions. bType \"$name\" encountered\n";
		}
		if($name ne "*"  && !defined $Btypespresent{$name}){	
			$string .="bType $name appears in .b file, but doesn't appear anywhere in .bif.  Likely a typo in your .b file.\n";
		}
	}

	foreach my $name(keys %Btypespresent){
		unless($name =~ /^[a-zA-Z0-9_]+$/){
			$string .="Only letters, numbers or _ can be used in a bType definition of a residue. bType \"$name\" encountered\n";
		}
		if(!defined $bondtypesused{$name} && !defined $bondtypesused{"*"}){	
			$string .="bType $name appears in .bif file, but doesn't appear to have parameters defined in the .b file.\n";
		}
	}

	return $string;
}

sub checkNONBONDnames
{
	my $string="";
	foreach my $name(keys %{$interactions->{"nonbonds"}}){
		unless($name =~ /^[a-zA-Z0-9_]+$|^\*$/){
			$string .="Only letters, numbers and _, or a solitary *, can appear in nobonded definitions. nbType \"$name\" encountered\n";
		}
		if($name ne "*"  && !defined $NBtypespresent{$name}){	
			$string .="nbType $name appears in .nb file, but doesn't appear anywhere in .bif.  Likely a typo in your .nb file.\n";
		}
	}

	foreach my $name(keys %NBtypespresent){
		unless($name =~ /^[a-zA-Z0-9_]+$/){
			$string .="Only letters, numbers or _ can appear in nbType name within a residue. nbType \"$name\" encountered\n";
		}
		if(!defined ${$interactions->{"nonbonds"}}{$name} && !defined ${$interactions->{"nonbonds"}}{"*"}){	
			$string .="nbType $name appears in .bif file, but parameters are not defined in the .nb file. \n";
		}
	}


	return $string;
}

sub checkPAIRnames
{
	my $string="";
	foreach my $name(keys %pairtypesused){
		unless($name =~ /^[a-zA-Z0-9_]+$|^\*$/){
			$string.="Only letters, numbers and _, or a solitary *, can appear in contact definitions. pairType \"$name\" encountered\n";
		}
		if($name ne "*"  && !defined $PAIRtypespresent{$name}){	
			$string .="pairType $name appears in .nb file, but doesn't appear anywhere in .bif.  Likely a typo in your .nb file.\n";
		}
	}

	foreach my $name(keys %PAIRtypespresent){
		unless($name =~ /^[a-zA-Z0-9_]+$/){
			$string.="Only letters, numbers or _ can appear in pairType names within a residue definition. pairType \"$name\" encountered\n";
		}
		if( !defined $pairtypesused{$name} && !defined $pairtypesused{"*"} ){	
			$string .="pairType $name appears in .bif file, but parameters are not defined in the .nb file.\n";
		}
	}

	return $string;
}

sub checkenergygroups
{
	my $messagestring="";
	my $string="It appears that dihedrals energy groups are only partially defined. Since dihedrals energy groups must be defined in the .bif (residue->bonds), .sif (energyGroup) and .b (bonds) files, in order for the interaction to be applied, a partial declaration is probably unintentional.  Specific issues listed below:\n\n";
	foreach my $II (sort keys %{$interactions->{"dihedrals"}}){
		if(! exists $EGinBif{$II}){
			$messagestring .="energyGroup \"$II\" defined in .b file, but is not used in .bif file.\n";
		}
		if(! exists $EGinSif{$II}){
			$messagestring .="energyGroup \"$II\" defined in .b file, but is not used in .sif file.\n";
		}
	}

	foreach my $II (sort keys %EGinBif){
		if(! exists $interactions->{"dihedrals"}->{$II}){
			$messagestring .="energyGroup \"$II\" defined in .bif file, but is not used in .b file.\n";
		}
		if(! exists $EGinSif{$II}){
			$messagestring .="energyGroup \"$II\" defined in .bif file, but is not used in .sif file.\n";
		}
	}

	foreach my $II (sort keys %EGinSif){
		if(! exists $interactions->{"dihedrals"}->{$II}){
			$messagestring .="energyGroup \"$II\" defined in .sif file, but is not used in .b file. \n";
		}
		if(! exists $EGinBif{$II}){
			$messagestring .="energyGroup \"$II\" defined in .sif file, but is not used in .bif file.\n";
		}
	}

# clear the hashes, in case we need them later.
	undef %EGinBif;
	undef %EGinSif;
	if($messagestring ne ""){
		$string.=$messagestring;
	}else{
		$string="";
	}
	return $string;
}

sub parseBif {
	my ($tempdir,$maponly)=@_;
	## Read .bif ##
	my $data = $xml->XMLin($bifxml,KeyAttr=>{residue=>"name",connection=>"name"},ForceArray=>1);
	
	## PARSE RESIDUES INTO A HASH ##
	## Hash is formatted as below
	## residue => 
	##			"residueType" => residue type (ie. amino,rna),
	##			"atoms" => hash of atoms with nbtype, btype,pairType
	##			"impropers" => list of 4 atom impropers
	##			"bonds" => hash of bonds info with key as "atomA-atomB"
	##
	
	## Obtain handle to loop through residue
	my $residueHandle = $data->{"residues"}->[0]->{"residue"};

	# generate a map for adjustPDB
	if(defined $maponly){
		print "\nWill only make the map $maponly file and then exit\n\n";
		open(GENMAP,">$maponly") or smog_quit("unable to open $maponly for writing");
		print GENMAP ";This is a mapping file that was generated by smog2\n";
		print GENMAP ";This defines the composition of atoms in the templates found in:\n";
		print GENMAP ";   $tempdir\n";
		my @mstring="";
		my $ti=0;
		foreach my $res ( keys %{$residueHandle} )
		{
			$mstring[$ti] = "residue $res";
			my @atomHandle = @{$residueHandle->{$res}->{"atoms"}->[0]->{"atom"}};
			my %list;
			foreach my $atom(@atomHandle)
			{
				my $AT= $atom->{"content"};
				$list{$AT}=1;
			}
			foreach my $ll(sort keys %list){
				$mstring[$ti] .= " $ll";
			}
			if(defined $residueHandle->{$res}->{"meta"}){
				$mstring[$ti] .= " %";
				$mstring[$ti] .= $residueHandle->{$res}->{"meta"};
			}
			$mstring[$ti] .= "\n";
			$ti++;
		}
		foreach my $tt(sort @mstring){
			print GENMAP $tt;
		}
		close(GENMAP);
		return;
	}

	
	## Loop through residues
	foreach my $res (keys %{$residueHandle} )
	{
		## CREATE ATOM HASH ##
		my %atoms; my $index = 0;
		# Obtain handle to loop through atoms
		my @atomHandle = @{$residueHandle->{$res}->{"atoms"}->[0]->{"atom"}};
		my %seen;
		foreach my $atom(@atomHandle)
		{
			my $AT= $atom->{"content"};
			if(exists $seen{"$AT"}){smog_quit("Error in .bif. Duplicate declaration of atom $AT in residue $res.")};
			$seen{"$AT"}=1;
		
			unless($atom->{"nbType"} =~ /^[a-zA-Z0-9_]+$/){
		 		my $T=$atom->{"nbType"};
		 		smog_quit("Only letters, numbers and _ can appear in nbType definitions. nbType \"$T\" found in residue $res");
			}
			unless($atom->{"bType"} =~ /^[a-zA-Z0-9_]+$/){
		 		my $T=$atom->{"bType"};
		 		smog_quit("Only letters, numbers and _ can appear in bType definitions. nbType \"$T\" found in residue $res");
			}
			unless($atom->{"pairType"} =~ /^[a-zA-Z0-9_]+$/){
		 		my $T=$atom->{"pairType"};
		 		smog_quit("Only letters, numbers and _ can appear in pairType definitions. nbType \"$T\" found in residue $res");
			}
			$NBtypespresent{$atom->{"nbType"}}=1;
			$Btypespresent{$atom->{"bType"}}=1;
			$PAIRtypespresent{$atom->{"pairType"}}=1;

		      	$atoms{$atom->{"content"}} = {"index"=>$index,"nbType" => $atom->{"nbType"},"bType" => $atom->{"bType"},
		      	"pairType" => $atom->{"pairType"}, "charge" => $atom->{"charge"}};
		      
		      	$index++;
		}
		undef %seen;
		
		## CREATE IMPROPER ARRAY ##
		my @impropers;my @improperHandle;
		# Obtain handle to loop through impropers
		if(exists $residueHandle->{$res}->{"impropers"})
			{@improperHandle = @{$residueHandle->{$res}->{"impropers"}->[0]->{"improper"}};}
		my %seenIMP;
		foreach my $improper(@improperHandle)
		{
		  ## [[A,B,C,D],[E,F,G,H],...]
		      	if(!exists $improper->{"atom"}->[0]){smog_quit("Declaration of residue $res has an improper that lacks atoms\n")};
			my %seenAtom;
		      	my $atomstring=$improper->{"atom"}->[0];
		      	my $atomstringRev=$improper->{"atom"}->[0];
		      	for(my $I=1;$I<4;$I++){
		      		my $T=$improper->{"atom"}->[$I];
		      		$atomstring=$atomstring . "-" . $T;
		      		$atomstringRev=$T . "-" . $atomstring;
	        		if($T =~ /[?^&!@#%()-]/){smog_quit ("Special characters not permitted in \"connection\" atom names: $T found.")};
		      		if(exists $seenAtom{"$T"}){
		      			smog_quit("Error in .bif.  Duplicate declaration of atom $T in improper dihedral for residue $res.");
		      		}else{
		      			$seenAtom{"$T"}=1;
		      		}
		      	}
		
		      	if(exists $seenIMP{"$atomstring"} or exists $seenIMP{"$atomstringRev"}){
		      		smog_quit("Error in .bif.  Duplicate declaration of improper dihedral $atomstring for residue $res.");
		      	}else{
		      		$seenIMP{"$atomstring"}=1;
		      		$seenIMP{"$atomstringRev"}=1;
		      	}
				undef %seenAtom;
			push(@impropers,$improper->{"atom"});
		}
		undef %seenIMP;
		## CREATE BOND HASH ##
		my %bonds; my %energyGroups; 
		my @bondHandle;
		# Obtain handle to loop through bonds
		if(exists $residueHandle->{$res}->{"bonds"})
		      {@bondHandle = @{$residueHandle->{$res}->{"bonds"}->[0]->{"bond"}};}
		foreach my $bond(@bondHandle)
		{
		  ## bonds{atomA-atomB} = bond info from XML
			my $atomA = $bond->{"atom"}->[0];
			my $atomB = $bond->{"atom"}->[1];
			if(exists $bonds{"$atomA-$atomB"}){
				smog_quit("Error in .bif.  Duplicate declaration of bond between atoms $atomA and $atomB in residue $res.");
			}
			$bonds{"$atomA-$atomB"} = $bond;
			
			$energyGroups{"$atomA-$atomB"} = $bond->{"energyGroup"};
			$energyGroups{"$atomB-$atomA"} = $bond->{"energyGroup"};
			# log what energyGroups have been used in .bif
			$EGinBif{$bond->{"energyGroup"}}=1;	
		}
		
		##atomCount !exists == -1, else atomCount
		if(!exists $residueHandle->{$res}->{"atomCount"})
		{
			$residueHandle->{$res}->{"atomCount"}=-1;
		}

		## Create residue hash containing all data
		my $interRes = {
		      "residueType" => $residueHandle->{$res}->{residueType},
		      "atoms" => \%atoms,
		      "impropers" => \@impropers,
		      "bonds" => \%bonds,
		      "energyGroups" => \%energyGroups,
		      "atomCount" => $residueHandle->{$res}->{"atomCount"},
		      "totalcharge" => $residueHandle->{$res}->{"totalcharge"},
		      "connect" => $residueHandle->{$res}->{"connect"},
		      "meta" => $residueHandle->{$res}->{"meta"}
		      };
		$residues{$res} = $interRes;
	  
	}

	## PARSE CONNECTIONS ##
	## Parse connections into an array of connection information.
	## The index will represent all the unique type of connections
	##
	
	## Obtain handle to loop through connections
	my $conHandle = $data->{"connections"}->[0]->{"connection"};
	my $resA; my $resB; 
	## Loop through connections
	foreach my $connname (keys %{$conHandle})
	{
		$resA = $conHandle->{$connname}->{"residueType1"};
		$resB = $conHandle->{$connname}->{"residueType2"};
		if(exists $connections{$resA}->{$resB}){
			smog_quit ("Multiple connections defined between residueTypes $resA and $resB");
		}
		$connections{$resA}->{$resB}=$conHandle->{$connname};
		foreach my $II(@{$conHandle->{$connname}->{"bond"}}){
			$EGinBif{$II->{"energyGroup"}}=1;
		} 
	}
}

sub parseSif {

	## Read .sif file ##
	$data = $xml->XMLin($sifxml,KeyAttr => ['name'],ForceArray=>1);
	if(!exists $data->{"version"}->[0]->{"min"}){
		smog_quit("Minimum required SMOG version is not defined in .sif file. This check is intended to ensure that one does not use new templates with an old version of SMOG.  However, newer version of SMOG should always work with old templates.  Accordingly, if you are using the newest release of SMOG, you can probably ignore this warning.",0);
	}else{
		my $minver=$data->{"version"}->[0]->{"min"};
		if(!exists $SMOGversions{"$minver"}){
			smog_quit("These templates require SMOG v$minver, or newer (defined in .sif file).");
		}
	}
	## Parse function data
	$functions = $data->{"functions"}->[0]->{"function"};
	foreach my $funcName(keys %{$functions}){
	
		if($functions->{$funcName}->{"directive"} eq "pairs" 
			&& !exists $functions->{$funcName}->{"exclusions"}){
			smog_quit( "Since $funcName is of directive \"pairs\", boolean element \"exclusions\" must be included in the function declaration in the .sif file.\n");
		}elsif(exists $functions->{$funcName}->{"exclusions"} && $functions->{$funcName}->{"exclusions"} ==1 
			&& $functions->{$funcName}->{"exclusions"} ==0){
			smog_quit("function $funcName element exclusions must be 0, or 1.");
		}elsif($functions->{$funcName}->{"directive"} ne "pairs"
	                && exists $functions->{$funcName}->{"exclusions"}){
			print "\nNOTE: Element \"exclusions\" is defined for function $funcName. This is likely unnecessary, since the \"exclusions\" element is only relevant for contacts.\n";
		}
	}
	## Parse settings data
	$settings = $data->{"settings"}->[0];
	our $energyGroups = $settings->{"Groups"}->[0]->{"energyGroup"};
	my $contactGroups = $settings->{"Groups"}->[0]->{"contactGroup"};
	my $groupRatios = $settings->{"Groups"}->[0]->{"groupRatios"}->[0];
	my $residueType; 
	my $intraRelativeStrength; my $normalize;
	my $total;
	
	## PARSE ENERGY GROUP INFORMATION ##
	## INFO PLACED IN termRatio HASH
	## $termRatio->"residueType"->"energyGroup" = 
	##	"relativeStrength" -> factor to determine ratios
	##  "normalize" -> if strength is to be normalized among bonds
	##
	
	##
	# Contact to Dihedral Group ratio is global
	# Contact ratio is global
	# Dihedral group ratio is residue dependent
	my $EG_NORM=0;
	foreach my $egName(keys %{$energyGroups})
	{
		$EGinSif{$egName}=1;
		$residueType = $energyGroups->{$egName}->{"residueType"};
		$intraRelativeStrength = $energyGroups->{$egName}->{"intraRelativeStrength"};
		$normalize = $energyGroups->{$egName}->{"normalize"};
		$termRatios->{$residueType}->{"energyGroup"}->{$egName}={"normalize"=>$normalize,"intraRelativeStrength"=>$intraRelativeStrength};
		if($normalize == 0 && exists $energyGroups->{$egName}->{"intraRelativeStrength"}){
			smog_quit("Issue in .sif, energy group $egName. intraRelativeStrength only supported when normalization is on.");
		}elsif($normalize == 1 && !exists $energyGroups->{$egName}->{"intraRelativeStrength"}){
			smog_quit("Issue in .sif, energy group $egName. intraRelativeStrength must be set if normalization is on.");
		}elsif($normalize != 1 && $normalize != 0){
			smog_quit("Issue in .sif, energy group $egName. normalization must be 0, or 1. Found $normalize");
		}elsif(exists $energyGroups->{$egName}->{"intraRelativeStrength"} && $energyGroups->{$egName}->{"intraRelativeStrength"} <=0 ){
	                smog_quit("intraRelativeStrength must be >= 0.  See energyGroup $egName in .sif.");
	        }
		
		if($normalize==1){
			$EG_NORM++;
		}
		
	}
	
	my $CG_NORM=0;
	$total=0;
	my $numCGs=0;
	foreach my $egName(keys %{$contactGroups})
	{
		$numCGs++;
		$intraRelativeStrength = $contactGroups->{$egName}->{"intraRelativeStrength"};
		$normalize = $contactGroups->{$egName}->{"normalize"};
		$termRatios->{"contactGroup"}->{$egName}={"normalize"=>$normalize,"intraRelativeStrength"=>$intraRelativeStrength};
		##if($normalize){$total+=$intraRelativeStrength;$CG_NORM++}
		if($normalize){$CG_NORM++}
	
		if($normalize == 0 && exists $contactGroups->{$egName}->{"intraRelativeStrength"}){
			smog_quit("Issue in .sif, contact group $egName. intraRelativeStrength only supported when normalization is on.");
		}elsif($normalize == 1 && !exists $contactGroups->{$egName}->{"intraRelativeStrength"}){
			smog_quit("Issue in .sif, contact group $egName. intraRelativeStrength must be set if normalization is on.");
		}elsif($normalize != 1 && $normalize != 0){
			smog_quit("Issue in .sif, contact group $egName. normalization must be 0, or 1. Found $normalize");
		}elsif(exists $contactGroups->{$egName}->{"intraRelativeStrength"} && $contactGroups->{$egName}->{"intraRelativeStrength"} <=0 ){
			smog_quit("intraRelativeStrength must be >= 0.  See contactGroup $egName .sif.");
		}
	}
	
	if(($CG_NORM > 0 and $EG_NORM ==0) or ($EG_NORM > 0 and $CG_NORM ==0)){
		smog_quit('Issue in .sif. Normalization only turned on for ContactGroups, or EnergyGroups. Normalization must be off, or on, for both.');
	}
	
	## NOTE:Contact Type is Global ##
	## Sum of contact scalings ##
	
	if($groupRatios->{"contacts"} <=0 || $groupRatios->{"contacts"} <=0){
		smog_quit("All values for groupRatios must be greater than zero. See .sif file.")
	}
	## contact/dihe relative Scale ##
	$termRatios->{"contactRelative"} = $groupRatios->{"contacts"};
	## dihe/contact relative Scale ##
	$termRatios->{"energyRelative"} = $groupRatios->{"dihedrals"};
	## Sum of total global scaling ##
	$termRatios->{"interRelativeTotal"} = $groupRatios->{"contacts"}+$groupRatios->{"dihedrals"};
	
	## PARSE CONTACT MAP SETTINGS ##
	$contactSettings = $data->{"settings"}->[0]->{"Contacts"}->[0];
	
	## check for consistency in the contact settings
	my $method = $contactSettings->{"method"};
	
	if($method =~ m/shadow/)
	{
	
		if(!exists $contactSettings->{"shadowRadiusBonded"}){
			smog_quit("When using contact method=shadow, you must supply a value for shadowRadiusBonded in the .sif file.");
		}
		if(!exists $contactSettings->{"shadowRadius"}){
			smog_quit("When using contact method=shadow, you must supply a value for shadowRadius in the .sif file.");
		}
	}
	elsif($method =~ m/cutoff/)
	{
		if(exists $contactSettings->{"shadowRadiusBonded"}){
			smog_quit("Contact method=cutoff can not use shadowRadiusBonded.  Either change the method, or remove shadowRadiusBonded in the .sif file.");
		}
		if(exists $contactSettings->{"shadowRadius"}){
			smog_quit("Contact method=cutoff can not use shadowRadius.  Either change the method, or remove shadowRadius in the .sif file.");
		}
	}else{
		smog_quit ("Contact map method $method is not supported.");
	}
	
	# Scaling Parameters #
	if(exists $contactSettings->{"contactScaling"}){
		my $contactScaling = $contactSettings->{"contactScaling"};
		my %seenScaling;
		foreach my $k(keys %{$contactScaling})
		{
		
			my $scale = $contactScaling->{$k}->{"scale"};
			my $deltaMin = $contactScaling->{$k}->{"deltaMin"};
			my $deltaMax = $contactScaling->{$k}->{"deltaMax"};
			my $atomListString = $contactScaling->{$k}->{"atomList"};
			$atomListString = trim($atomListString);
			my @atomList = split(/\s+/,$atomListString);
			if(scalar(@atomList) == 0){smog_quit("No atom list at contact scaling");}
			my %atomListHash = map {$_=>1} @atomList;
			
			my $A = $contactScaling->{$k}->{"residueType1"};
			my $B = $contactScaling->{$k}->{"residueType2"};
			if(exists $seenScaling{"$A-$B"}  || exists $seenScaling{"$B-$A"} ){
				smog_quit("In .sif file, contactScaling given twice for residueType pair $A-$B.");
			}else{
				$seenScaling{"$A-$B"}=1;
				$seenScaling{"$B-$A"}=1;
			}
			
			delete $contactScaling->{$k};
			$contactScaling->{$A}->{$B} 
			= {"deltaMin"=>$deltaMin,"deltaMax"=>$deltaMax,"scale"=>$scale,"atomList"=>\%atomListHash};
			$contactScaling->{$B}->{$A} 
			= {"deltaMin"=>$deltaMin,"deltaMax"=>$deltaMax,"scale"=>$scale,"atomList"=>\%atomListHash};
		}
	}
	
	## PARSE INTERACTION THRESHOLD SETTINGS ##
	my $bondsThreshold = $data->{"settings"}->[0]->{"bondsThreshold"}->[0];
	my $anglesThreshold = $data->{"settings"}->[0]->{"anglesThreshold"}->[0];
	my $contactsThreshold = $data->{"settings"}->[0]->{"contactsThreshold"}->[0];
	my $distanceThreshold = $data->{"settings"}->[0]->{"distanceThreshold"}->[0];
	$interactionThreshold->{"bonds"}={"shortBond"=>$bondsThreshold->{"shortBond"},
									  "longBond"=>$bondsThreshold->{"longBond"}};
	$interactionThreshold->{"angles"}={"smallAngles"=>$anglesThreshold->{"smallAngles"},"largeAngles"=>$anglesThreshold->{"largeAngles"}};
	$interactionThreshold->{"contacts"}={"shortContacts"=>$contactsThreshold->{"shortContacts"}};
	$interactionThreshold->{"distance"}={"tooShortDistance"=>$distanceThreshold->{"tooShortDistance"}};
	if(exists $data->{"settings"}->[0]->{"dihedralNormalization"}->[0]->{"dihedralCounting"}){
	 $countDihedrals=$data->{"settings"}->[0]->{"dihedralNormalization"}->[0]->{"dihedralCounting"};
        }else{
         # by default, we will count dihedrals.
         $countDihedrals=1;
        }
# depricated options.  We have left them in the schemas so that the error that would be thrown when using an old template will not be cryptic.  However, the entries are optional 
	if( $interactionThreshold->{"distance"}->{"tooShortDistance"} )
	{
		smog_quit("tooShortDistance found in .sif file. The use of tooShortDistance has been replaced with bondsThreshold. Please remove tooShortDistance from your .sif");
	} 



}

sub funcToInt
{
	my($interType,$func,$eG) = @_;
	
	if($interType eq "dihedrals")
 	{
  
		if(exists $funcTable{$interType} && exists $funcTable{$interType}->{$eG} && exists $funcTable{$interType}->{$eG}->{$func})
		{
			return $funcTable{$interType}->{$eG}->{$func};
		}
		else
		{
			return -1;
		}
 
	}
	elsif($interType eq "contacts")
	{
		if(exists $funcTable{"contacts"} && exists $funcTable{"contacts"}->{"func"}->{$func})
		{return $funcTable{"contacts"}->{"func"}->{$func};}
		else {return -1;}
	}
	else
	{
		if(exists $funcTable{$interType} && exists $funcTable{$interType}->{$func}) 
 		{
			return  $funcTable{$interType}->{$func};
		}
		else {
			return -1;
		}
	}
}

sub intToFunc
{
	my($interType,$int,$eG) = @_;
	if($interType eq "dihedrals")
	{
		$eG = $eGTable{$eG};
		if(exists $funcTableRev{$interType} && exists $funcTableRev{$interType}->{$eG} && exists $funcTableRev{$interType}->{$eG}->{$int})
		{
			return $funcTableRev{$interType}->{$eG}->{$int};
		}
		else{
			return -1;
		}
	}
	else
	{
		if(exists $funcTableRev{$interType} && exists $funcTableRev{$interType}->{$int})
		{
			return $funcTableRev{$interType}->{$int};
		}
		else{
			smog_quit ("NO EXISTENCE OF INTERACTION $interType,$int");
		}
	}
}

# getEnergyGroup: Return the energy group for both connected, and internal dihedrals
sub getEnergyGroup
{
	my($residuea,$residueb,$atoma,$atomb) = @_;
	my $residueIn=$residuea;
	my $residueTypea;my $residueTypeb;
	
	
 	## If Bond is internal ##
 	if(($atoma =~/(.*)\?/ && $atomb =~/(.*)\?/)
 	|| ($atoma !~/(.*)\?/ && $atomb !~/(.*)\?/))
	{
	 
		$residueIn = $residueb if($atoma =~ /\?/|| $atomb =~ /\?/);
		$atoma =~ s/\?//;$atomb =~ s/\?//;
		if(exists $residues{$residueIn}->{"energyGroups"}->{"$atoma-$atomb"})
			{return $residues{$residueIn}->{"energyGroups"}->{"$atoma-$atomb"};}
		else{smog_quit("A specified energy group for $residuea:$atoma, $residueb:$atomb doesn't exists");}
	}
 	## If Bond is between two residues ##
	elsif($atomb =~/(.*)\?/)
	{
		$residueTypea =$residues{$residuea}->{"residueType"};
		$residueTypeb =$residues{$residueb}->{"residueType"};
		if(!exists $connections{$residueTypea}->{$residueTypeb}->{"bond"}->[0]->{"energyGroup"}){
			smog_quit("Connection not defined for resTypes $residueTypea-$residueTypeb (residues $residuea $residueb)");
		}
		return $connections{$residueTypea}->{$residueTypeb}->{"bond"}->[0]->{"energyGroup"};
	}elsif($atoma =~/(.*)\?/)
	{
		my $res1=$residuea;
		$residuea=$residueb;
		$residueb=$res1;
		$residueTypea =$residues{$residuea}->{"residueType"};
		$residueTypeb =$residues{$residueb}->{"residueType"};
		if(!exists $connections{$residueTypea}->{$residueTypeb}->{"bond"}->[0]->{"energyGroup"}){
			smog_quit("Connection not defined for resTypes $residueTypea-$residueTypeb (residues $residuea $residueb)");
		}
		return $connections{$residueTypea}->{$residueTypeb}->{"bond"}->[0]->{"energyGroup"};
	}
}

sub parseBonds {
	$data = $xml->XMLin($bondxml,KeyAttr=>{function=>"name"},ForceArray=>1);
	
	## Obtain bond handle
	my @interHandle = @{$data->{"bonds"}->[0]->{"bond"}};
	
	## Loop over bonds, save in $interaction{bonds}{typeA}{typeB} = func info.
	my $counter = 0;
	foreach my $inter(@interHandle)
	{
		my ($typeA,$typeB) = @{$inter->{"bType"}}[0..1];
 		$bondtypesused{$typeA}=1;
 		$bondtypesused{$typeB}=1;
		my $func = $inter->{"func"};
		&checkFunction($func);
		if(exists $interactions->{"bonds"}->{$typeA}->{$typeB} || 
	               exists $interactions->{"bonds"}->{$typeA}->{$typeB}){
			smog_quit ("bond type between bType $typeA and bType $typeB defined more than once. Check .b file.");
		}
		$interactions->{"bonds"}->{$typeA}->{$typeB} = $func;
		$interactions->{"bonds"}->{$typeB}->{$typeA} = $func;
		$funcTable{"bonds"}->{$func} = $counter;
		$funcTableRev{"bonds"}->{$counter} = $func;
		$counter++;
	}
	
	## Obtain dihedral handle
	@interHandle = @{$data->{"dihedrals"}->[0]->{"dihedral"}};
	
	## Loop over dihedrals, save $interaction{dihedrals}{typeA}{typeB}{typeC}{typeD}
	## = function info. 
	$counter=0;
	foreach my $inter(@interHandle)
	{
		my ($typeA,$typeB,$typeC,$typeD) = @{$inter->{"bType"}}[0..3];

		foreach my $name($typeA, $typeB, $typeC, $typeD){
 			$bondtypesused{$name}=1;
		}

		my $func = $inter->{"func"};
		&checkFunction($func);

		my $eG;
		if(exists $inter->{"energyGroup"})
		{	
			$eG = $inter->{"energyGroup"};
		}
		
		my $keyString = "$typeA-$typeB-$typeC-$typeD";
		if(exists $interactions->{"dihedrals"}->{$eG}->{$keyString}){
			smog_quit ("dihedral type between bTypes $typeA-$typeB-$typeC-$typeD and energy group $eG defined more than once. Check .b file.");
		}
		$interactions->{"dihedrals"}->{$eG}->{$keyString} = $func;
		
		$keyString = "$typeD-$typeC-$typeB-$typeA";
	
		$interactions->{"dihedrals"}->{$eG}->{$keyString} = $func;
		
		$funcTable{"dihedrals"}->{$eG}->{$func} = $counter;
		$funcTableRev{"dihedrals"}->{$eG}->{$counter} = $func;
		
		$eGTable{$counter} = $eG;
		$eGRevTable{$eG} = $counter;
		$counter++;
	}

		
	## Obtain improper handle
	if(exists $data->{"impropers"}->[0]->{"improper"}){
		@interHandle = @{$data->{"impropers"}->[0]->{"improper"}};
		
		## Loop over dihedrals, save $interaction{dihedrals}{typeA}{typeB}{typeC}{typeD}
		## = function info. 
		$counter=0;
		foreach my $inter(@interHandle)
		{
			my ($typeA,$typeB,$typeC,$typeD) = @{$inter->{"bType"}}[0..3];

			foreach my $name($typeA, $typeB, $typeC, $typeD){
 				$bondtypesused{$name}=1;
			}

			my $func = $inter->{"func"};
			&checkFunction($func);
			my $keyString = "$typeA-$typeB-$typeC-$typeD";
			if(exists $interactions->{"impropers"}->{$keyString}){
				smog_quit ("improper type between bTypes $typeA-$typeB-$typeC-$typeD defined more than once. Check .b file.");
			}
			$interactions->{"impropers"}->{$keyString} = $func;
			
			$keyString = "$typeD-$typeC-$typeB-$typeA";
		
			$interactions->{"impropers"}->{$keyString} = $func;
			
			$funcTable{"impropers"}->{$func} = $counter;
			$funcTableRev{"impropers"}->{$counter} = $func;
			
			
			$counter++;
		}
	}	
	## Obtain angles handle
	@interHandle = @{$data->{"angles"}->[0]->{"angle"}};
	## Loop over angles, save $interaction{angles}{A}{B}{C} = function info.
	$counter = 0;
	foreach my $inter(@interHandle)
	{
		my ($typeA,$typeB,$typeC) = @{$inter->{"bType"}}[0..2];

		foreach my $name($typeA, $typeB, $typeC){
 			$bondtypesused{$name}=1;
		}

		my $func = $inter->{"func"};
		&checkFunction($func);
		my $keyString = "$typeA-$typeB-$typeC";
		## NOTE THE ORDER OF CENTRAL TYPE LISTED IN XML FILE MATTERS ##
		if(exists $interactions->{"angles"}->{$keyString}){
			smog_quit ("bond angle type between bTypes $typeA-$typeB-$typeC defined more than once. Check .b file.");
		}
		$interactions->{"angles"}->{$keyString} = $func;
		$keyString = "$typeC-$typeB-$typeA";
		$interactions->{"angles"}->{$keyString} = $func;
		$funcTable{"angles"}->{$func} = $counter;
		$funcTableRev{"angles"}->{$counter} = $func;
		$counter++;
	}

}

## PARSE NONBOND FILE ##
sub parseNonBonds {
	$data = $xml->XMLin($nbondxml,KeyAttr => ['name'],ForceArray=>1);

	my @interHandle;
	## Obtain default options ##
	if(exists $data->{"defaults"}){
		# defaults found in the templates, so use the values
		@interHandle = @{$data->{"defaults"}};
	}
	if(exists $interHandle[0]->{"gen-pairs"}) {
		if($interHandle[0]->{"gen-pairs"}==0){
			$interactions->{"gen-pairs"}="no";
		}else{
			$interactions->{"gen-pairs"}="yes";
		}
	}else{
		print "gen-pairs not defined in templates.  Will assume no.\n";
		$interactions->{"gen-pairs"}="no";
	}

	if(exists $interHandle[0]->{"nbfunc"}) {
		$interactions->{"nbfunc"} = $interHandle[0]->{"nbfunc"};
	}else{
		print "nbfunc not defined in templates.  Will assume a value of 1.\n";
		$interactions->{"nbfunc"}=1;
	}
	if(exists $interHandle[0]->{"gmx-combination-rule"}) {
		$interactions->{"gmx-combination-rule"} = $interHandle[0]->{"gmx-combination-rule"};
	}else{
		print "gmx-combination-rule not defined in templates.  Will assume a value of 1.\n";
		$interactions->{"gmx-combination-rule"}=1;
	}

	if(exists $interHandle[0]->{"fudgeLJ"}) {
		$interactions->{"fudgeLJ"} = $interHandle[0]->{"fudgeLJ"};
		if($interactions->{"fudgeLJ"} < 0 ){
			my $tmp=$interactions->{"fudgeLJ"};
			smog_quit("fudgeLJ must be greater than, or equal to, 0.  Found $tmp.");
		}

		if($interactions->{"gen-pairs"} eq "no"){
			smog_quit("fudgeLJ has no effect when gen-pairs=no.");
		}
	}else{
		print "fudgeLJ not defined in templates. Will assume a value of 1.\n";
		# we are setting to -1, to indicate that we don't want to write it later.
		$interactions->{"fudgeLJ"}=-1;
	}

	if(exists $interHandle[0]->{"fudgeQQ"}) {
		$interactions->{"fudgeQQ"} = $interHandle[0]->{"fudgeQQ"};
		if($interactions->{"fudgeQQ"} < 0 ){
			my $tmp=$interactions->{"fudgeQQ"};
			smog_quit("fudgeQQ must be greater than, or equal to, 0.  Found $tmp.");
		}
		if($interactions->{"fudgeLJ"}==-1){
			# since we gave a value for fudgeQQ and not fudgeLJ, we need to set the latter to 1, 
			# in order to write it properly, later
			$interactions->{"fudgeLJ"}=1;
			if($interactions->{"gen-pairs"} ne "no"){
				print "\nNote: Since fudgeQQ was set and gen-pairs=no, fudgeLJ will be written to the top file.\nHowever, this parameter has no effect, and is only written as a placeholder\n\n.";
			}
		}
	}else{
		print "fudgeQQ not defined in templates. Will assume a value of 1.\n";
		# we are setting to -1, to indicate that we don't want to write it later.
		$interactions->{"fudgeQQ"}=-1;
	}

	## Loop over nonbonds, save in $interaction{nonbonds}{typeA} = func info.
	@interHandle = @{$data->{"nonbond"}};
	my $counter = 0;
	foreach my $inter(@interHandle)
	{
		my $typeA = $inter->{"nbType"}->[0];
		my $func;
		if($inter->{"mass"} <=0){
			my $M=$inter->{"mass"};
			smog_quit("The mass of each atom must be positive. $M given for nbType=$typeA.");
		}

		if($interactions->{"gmx-combination-rule"}==1){
			if(!exists $inter->{"c6"} || !exists $inter->{"c12"}){
				smog_quit ("nonbonded parameters issue for nbType $typeA. c6, c12 must be provided when using gmx-combination-rule=1. Check .nb file.");
			}
			if(exists $inter->{"sigma"} || exists $inter->{"epsilon"}){
				smog_quit ("nonbonded parameters issue for nbType $typeA. sigma and epsilon can not be provided when using gmx-combination-rule=1. Check .nb file.");
			}	
			$func = {"mass" => $inter->{"mass"},"charge" => $inter->{"charge"},
			"ptype"=>$inter->{"ptype"},"c6"=>$inter->{"c6"},"c12"=>$inter->{"c12"}};
		}

		if($interactions->{"gmx-combination-rule"}==2){
			if(!exists $inter->{"sigma"} || !exists $inter->{"epsilon"}){
				smog_quit ("nonbonded parameters issue for nbType $typeA. sigma and epsilon must be provided when using gmx-combination-rule=2. Check .nb file.");
			}
			if(exists $inter->{"c6"} || exists $inter->{"c12"}){
				smog_quit ("nonbonded parameters issue for nbType $typeA. c6, c12 sigma and epsilon can not be provided when using gmx-combination-rule=2. Check .nb file.");
			}	
			$func = {"mass" => $inter->{"mass"},"charge" => $inter->{"charge"},
			"ptype"=>$inter->{"ptype"},"sigma"=>$inter->{"sigma"},"epsilon"=>$inter->{"epsilon"}};
		}

		if(exists $interactions->{"nonbonds"}->{$typeA}){
			smog_quit ("nonbonded parameters defined multiple times for nbType $typeA. Check .nb file.");
		}
		$interactions->{"nonbonds"}->{$typeA} = $func;
		$funcTable{"nonbonds"}->{$func} = $counter;
		$funcTableRev{"nonbonds"}->{$counter} = $func;
		$counter++;
	}
	
	@interHandle = @{$data->{"contact"}};
	## Loop over contacts, save in $interaction{contacts}{typeA}{typeB} = func info.
	$counter = 0;
	foreach my $inter(@interHandle)
	{
		my $nbtype;my $type;
	 	$nbtype = $inter->{"pairType"};$type = "pairType";
		my $typeA = $inter->{$type}->[0];
		my $typeB = $inter->{$type}->[1];
	 	$pairtypesused{$typeA}=1;
	 	$pairtypesused{$typeB}=1;
		my $func = $inter->{"func"};
		&checkFunction($func);
	 	my $cG = $inter->{"contactGroup"};
		if(exists $interactions->{"contacts"}->{"func"}->{$typeA}->{$typeB}){
			smog_quit ("contact parameters defined multiple times for types $typeA-$typeB. Check .nb file.");
		}
		$interactions->{"contacts"}->{"func"}->{$typeA}->{$typeB} = $func;
		$interactions->{"contacts"}->{"func"}->{$typeB}->{$typeA} = $func;
		$interactions->{"contacts"}->{"contactGroup"}->{$typeA}->{$typeB} = $cG;
		$interactions->{"contacts"}->{"contactGroup"}->{$typeB}->{$typeA} = $cG;
		$funcTable{"contacts"}->{"func"}->{$func} = $counter;
		$funcTable{"contacts"}->{"contactGroup"}->{$cG} = $counter;
		$funcTableRev{"contacts"}->{$counter}->{"func"} = $func;
		$funcTableRev{"contacts"}->{$counter}->{"contactGroup"} = $cG;
		$counter++;
	}
	
	

}

######################
## BOND FUNCTIONALS ##
######################
sub createBondFunctionals {
foreach my $res (keys %residues)
{
	$residueHandle = $residues{$res};
	%bondHandle = %{$residueHandle->{"bonds"}};
	my @inputString; my @functionString;
	my %adjList; 

# check that every atom in an improper is actually present in the residue
	@improperHandle = @{$residueHandle->{"impropers"}};

	foreach my $imp(@improperHandle){
		my $ats = $imp->[0] . "-" . $imp->[1] . "-" . $imp->[2] . "-" . $imp->[3];
		for (my $ii=0;$ii<4;$ii++){
			if(!exists $residueHandle->{"atoms"}->{$imp->[$ii]}){
				smog_quit ("$imp->[$ii] doesn't exists in $res but it is used in improper: $ats");
			}
		}
	}

	foreach my $bondInfo(keys %bondHandle)
	{
		my ($atomA,$atomB) = $bondInfo =~ /(.*)\-(.*)/;
  ## Check if atoms exists in declaration ##
  	    if(!exists $residueHandle->{"atoms"}->{"$atomA"}) 
		{smog_quit ("$atomA doesn't exists in $res but a bond was defined $bondInfo"); }
		if(!exists $residueHandle->{"atoms"}->{"$atomB"}) 
		{smog_quit ("$atomB doesn't exists in $res but a bond was defined $bondInfo"); }

		my $typeA = $residueHandle->{"atoms"}->{"$atomA"}->{"bType"};
		my $typeB = $residueHandle->{"atoms"}->{"$atomB"}->{"bType"};
		## WILD CARD MATCHING CONDITIONALS ##

		## If both bond types exists ##
		my $funct="";
		if( exists $interactions->{"bonds"}->{$typeA}->{$typeB})
		{$funct = $interactions->{"bonds"}->{$typeA}->{$typeB};}
			
		elsif ($typeA ne $typeB && (exists $interactions->{"bonds"}->{$typeA}->{"*"} 
                                 && exists $interactions->{"bonds"}->{$typeB}->{"*"})){
			smog_quit ("Wildcard conflict in bonds $typeA-$typeB. 
			Both $typeA-\* and $typeB-\* are defined in .b file. Can not unambiguously assign a function...\n\n");
 		}
		## If typeA exists while TypeB is a wildcard ##
		elsif (exists $interactions->{"bonds"}->{$typeA}->{"*"})
		{$funct = $interactions->{"bonds"}->{$typeA}->{"*"};}

		## If typeB exists while TypeA is a wildcard ##
		elsif (exists $interactions->{"bonds"}->{$typeB}->{"*"})
		{$funct = $interactions->{"bonds"}->{$typeB}->{"*"};}
	
		if(!defined $funct || $funct eq ""){
			if(exists $interactions->{"bonds"}->{"*"}->{"*"})
            		{$funct = $interactions->{"bonds"}->{"*"}->{"*"};}
		     	else{
			smog_quit ("Unable to unambiguously assign bond types to all bonds in a residue\n Offending btypes are $typeA $typeB");
			}
		}
		my $indexA = $residueHandle->{"atoms"}->{"$atomA"}->{"index"};
		my $indexB = $residueHandle->{"atoms"}->{"$atomB"}->{"index"};
		push(@inputString,"$indexA-$indexB"); ## MIGHT CHANGE
		push(@functionString,"$funct");
		
		push(@{$adjList{$atomA}},$atomB);
		push(@{$adjList{$atomB}},$atomA);
	   }
	$bondFunctionals{$res} = {"bonds"=>\@inputString, "functions"=>\@functionString}; 
	$dihedralAdjList{$res} = \%adjList;
	}
}

sub getContactFunctionals
{
	my($typeA,$typeB) = @_;
	## WILD CARD MATCHING CONDITIONALS ##
	my $funct=""; my $cG = "";
	my $assigned=0; 
	## If contact type is specifically defined ##
	if(exists $interactions->{"contacts"}->{"func"}->{$typeA}
		&& exists $interactions->{"contacts"}->{"func"}->{$typeA}->{$typeB})
	{
		$funct = $interactions->{"contacts"}->{"func"}->{$typeA}->{$typeB};
		$cG = $interactions->{"contacts"}->{"contactGroup"}->{$typeA}->{$typeB};
		$assigned++;
	}else{
	
		## If typeA matches and TypeB matches a wildcard ##
		if (exists $interactions->{"contacts"}->{"func"}->{$typeA}
			&& (exists $interactions->{"contacts"}->{"func"}->{$typeA}->{"*"}))
		{
			$funct = $interactions->{"contacts"}->{"func"}->{$typeA}->{"*"};
			$cG = $interactions->{"contacts"}->{"contactGroup"}->{$typeA}->{"*"};
			$assigned++;
		}
		if($typeB ne $typeA){
			# only check typeB if it is different from typeA
	
			## If typeB matches and TypeA matches a wildcard ##
			if((exists $interactions->{"contacts"}->{"func"}->{$typeB}) &&
					exists $interactions->{"contacts"}->{"func"}->{$typeB}->{"*"})
			{
				my $functB = $interactions->{"contacts"}->{"func"}->{$typeB}->{"*"};
				my $cGB = $interactions->{"contacts"}->{"contactGroup"}->{$typeB}->{"*"};
				if($assigned ==0) {
					#typeA did not match earlier, so set B.
					$cG=$cGB;
					$funct=$functB;
					$assigned++;
				}else{
					#already matched typeA, and now typeB.  If the interactions are identical, there is no ambiguity.  If they are different function types, or groups, then give an error.
					if($cG ne $cGB || $funct ne $functB){
						smog_quit("Can\'t unambiguously assign a contact interaction between atoms of type $typeA and $typeB. Matched the following definitions equally well:\n\nfunc contactGroup\n$funct $cG\n$functB $cGB\n\nSee .nb for contact group definitions.\n");
					}
				}
			}
		}	
		## If nothing has matched, and double wildcard interaction is defined ## 
		if($assigned==0 && exists $interactions->{"contacts"}->{"func"}->{"*"}->{"*"}){

			$funct = $interactions->{"contacts"}->{"func"}->{"*"}->{"*"};
			$cG = $interactions->{"contacts"}->{"contactGroup"}->{"*"}->{"*"};
			$assigned++;
		}
	}
	
	if($assigned ==0){
		smog_quit("Can\'t find a contact interaction that matches atomtype pair $typeA and $typeB.  See .nb for contact group definitions.\n");
	}elsif($assigned >1){
		smog_quit("Internal Error 10. Please contact SMOG team to report.")
	}

	return ($funct,$cG);
}


###################################
## DIHEDRALS & ANGLE FUNCTIONALS ##
###################################

## FINDING DIHEDRALS & ANGLES ##
## adjListTraversalHelper & adjListTraversal
## Takes a adjacency list of bonds and finds the dihedrals
## by traversing through the adj list (graph traversal).
## After finding four (3 for angles) consecutive atoms a 
## dihedralString (angleString) is created. 
## Post pruning is done to remove duplicate dihedrals/angles.
##

## Traverses through an adj list; works with adjListTraversal 
sub adjListTraversalHelper
{
	my($listHandle,$atomParent,$visitedString,$diheStringList,$angleStringList,$visitedList,$counter) = @_;
	my $limitCounter = $counter;
 	## End reached don't need to search further
 	if($counter >=4){return;}
 	## Given an atom loop through all the atoms it is bonded to
 	foreach my $atomIn(@{$listHandle->{$atomParent}})
 	{
		## If bond was already considered from atomIn goto next
		if(exists($visitedList->{$atomIn})){next;}
		## Stitch dihedral string with atomIn
		my $visitedStringNew = "$visitedString-$atomIn";
		## Update counters, visitedList hash
		$counter=$limitCounter+1;
		my %sendHash = %{$visitedList};
		$sendHash{$atomIn} = 1;	
		## Traverse through the child of atomIn
		adjListTraversalHelper($listHandle,$atomIn,$visitedStringNew,$diheStringList,$angleStringList,\%sendHash,$counter);
	
		## 3 atom count is reached, creating an angle
		## add to angleStringList; continue
		if($counter == 3) {push(@{$angleStringList},$visitedStringNew);$counter--;}
	
		## 4 atom count is reached, creating a dihedral;
		## add to the diheStringList; exit out of traversing this branch
		if($counter == 4) {push(@{$diheStringList}, $visitedStringNew);$counter--;}
	
		#if($counter == 4) {push(@{$diheStringList}, $visitedStringNew);last;}
	}
}

## Traverse through the adjacency list of bonds
sub adjListTraversal
{
	#print "START ROUTINE\n";
	my($listHandle) = @_;
	my @diheStringList; ## List of Dihedral atoms
	my @angleStringList; ## List of Angle atoms
	my @oneFourStringList; ## List of 1-4 atoms

	## Loop through all the atoms in the residue
	foreach my $atomOut(keys(%{$listHandle}))
	{
       		my %visitedList;
		my $visitedString;
   		## Add the current atom as the first atom in the dihedral
       		$visitedString="$atomOut";
       		$visitedList{$visitedString}=1; ## Set visited flag
       		adjListTraversalHelper($listHandle,$atomOut,$visitedString,\@diheStringList,\@angleStringList,\%visitedList,1); ## Traverse through bond graph
	}

	## Remove duplicate dihedrals and angles
 	my @uniqueD = ();
	my %seen  = ();
 	my @uniqueA = ();
 	foreach my $elem (@angleStringList)
 	{
    		my @orgsplit = split('-', $elem);
    		my $testKey = join('-', reverse(@orgsplit));
    		if(exists $seen{$testKey}){next;}
        	$seen{$elem} = 1;
       		push (@uniqueA,$elem)
 	}
	%seen  = ();

 	foreach my $elem (@diheStringList)
 	{
    		my @orgsplit = split('-', $elem);
    		my $testKey = join('-', reverse(@orgsplit));
    		if(exists $seen{$testKey}){next;}
        	$seen{$elem} = 1;
       		push (@uniqueD,$elem);
 	}

 	## Reset diheStringList, angleStringList after removing duplicates
 	@diheStringList = @uniqueD;
 	@angleStringList = @uniqueA;
 	return (\@diheStringList,\@angleStringList);
}

sub createDihedralAngleFunctionals {
## CALCULATE DIHEDRALS AND ANGLES FOR ALL RESIDUES
## AND MATCH IT WITH THE SPECIFIC ATOM TYPE FUNCTIONS
##

	foreach my $res(keys %dihedralAdjList)
	{
		my ($dihes,$angles) = adjListTraversal($dihedralAdjList{$res});
		## MATCH DIHEDRALS WITH FUNCTIONS ##
		my $diheHandle = $interactions->{"dihedrals"};
		my @allAngles; 
		my @allDihes;
		my @allAnglesFunct; 
		my @allDihesFunct;
		my $funct;
		 foreach my $dihs(@{$dihes})
		{
			my @atoms = split("-",$dihs);
			my $atomHandle = $residues{$res}->{"atoms"};
			my $a = $atomHandle->{$atoms[0]}->{"bType"};
			my $b = $atomHandle->{$atoms[1]}->{"bType"};
			my $c = $atomHandle->{$atoms[2]}->{"bType"};
			my $d = $atomHandle->{$atoms[3]}->{"bType"};
			
			my $eG = getEnergyGroup($res,$res,$atoms[1],$atoms[2]);
			
			## WILD CARD MATCHING CONDITIONALS ##
			my $matchScore = 0; 
			my $saveScore = 0;
			my $matchScoreCount=0; 
			my $symmatch=0;
			my $Nd=0;
			foreach my $matches(keys %{$diheHandle->{$eG}})
			{
			$Nd++;
				$matchScore = 0;
				my ($aM,$bM,$cM,$dM) = split("-",$matches);
				unless(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
					|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
					|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)
					|| ($d !~ /\Q$dM\E/ && $dM !~ /\Q*\E/)){
					if($a =~ /\Q$aM\E/) {$matchScore+=2;} else {$matchScore+=1;}
					if($b =~ /\Q$bM\E/) {$matchScore+=2;} else {$matchScore+=1;}
					if($c =~ /\Q$cM\E/) {$matchScore+=2;} else {$matchScore+=1;}
					if($d =~ /\Q$dM\E/) {$matchScore+=2;} else {$matchScore+=1;}
					if($matchScore >= $saveScore){
	
						if(($aM eq $dM and $bM eq $cM) || ($aM eq $bM and $bM eq $cM and $cM eq $dM)){
							$symmatch=1;
						}else{
							$symmatch=0;
						}
						## this to make sure that the highest scoring angle is unique
						if($matchScore == $saveScore){
							if($saveScore != 0){
							$matchScoreCount++;
							}
						}else{
							$matchScoreCount=0;
						}
						$saveScore = $matchScore;
						$funct = $diheHandle->{$eG}->{$matches};
					}
				}
			}
	
			if($Nd ==0){
				smog_quit ("energy group $eG is used in .bif file, but it is not defined in .b file.");
			}
			
			my $sym=0;
			if(($a eq $d and $b eq $c) || ($a eq $b and $b eq $c and $c eq $d)){
				$sym=1;
			}
			if(($symmatch ==0 && $sym == 1 && $matchScoreCount > 1)  || ($symmatch ==0 && $sym == 0 && $matchScoreCount > 0) || ($symmatch ==1 && $sym == 0 && $matchScoreCount > 0) || ($symmatch ==1 && $sym == 1 && $matchScoreCount > 0)){
	
				smog_quit ("$symmatch  $sym $matchScoreCount Multiple possible angles match $a-$b-$c-$d, and energyGroup $eG equally well. Can not determine function based on .b file.");
			}
			if($saveScore == 0){
				smog_quit ("Dihedral Angle between bTypes $a-$b-$c-$d and energyGroup $eG: Unable to match to a function in .b file.");
			}
	
			my $indexA = $residues{$res}->{"atoms"}->{$atoms[0]}->{"index"};
			my $indexB = $residues{$res}->{"atoms"}->{$atoms[1]}->{"index"};
			my $indexC = $residues{$res}->{"atoms"}->{$atoms[2]}->{"index"};
			my $indexD = $residues{$res}->{"atoms"}->{$atoms[3]}->{"index"};
			my $inputString = "$indexA-$indexB-$indexC-$indexD";
			push(@allDihes,$inputString);
			push(@allDihesFunct,$funct);
		
		}
		$dihedralFunctionals{$res} = {"dihedrals"=>\@allDihes,"functions"=>\@allDihesFunct};
	
	#######################
	## ANGLE FUNCTIONALS ##
	#######################
	
		## MATCH ANGLES WITH FUNCTIONS ##
		my $angHandle = $interactions->{"angles"};
		foreach my $angs(@{$angles})
		{
			my @atoms = split("-",$angs);
			my $atomHandle = $residues{$res}->{"atoms"};
			my $a = $atomHandle->{$atoms[0]}->{"bType"};
			my $b = $atomHandle->{$atoms[1]}->{"bType"};
			my $c = $atomHandle->{$atoms[2]}->{"bType"};
			#print $angHandle->{$atoms[0]}->{$atoms[1]}->{$atoms[2]},"\n";
			
			## WILD CARD MATCHING CONDITIONALS ##
			my $matchScore = 0; my $saveScore = 0; my $matchScoreCount=0; my $symmatch=0;
			foreach my $matches(keys %{$angHandle})
			{
				$matchScore = 0;
				my ($aM,$bM,$cM) = split("-",$matches);
				unless(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
					|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
					|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)){
				if("$a" eq "$aM") {$matchScore+=2;} else {$matchScore+=1;}
				if("$b" eq "$bM") {$matchScore+=2;} else {$matchScore+=1;}
				if("$c" eq "$cM") {$matchScore+=2;} else {$matchScore+=1;}
				if($matchScore >= $saveScore)
					{
					if($aM eq $cM || ($aM eq $bM and $bM eq $cM)){
						$symmatch=1;
					}else{
						$symmatch=0;
					}
					## this to make sure that the highest scoring angle is unique
					if($matchScore == $saveScore){
						if($saveScore != 0){
						$matchScoreCount++;
						}
					}else{
						$matchScoreCount=0;
					}
					$saveScore = $matchScore;$funct = $angHandle->{$matches};
				}
			    }
			}
			my $sym=0;
			if($a eq $c || ($a eq $b and $b eq $c)){
				$sym=1;
			}
			if(($symmatch ==0 && $sym == 1 && $matchScoreCount != 1)  || ($symmatch ==0 && $sym == 0 && $matchScoreCount != 0) || ($symmatch ==1 && $sym == 0 && $matchScoreCount != 0) || ($symmatch ==1 && $sym == 1 && $matchScoreCount != 0)){
				smog_quit ("Multiple possible angles match $a-$b-$c equally well. Unclear assignment of function type");
			}
			my $indexA = $residues{$res}->{"atoms"}->{$atoms[0]}->{"index"};
			my $indexB = $residues{$res}->{"atoms"}->{$atoms[1]}->{"index"};
			my $indexC = $residues{$res}->{"atoms"}->{$atoms[2]}->{"index"};
			my $inputString = "$indexA-$indexB-$indexC";
			push(@allAngles,$inputString);
			push(@allAnglesFunct,$funct);
		}
		$angleFunctionals{$res} = {"angles"=>\@allAngles,"functions"=>\@allAnglesFunct};
	}

}


1;
