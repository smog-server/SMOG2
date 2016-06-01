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
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#########################################################################################

##############################################################################
# PDB_Bonded: A set of routines for parsing the PDB file and generating 
# bonded interaction information, residue and coordinate info.
##############################################################################
package PDB_Bonded;

use templateParser;
use setRatios;
use strict;
use warnings;
use Exporter;
use PDL; ## LOAD PDL MODULE
use Storable qw(dclone);

use Time::HiRes qw( gettimeofday);

## DELETE LATER ##
#use Devel::Size qw(size total_size);

## DECLEARATION TO SHAR DATA STRUCTURES ##
our @ISA = 'Exporter';
our @EXPORT = 
qw(%eGTable $energyGroups $interactionThreshold %fTypes %residues $termRatios %allAtoms parseCONTACT $contactPDL catPDL $totalAtoms returnFunction intToFunc funcToInt %bondFunctionals %AngleData %DihedralData %BondData %resPDL %bondPDL %dihedralFunctionals %angleFunctionals setInputFileName parseBif parseSif parseBonds createBondFunctionals createDihedralAngleFunctionals parseNonBonds getContactFunctionals $contactSettings $interactions clearPDBMemory clearBifMemory parsePDBATOMS);

my @vector;
my $coorPDL;
my %results;
my %residuePDL=();

our %tempPDL = ();
our %resPDL;
our %bondPDL;
our %indexMap;
my $angToNano = 0.1;


our %AngleData;
our %DihedralData;
our %BondData;

my @consecResidues;


our $totalAtoms;
our $contactPDL;
our %allAtoms;
our %allAtomsBackup;

our %extContacts;

###########################
## CLEAR VARIABLE MEMORY ##
###########################
sub clearPDBMemory {
undef %tempPDL;undef %resPDL;undef %bondPDL;
undef %AngleData;undef %DihedralData;
undef %BondData;undef $totalAtoms;
#undef $contactPDL; 
undef %indexMap;
#for coarse graining contact maps we need to know the atomNum<->resNum mapping
%allAtomsBackup = %{ dclone (\%allAtoms) };
undef %allAtoms; 
}


####################################################################
# parseExternalContacts($filename)
####################################################################
sub parseExternalContacts
{
  my($file) = @_;
}


####################################################################
# parsePDBATOMS
####################################################################
sub parsePDBATOMS
{
	
	my ($fileName,$CGenabled) = @_;
	
	## INTERNAL VARIABLES ##
	my $counter = 0;
	my @temp; my @connResA; my @connResB; my @union;
	my @tempBond;
	my @consecResidues;
	my $x;my $y;my $z;
	my $residue; my $interiorResidue; my $atom;my $atomSerial;
	my $atomsInRes; my $lineEnd;
	my $i; my $putIndex=0; my $strLength;
	my $resType;
	my $angH; my $diheH;
	my $bondStrA;my $bondStrB;my $typeA;my $typeB;
	my $endFlag=0; my $headFlag=1;my $outLength;
	$totalAtoms = 0;my $nbType;my $residueType; my $pairType;
	my $atomCounter=0;my $singleFlag = 1;
	my $chainNumber = 0;my $linkFlag = 0;
	my $residueIndex=1;
	my $interiorPdbResidueIndex=0;
	my $lineNumber = 0;
	my %connectedatom;
	my $lastchainstart=0;
	my $endfound=0;
	my $residueSerial=0;
	## OPEN .PDB FILE ##

	unless (open(PDBFILE, $fileName)) {
		smog_quit ("Cannot read from '$fileName'.");
	}

	my $lastresindex="null";
	 ## LOOP THROUGH EACH LINE ##
	while(my $record = <PDBFILE>)
	{
		$lineNumber++;
	
		my @impAtoms = ();
		## PARSE BOND LINES ##
	
	if($record =~ m/^COMMENT/){
		next;
	# make sure BOND appears after END
	}elsif($record !~ m/^BOND/ && $endfound ==1){
		smog_quit("PDB format issue: Only user-defined bonds given by BOND, or COMMENT lines, may be listed after END.");
	}

 	if($record =~ m/^BOND/)
 	{

		if($CGenabled==1){
			smog_quit("User-defined bonds, via BOND declaration, are not supported with Coarse-Grained models. Remove BOND lines and try again.");
			next;
  		}elsif($endfound ==0){
   			smog_quit("PDB format issue: User-defined bonds given by BOND should be listed immediately after END.");
  		}

    		chomp($record);
   
		my @TMP = split(/\s+/,$record);
		if(@TMP <= 5){
			smog_quit("Directive BOND must have 5 arguments. Offending line:\n$record");
		}
		my($trig,$chaina,$atoma,$chainb,$atomb,$eG) = split(/\s+/,$record);
		
		#internally, chains are indexed 0,1...
		$chaina--;
		$chainb--;
		if(!exists $indexMap{"$chaina-$atoma"}){
			my $chaina1=$chaina+1;
			smog_quit("Can not find atom $atoma in chain $chaina1");
		}
		if(!exists $indexMap{"$chainb-$atomb"}){
			my $chainb1=$chainb+1;
			smog_quit("Can not find atom $atomb in chain $chainb1");
		}
		my $idxA = $indexMap{"$chaina-$atoma"};
		my $idxB = $indexMap{"$chainb-$atomb"};
		my $resA = $allAtoms{$idxA}->[5];
		my $resB = $allAtoms{$idxB}->[5];
		my $atomA = $allAtoms{$idxA}->[3];
		my $atomB = $allAtoms{$idxB}->[3];
		my $resAIdx = $allAtoms{$idxA}->[2];
		my $resBIdx = $allAtoms{$idxB}->[2];
		
		my $sizeA = scalar(keys %{$residues{$resA}->{"atoms"}});
		my $union;
		$union=($tempPDL{$resA}->{$resAIdx})->glue(1,$tempPDL{$resB}->{$resBIdx});
		print "\nNOTE:";
		my $chaina1=$chaina+1;
		my $chainb1=$chainb+1;
		print "Generating user-specified bonded interaction between chain-atom pair $chaina1-$atoma,$chainb1-$atomb.\nWill assign to energy group $eG.\n";
		if(exists $connectedatom{$idxA}){ 
			smog_quit("Currently, including a BOND with an atom that is also declared in \"connections\" is not supported.\nOffending atom ($atomA, in $resA$resAIdx) and line:$record");
		}
		if(exists $connectedatom{$idxB}){ 
			smog_quit("Currently, including a BOND with an atom that is also declared in \"connections\" is not supported.\nOffending atom ($atomB, in $resB$resBIdx) and line:$record");
		}
		$bondPDL{$counter}=$union;
		## Check if improper directive is present ##
		if($record =~ m/IMPROPER/)
		{
			my($left,$right) = split(/IMPROPER/,$record);
			$right =~ s/^\s+|\s+$//g;
			@impAtoms = split(/\s+/,$right);
			print "IMPROPER DETECTED @impAtoms\n";
		}
		connCreateInteractionsSingleBOND([$resA,$resB],$sizeA,$counter,$atomA,$atomB,$resAIdx,$resBIdx,$eG,\@impAtoms); 
		$counter++;
		next;
	}

	## IF TER LINE  ##
	if($record =~ m/TER|END/)
	{
 		$lastresindex="null";
		$chainNumber++; ## INCREMENT CHAIN NUMBER ##
		## CREATE INTERACTION ##
		my $chainlength=$atomSerial-$lastchainstart;
		print "Building covalent geometry for chain $chainNumber\n";
        	my $connset=GenerateBondedGeometry(\@consecResidues,$counter,$chainNumber,$chainlength);
		my @connset=@{$connset};
		foreach my $I(@connset){
			my $T=$I+$lastchainstart+1;
			$connectedatom{$T}=1;
		}
	   	$bondPDL{$counter}=pdl(@union);
		@union = ();$counter++;
        	@consecResidues = ();
		$lastchainstart=$atomSerial;
		if($record =~ m/END/){$endfound=1;}
		next;
	} 
	
	## ONLY WORK WITH ATOM LINES ##
	if($record =~ m/ATOM/ || $record =~ m/HETATM/)
	{
        	$lineNumber--;
		$outLength = length($record);
	 	## OBTAIN RESIDUE NAME ##
		$residue = substr($record,17,4);
		$residue =~ s/^\s+|\s+$//g;
		if(!exists $residues{$residue}){smog_quit (" \"$residue\" doesn't exist in .bif. See line $lineNumber of PDB file.");}
		## if first iteration, save residueBackup, and use residues
		#if(exists $residueBackup{$residue}){
		if($CGenabled == 1){
			$atomsInRes = scalar(keys(%{$residueBackup{$residue}->{"atoms"}}));
		}else{
			$atomsInRes = scalar(keys(%{$residues{$residue}->{"atoms"}}));
		}
		my $atomsInBif=scalar(keys(%{$residues{$residue}->{"atoms"}}));
		if($atomsInBif != 1 && $CGenabled ==1)
                {
			smog_quit ("When using CG, each residue can only have one atom in the CG template. Check .bif definition for $residue");
		}
		my $atomsmatch=0;
	 	seek(PDBFILE, -$outLength, 1); # place the same line back onto the filehandle
		my $resname=$residue;
        	my $resindex = substr($record,22,5);
		if ($lastresindex ne "null" && $resindex-$lastresindex ==0 ){
			$lineNumber++;
			smog_quit("Extra atoms found in residue $resname. Check near line $lineNumber.");
		}elsif ($lastresindex ne "null" && $resindex-$lastresindex != 1){
			smog_quit("Non-sequential residue numbers ($lastresindex,$resindex) appear at line $lineNumber.");
		}
		$lastresindex=$resindex;
		my %uniqueAtom;
		$residueSerial++;
		for($i=0;$i<$atomsInRes;$i++)
		{
 			$lineNumber++;
			$record = <PDBFILE>;
			if($record !~ m/^ATOM|^HETATM/)
			{smog_quit("Expected ATOM or HETATM at line $lineNumber. Residue $residue might have been truncated at $lineNumber");}

			$interiorResidue = substr($record,17,4);
			$interiorResidue =~ s/^\s+|\s+$//g;
	   		$residue = substr($record,17,4);
        	        $residue =~ s/^\s+|\s+$//g;
	   		my $altlocator = substr($record,16,1);
			if($altlocator ne " "){
				smog_quit("Alternate location indicator found at line $lineNumber.  Alt. Loc. Indic. not supported by SMOG.");
			}

	        	if(!exists $residues{$residue}){smog_quit (" \"$residue\" doesn't exist in .bif. See line $lineNumber of PDB file.");}
            		$interiorPdbResidueIndex = substr($record,22,5);  
			if($resname ne $residue or $resindex ne $interiorPdbResidueIndex){
				$lineNumber--;
				my $missingatoms="";
				$uniqueAtom{$atom}=1;
			        if($CGenabled == 1){
					foreach my $atomcheck(keys (%{$residueBackup{$resname}->{"atoms"}})){
                        			if(!exists $uniqueAtom{$atomcheck}){
							$missingatoms=$missingatoms . "$atomcheck ";
						}
					}		
                		}else{
					foreach my $atomcheck(keys (%{$residues{$resname}->{"atoms"}})){
                        			if(!exists $uniqueAtom{$atomcheck}){
							$missingatoms=$missingatoms . "$atomcheck ";
						}
					}		
                		}
				smog_quit("It appears that a residue in the PDB file does not contain all of the atoms defined in the .bif file.\nOffending residue: $resname (ending at line $lineNumber).  Missing atoms: $missingatoms");	
			}

			$interiorPdbResidueIndex =~ s/^\s+|\s+$//g;
			unless($interiorPdbResidueIndex =~ /^\d+$/){;
				smog_quit ("Residue $residue$interiorPdbResidueIndex contains non integer value for the index, or an insertion code.");
			}

			## CHECK IF ALL ATOMS CONFORM TO BIF RESIDUE DECLARATION ##
			$atom = substr($record, 12, 4);
			$atom =~ s/^\s+|\s+$//g;

			if(exists $uniqueAtom{$atom})
			{
				smog_quit("$atom appears twice in $residue at line $lineNumber\n");i
			}
			else {
				$uniqueAtom{$atom}=1;
			}
				
			if($CGenabled == 0 && !exists $residues{$residue}->{"atoms"}->{$atom})
			{
				smog_quit ("\"$atom\" doesn't exist in .bif declaration of \"$residue\"");
			}
			
			## CHECK IF ATOM EXISTS IN MODEL ##
			if(!exists $residues{$residue}->{"atoms"}->{$atom}){next;}
			$atomsmatch++;
			$x = substr($record, 30, 8);
			$y = substr($record, 38, 8);
			$z = substr($record, 46, 8);
			$atomCounter++;
			$atomSerial=$atomCounter;
			$atomSerial =~ s/^\s+|\s+$//g;	
			
			$putIndex = $residues{$residue}->{"atoms"}->{$atom}->{"index"};
			$nbType = $residues{$residue}->{"atoms"}->{$atom}->{"nbType"};
			$pairType = $residues{$residue}->{"atoms"}->{$atom}->{"pairType"};
			$residueType = $residues{$residue}->{"residueType"};
			
			$allAtoms{$atomSerial}=[$nbType,$residueType,$residueIndex,$atom,$chainNumber,$residue,$x,$y,$z,$residueSerial,$pairType]; 
			my $pdbIndex;
			if($CGenabled==1){
				$pdbIndex = $interiorPdbResidueIndex;
			}else{
				$pdbIndex = substr($record,6,5);
			}
			$pdbIndex =~ s/^\s+|\s+$//g;
			if(exists $indexMap{"$chainNumber-$pdbIndex"}){
				my $chainID=$chainNumber+1;
				smog_quit("Atom/Residue numbers must be unique within each chain. Offending line:\n$record");
			}
			$indexMap{"$chainNumber-$pdbIndex"}=$atomSerial;
			$temp[$putIndex]=[$x,$y,$z,$atomSerial];
			$tempBond[$putIndex]=[$x,$y,$z,$atomSerial];
			$totalAtoms++;
		}



		if($atomsmatch != $atomsInBif){
			smog_quit ("Not all atoms in the bif appear in the PDB. See line $lineNumber");
		}

		if($i != $atomsInRes){smog_quit ("Total number of atoms of $residue doesn't match with .bif declaration");}
		## CONCAT RESIDUE ##
	  	@union = (@union,@temp);@temp=();
		push(@consecResidues,$residue);
		$headFlag = 0;		
        	$tempPDL{$residue}->{$residueIndex}=pdl(@tempBond);
		@tempBond = ();
		$residueIndex++;
				
		}else{smog_quit(" Expected ATOM or HETATM line at Line $lineNumber. Residue $residue might have been truncated at $lineNumber");}
		$record = "";
	
 	}
}

# returnFunction: Return the fType and directive field for a specified function
sub returnFunction
{
	my($funcString) = @_;
	my $addExclusions;
	if(!exists $fTypes{"$funcString"}){smog_quit ("$funcString is not a supported function type in SMOG");}
	if(!exists $functions->{$funcString}){smog_quit ("Function $funcString is being used, but is not defined in .sif file");}
	#Sometimes exclusions are not defined for contacts that go under other directives. Need to set it to zero.
	if(!exists $functions->{$funcString}->{"exclusions"}){ $addExclusions = 0; }
	else { $addExclusions = $functions->{$funcString}->{"exclusions"}; }
	return ($fTypes{"$funcString"},$functions->{$funcString}->{"directive"},$addExclusions);
}

##
# getAtomIndexInResidue: Return the index of an atom from storef indexing
sub getAtomIndexInResidue
{
	my($residue,$atom) = @_;
	if(!exists $residues{$residue}){smog_quit ("$residue wasn't defined in bif");}
	if(!exists $residues{$residue}->{"atoms"}->{$atom}){smog_quit ("$atom wasn't defined in $residue in the bif");}
	return $residues{$residue}->{"atoms"}->{$atom}->{"index"};
}

##
# getAtomBType: Return the bondType of an atom
sub getAtomBType
{
	my($residue,$atom) = @_;
	return $residues{$residue}->{"atoms"}->{$atom}->{"bType"};
}


sub singleCreateInteractions
{
	my($residue,$counter) = @_;

	## CREATE SINGLE INTERACTIONS ##
	my $union = $dihedralAdjList{$residue};	
	my ($diheH,$angH,$oneFour)=adjListTraversal($union); 

	## CREATE INTERACTION PDLs ##
	## ANGLES ##
	my @tempArr;
	foreach my $angs(@{$angH})
	{
		my($a,$b,$c) = split("-",$angs);
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
	 	($ia,$ta) = (getAtomIndexInResidue($residue,$a),getAtomBType($residue,$a));
		($ib,$tb) = (getAtomIndexInResidue($residue,$b),getAtomBType($residue,$b));
		($ic,$tc) = (getAtomIndexInResidue($residue,$c),getAtomBType($residue,$c));
   		my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
   		push(@tempArr,pdl($ia,$ib,$ic,$if));		
	}
		$AngleData{$counter} = cat(@tempArr);
		@tempArr = ();
				
	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
 		($ia,$ta) = (getAtomIndexInResidue($residue,$a),getAtomBType($residue,$a));
		($ib,$tb) = (getAtomIndexInResidue($residue,$b),getAtomBType($residue,$b));
		($ic,$tc) = (getAtomIndexInResidue($residue,$c),getAtomBType($residue,$c));
		($id,$td) = (getAtomIndexInResidue($residue,$d),getAtomBType($residue,$d));
				
		my $eG = getEnergyGroup($consecResidues[0],$consecResidues[1],$b,$c);
		my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);
		$eG = $eGRevTable{$eG};
		
		## [x,y,z,func,countDihedrals,energyGroup]
		push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);	
	}
	    
		## Manually add Improper dihedrals ##
 	my $resBIp = $residues{$residue}->{"impropers"};
 	foreach my $ips(@{$resBIp})
 	{
		my $ta;my $tb;my $tc;my $td;
   		if(! (defined $ips) ) {next;}
		my($ia,$ib,$ic,$id) = @{$ips};
		($ia,$ta) = (getAtomIndexInResidue($residue,$ia),getAtomBType($residue,$ia));
		($ib,$tb) = (getAtomIndexInResidue($residue,$ib),getAtomBType($residue,$ib));
		($ic,$tc) = (getAtomIndexInResidue($residue,$ic),getAtomBType($residue,$ic));
		($id,$td) = (getAtomIndexInResidue($residue,$id),getAtomBType($residue,$id));
		my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
	## [a,b,c,d,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
	push(@tempArr,[$ia,$ib,$ic,$id,$if,1,-1]);	
	}
		
	$DihedralData{$counter} = pdl(@tempArr);
	@tempArr = ();
}

sub connectivityHelper
{
	my($listHandle,$atomParent,$visitedList) = @_;
	## Given an atom loop through all the atoms it is bonded to
	my @newatoms=();
	foreach my $atomIn(@{$listHandle->{$atomParent}})
	{
	 	## If atom has not already considered add to list for next check
	   	if(!exists($visitedList->{$atomIn})){
	       		push(@newatoms,$atomIn)
		}
	}
	$visitedList->{"$atomParent"} = 1;
	return  (\@newatoms);
}

sub connectivityCheck
{
# this routine will see if the atom listed in $unionref are all connected via bonds.
# This will help catch mistakes in .bif files, where a bond may be omitted.
	my ($unionref,$chid)=@_;
	my %union=%{$unionref};
        my %visitedList; my $visitedString;
	my @nextround;
	my $size =keys %union;
	if($size == 0){
		smog_quit("Found 0 atoms in chain $chid.  Perhaps TER appears on consecutive lines, or TER is immediately followed by END.");
	}

	$nextround[0]=0;
	while( $#nextround >= 0){
		my @newlist;
		foreach my $atomIn(@nextround){
			my $tmp=connectivityHelper(\%union,$atomIn,\%visitedList); ## Traverse through bond graph
		push(@newlist,@{$tmp});
		}
		@nextround=@newlist
    	}

	my $found=0;
	foreach my $atom(keys %visitedList){
		$found++;
	}

	my $missing=0;
	foreach my $atom(keys %union){
		if(!exists $visitedList{$atom}){
			$missing++;
			print "\n!!!Unable to connect the atom at position $atom of chain $chid to the rest of the chain!!!\n";
		}
	}
	return($found,$missing);
}

sub GenerateBondedGeometry {

	my ($connect,$counter,$chid,$chainlength) = @_;
	## $connect is a list of connected residues ##
   	my($bH,$angH,$diheH,$map,$bondMapHashRev,$union,$ConnectedAtoms) = GenAnglesDihedrals($connect,$chainlength);
	my %union=%{$union};

	if($chainlength == 0){
		smog_quit("Found 0 atoms in chain $chid.  Perhaps TER appears on consecutive lines, or TER is immediately followed by END.");
	}elsif($chainlength != 1){
		print "Attempting to connect all atoms in chain $chid to the first atom: ";
		my ($connected,$missed)=connectivityCheck(\%union,$chid);
		if($missed==0 && $connected == $chainlength){
			print "All $connected atoms connected via covalent bonds \n"; 
		}else{
			smog_quit("Only connected $connected of $chainlength atoms in chain $chid.  Unable to connect atoms to the rest of the chain using covalent bond definitions.\nThere may be a missing bond definition in the .bif file.\nSee messages above. ")
		}
	}elsif($chainlength == 1){
		print "Only 1 atom in chain $chid.  Will not perform connectivity checks.";
	}

	# convert and save the connected atoms' numbering, so that we can avoid trouble if we include BONDs later
	my @ConnectedAtoms=@{$ConnectedAtoms};
	my @ConnectedAtoms2;
	foreach my $I(@ConnectedAtoms){
		my $bondStrA = $map->{$I}->[0];
		my $sizeA=$map->{$I}->[2];
		my $ra=$connect->[$map->{$I}->[1]];
		my $ia=$sizeA+getAtomIndexInResidue($ra,$bondStrA);
		push(@ConnectedAtoms2,$ia);
	}

	print "Generating bonds for chain $chid.\n";

    	my @tempArr=();
	## BOND ##
    	for(my $i=0;$i<scalar(@{$bH})-1;$i+=2) {	
		my $bondStrA = $bH->[$i];
		$bondStrA = $map->{$bondStrA}->[0];
		my $bondStrB = $bH->[$i+1];
		$bondStrB = $map->{$bondStrB}->[0];
		my $sizeA = $map->{$bH->[$i]}->[2];my $sizeB = $map->{$bH->[$i+1]}->[2];
		my $ra=$connect->[$map->{$bH->[$i]}->[1]];my $rb=$connect->[$map->{$bH->[$i+1]}->[1]];
		my ($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$bondStrA)
			       ,getAtomBType($ra,$bondStrA));
		my ($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$bondStrB)
			       ,getAtomBType($rb,$bondStrB));
		my $if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");	
		push(@tempArr,pdl($ia,$ib,$if));
	}
	if(@tempArr){
		$BondData{$counter}=cat(@tempArr);
	}
	@tempArr=();
	## ANGLES ##
	print "Generating bond angles for chain $chid.\n";
	foreach my $angs(@{$angH})
	{
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
                my $na;my $nb;my $nc;
		my $ra;my $rb;my $rc;
		my $sizeA; my $sizeB;my $sizeC;
		my($a,$b,$c) = split("-",$angs);
		$na = $map->{$a}->[0];$ra = $connect->[$map->{$a}->[1]];
 		$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
		$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
		$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
		$sizeC=$map->{$c}->[2];

		($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
		($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
		($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));	
        	my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
        	push(@tempArr,pdl($ia,$ib,$ic,$if));		
	}
	if(@tempArr){
		$AngleData{$counter} = cat(@tempArr);
	}
	@tempArr = ();


	## DIHEDRALS ##
	print "Generating dihedral angles for chain $chid.\n";
	foreach my $dihes(@{$diheH})
	{
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
                my $na;my $nb;my $nc;my $nd;
		my $ra;my $rb;my $rc;my $rd;
		my $sizeA; my $sizeB;my $sizeC;my $sizeD;
		my($a,$b,$c,$d) = split("-",$dihes);
		##[AtomName,ResidueIndex,prevSize]##
		$na = $map->{$a}->[0];$ra = $connect->[$map->{$a}->[1]];
 		$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
		$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
		$nd = $map->{$d}->[0];$rd = $connect->[$map->{$d}->[1]];
		$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
		$sizeC=$map->{$c}->[2];$sizeD=$map->{$d}->[2];

		($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
		($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
		($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));
		($id,$td) = ($sizeD+getAtomIndexInResidue($rd,$nd),getAtomBType($rd,$nd));	
		## Adjust args for getEnergyGroup() ##
		if($map->{$b}->[1]-$map->{$c}->[1]>0){
			($nb,$nc)=("$nb?",$nc);	
		}elsif($map->{$b}->[1]-$map->{$c}->[1]<0){
			($nb,$nc)=($nb,"$nc?");	
		}
		my $eG = getEnergyGroup($rb,$rc,$nb,$nc);
		my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);	
        	$eG = $eGRevTable{$eG};
		## [x,y,z,func,countDihedrals,energyGroup]
		push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);

	}

	print "Generating improper angles for chain $chid.\n";

	my ($seconds, $microseconds) = gettimeofday;
	appendImpropers($map,$connect,$bondMapHashRev,\@tempArr,\%union);
  	my ($seconds1, $microseconds1) = gettimeofday;
	$seconds1-=$seconds;
	$microseconds1-=$microseconds;
	$seconds1+=$microseconds/1000000;
	print "dT=$seconds1\n";

	print "Storing dihedral info for chain $chid.\n";
	$DihedralData{$counter} = pdl(@tempArr);
	@tempArr = ();
	return(\@ConnectedAtoms2);
}


sub returnBondTypeFromIndex
{
	my($idx) = @_;
	my $residue = $allAtoms{$idx}->[5];
	my $atom = $allAtoms{$idx}->[3];
	if(!$residue || !$atom)
	{
		smog_quit("Error finding the residue for atom $idx. Perhaps your indices are wrong?");
	}
	if(!$residues{$residue}->{"atoms"}->{$atom})
		{smog_quit("$atom is not part of $residue");}
	return $residues{$residue}->{"atoms"}->{$atom}->{"bType"};
}


sub returnAtomFromIndex
{
	my($idx) = @_;
	return $allAtoms{$idx}->[3];
}
sub returnResidueIndexFromIndex
{
	my($idx) = @_;
	return $allAtoms{$idx}->[2];
}

sub appendImpropersBOND
{

 	my($resA,$resB,$resIDA,$resIDB,$sizeA,$ips,$tempArr) = @_;
 	my ($a,$b,$c,$d) = @{$ips};
 	my $ta=returnBondTypeFromIndex($a);
 	my $tb=returnBondTypeFromIndex($b);
 	my $tc=returnBondTypeFromIndex($c);
 	my $td=returnBondTypeFromIndex($d);
 	my $iia = returnAtomFromIndex($a);
 	my $iib = returnAtomFromIndex($b);
 	my $iic = returnAtomFromIndex($c);
 	my $iid = returnAtomFromIndex($d);
 	if(returnResidueIndexFromIndex($a)!=$resIDA && returnResidueIndexFromIndex($a)!=$resIDB)
 	{smog_quit("Ad-hoc Improper Create Error: Atom $a is part of neither residue $resIDA or $resIDB");}
 	if(returnResidueIndexFromIndex($b)!=$resIDA && returnResidueIndexFromIndex($b)!=$resIDB)
 	{smog_quit("Ad-hoc Improper Create Error: Atom $b is part of neither residue $resIDA or $resIDB");}
 	if(returnResidueIndexFromIndex($c)!=$resIDA && returnResidueIndexFromIndex($c)!=$resIDB)
 	{smog_quit("Ad-hoc Improper Create Error: Atom $c is part of neither residue $resIDA or $resIDB");}
 	if(returnResidueIndexFromIndex($d)!=$resIDA && returnResidueIndexFromIndex($d)!=$resIDB)
 	{smog_quit("Ad-hoc Improper Create Error: Atom $d is part of neither residue $resIDA or $resIDB");}
 	
	my($ia)= (returnResidueIndexFromIndex($a)==$resIDA ? (getAtomIndexInResidue($resA,$iia)) 
	: ($sizeA+getAtomIndexInResidue($resB,$iia)));
 	my($ib)= (returnResidueIndexFromIndex($b)==$resIDA ? (getAtomIndexInResidue($resA,$iib)) 
	: ($sizeA+getAtomIndexInResidue($resB,$iib)));
	my($ic)= (returnResidueIndexFromIndex($c)==$resIDA ? (getAtomIndexInResidue($resA,$iic)) 
	: ($sizeA+getAtomIndexInResidue($resB,$iic)));
	my($id)= (returnResidueIndexFromIndex($d)==$resIDA ? (getAtomIndexInResidue($resA,$iid)) 
	: ($sizeA+getAtomIndexInResidue($resB,$iid)));
 	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
 	## [x,y,z,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
 	push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);
  	
}




sub connCreateInteractionsSingleBOND
{

    	my($consecResiduesH,$sizeA,$counter,$atomA,$atomB,$resAIdx,$resBIdx,$bEG,$imp) = @_;
	my @consecResidues = @{$consecResiduesH};
	my $residue = $consecResidues[1];

    	## AD-HOC BONDS ##
	my($angH,$diheH,$adjList,$bondStrA,$bondStrB)=createConnection($consecResiduesH,0,$atomA,$atomB);
	## BOND ##
	my ($ia,$ta) = (getAtomIndexInResidue($consecResidues[0],$bondStrA)
							,getAtomBType($consecResidues[0],$bondStrA));
	my ($ib,$tb) = ($sizeA+getAtomIndexInResidue($consecResidues[1],$bondStrB)
							,getAtomBType($consecResidues[1],$bondStrB));
	
	
	my $if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");
	
    	$BondData{$counter}=pdl($ia,$ib,$if);
			
	## ANGLES ##
	my @tempArr;
	foreach my $angs(@{$angH})
	{
		my($a,$b,$c) = split("-",$angs);
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
		($ia,$ta) = ($a =~ /(.*)\?/ )
			? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomIndexInResidue($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
				
		($ib,$tb) = ($b =~ /(.*)\?/)
			? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomIndexInResidue($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
				
		($ic,$tc) = ($c =~ /(.*)\?/) 
			? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomIndexInResidue($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
				
			
        	my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
        	push(@tempArr,pdl($ia,$ib,$ic,$if));

	}
        if(@tempArr)
        {
		$AngleData{$counter} = cat(@tempArr);
	}else{
		warn("PDB PARSE WARN:: There are no angles between ",$consecResidues[0]," and ",$consecResidues[1]);
        }
	@tempArr = ();
			
			
	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
		my $eG;
	
		# only save the dihedral if it is centered about the bond we just added
		if(($b =~ m/^$atomA$/ && $c =~ m/^$atomB\?$/) || ($b =~ m/^$atomA\?$/ && $c =~ m/^$atomB$/) 
                    || ($c =~ m/^$atomA$/ && $b =~ m/^$atomB\?$/) || ($c =~ m/^$atomA\?$/ && $b =~ m/^$atomB$/)){
			($ia,$ta) = ($a =~ /(.*)\?/)
				? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
				: (getAtomIndexInResidue($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
					
			($ib,$tb) = ($b =~ /(.*)\?/ )
				? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
				: (getAtomIndexInResidue($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
					
			($ic,$tc) = ($c =~ /(.*)\?/ )
				? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
				: (getAtomIndexInResidue($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
				
			($id,$td) = ($d =~ /(.*)\?/) 
				? ($sizeA+getAtomIndexInResidue($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
				: (getAtomIndexInResidue($consecResidues[0],$d),getAtomBType($consecResidues[0],$d));
		
			if(!$bEG || ($b =~/.*\?/ && $c =~/.*\?/)|| ($b !~/.*\?/ && $c !~/.*\?/))
			{$eG=getEnergyGroup($consecResidues[0],$consecResidues[1],$b,$c);}
			else{$eG=$bEG;}
			my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);
			$eG = $eGRevTable{$eG};
			push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);	
		}		
	}
	     
	   
	## Manually add Improper dihedrals ##
	if(scalar(@{$imp})!=0){
	appendImpropersBOND($consecResidues[0],$consecResidues[1],$resAIdx,$resBIdx,$sizeA,$imp,\@tempArr);
	}
		
        if(@tempArr)
        {$DihedralData{$counter} = pdl(@tempArr);
        }else{
		warn("PDB PARSE WARN:: There are no dihedrals between ",$consecResidues[0]," and ",$consecResidues[1]);
	}
	@tempArr = ();
}

sub appendImpropers
{
my($map,$connect,$bondMapHashRev,$tempArr,$union) = @_;
my %union=%{$union};
my @connImproper; my $connHandle;
my %bondMapHashRev=%{$bondMapHashRev};
#loop through the residues in the chain
	for(my $resIndA=0;$resIndA<=$#$connect;$resIndA++){
		my $resA=$connect->[$resIndA];
		my $resAIp = $residues{"$resA"}->{"impropers"};
		# if not terminal, then also check the next residue
		if($resIndA != $#$connect){
			my $resIndB=$resIndA+1;
			my $resB=$connect->[$resIndB];
			my $resBIp = $residues{"$resB"}->{"impropers"};
	   		## WORK RESIDUE B ##
	   		foreach my $ips(@{$resBIp})
	   		{
	      			if(! (defined $ips) ) {next;}
	  
	  			my $ia;my $ib;my $ic;my $id;
	  			my $ta;my $tb;my $tc;my $td;
	        	        my $na;my $nb;my $nc;my $nd;
	  			my $ra;my $rb;my $rc;my $rd;
	  			my $sizeA; my $sizeB;my $sizeC;my $sizeD;
	  			my($a,$b,$c,$d) = @{$ips};
	  
	  
	   			$a=$bondMapHashRev{"$a-$resIndB"};
	   			$b=$bondMapHashRev{"$b-$resIndB"};
	   			$c=$bondMapHashRev{"$c-$resIndB"};
	   			$d=$bondMapHashRev{"$d-$resIndB"};
	  
	  			my $IMPFLAG1=0;
	  			my $IMPFLAG2=0;
	  			my @TMPARR2 = ($a,$b,$c,$d);
	  			# checking that one of the atoms is bonded to the other three
	  			for(my $I=0;$I<4;$I++){
	  				foreach my $VAL(@{$union->{$TMPARR2[$I]}}){
	  					if($VAL == $TMPARR2[0] || $VAL == $TMPARR2[1] ||$VAL == $TMPARR2[2] ||$VAL == $TMPARR2[3] ){
	  						$IMPFLAG1++;
	  					}
	  				}
	  				if($IMPFLAG1==3){$IMPFLAG2=1;}
	  				$IMPFLAG1=0;
	  			}
	  
	  			##[AtomName,ResidueIndex,prevSize]##
	  			$na = $map->{$a}->[0];
	  			$ra = $connect->[$map->{$a}->[1]];
	   			$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
	  			$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
	  			$nd = $map->{$d}->[0];$rd = $connect->[$map->{$d}->[1]];
	  			$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
	  			$sizeC=$map->{$c}->[2];$sizeD=$map->{$d}->[2];
	  
	  			($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
	  			($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
	  			($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));
	  			($id,$td) = ($sizeD+getAtomIndexInResidue($rd,$nd),getAtomBType($rd,$nd));	
	        	  	if($IMPFLAG2==0){
	  				smog_quit("There is an incorrectly formed improper dihedral. Three atoms must be bonded to a central atom. Improper defined by atoms $ia-$ib-$ic-$id.\nThere may be a missing bond, or incorrectly defined improper in the .bif file.\n");
	  			}
	  
	  		
	  			## Adjust args for getEnergyGroup() ##
	        	  	($nb,$nc) =  ($map->{$b}->[1]-$map->{$c}->[1]==0)?($nb,$nc):("nb?",$nc);
	  			my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");	
	  			push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
	   		}
	  
	  
	   		## WORK ON INTER-RESIDUAL IMPROPERS ##
	   		### CHANGE THIS, ONLY HANDLES SINGLE IMPROPERS ###
	   		$connHandle = $connections{$residues{$resA}->{"residueType"}}->{$residues{$resB}->{"residueType"}};
	   		#@connImproper = @{$connHandle->{"improper"}};
	   		foreach my $ips(@{$connHandle->{"improper"}})
	   		{
	  			if(exists $ips->{"atom"}){ 
	    			my ($a,$b,$c,$d) = @{$ips->{"atom"}}; 
	  			my $ia;my $ib;my $ic;my $id;
	  			my $ta;my $tb;my $tc;my $td;
	        	        my $na;my $nb;my $nc;my $nd;
	  			my $ra;my $rb;my $rc;my $rd;
	  			my $sizeA; my $sizeB;my $sizeC;my $sizeD;
	  			my($an,$bn,$cn,$dn) = @{$ips->{"atom"}};
	  
	  
	        		if($a =~ /[*?^&!@#%()-]/){smog_quit ("Special characters not permitted in connection atoms: $a found.")};
	        		if($b =~ /[*?^&!@#%()-]/){smog_quit ("Special characters not permitted in connection atoms: $b found.")};
	        		if($c =~ /[*?^&!@#%()-]/){smog_quit ("Special characters not permitted in connection atoms: $c found.")};
	        		if($d =~ /[*?^&!@#%()-]/){smog_quit ("Special characters not permitted in connection atoms: $d found.")};
	  
	  
	  			if( $a =~ s/\+$//g ){
	   				$a=$bondMapHashRev{"$a-$resIndB"};
	  			}else{
	   				$a=$bondMapHashRev{"$a-$resIndA"};
	  			}
	  			if( $b =~ s/\+$//g ){
	   				$b=$bondMapHashRev{"$b-$resIndB"};
	  			}else{
	   				$b=$bondMapHashRev{"$b-$resIndA"};
	  			}
	  			if( $c =~ s/\+$//g ){
	   				$c=$bondMapHashRev{"$c-$resIndB"};
	  			}else{
	   				$c=$bondMapHashRev{"$c-$resIndA"};
	  			}
	  			if( $d =~ s/\+$//g ){
	   				$d=$bondMapHashRev{"$d-$resIndB"};
	  			}else{
	   				$d=$bondMapHashRev{"$d-$resIndA"};
	  			}
	  			my $IMPFLAG1=0;
	  			my $IMPFLAG2=0;
	  			my @TMPARR2 = ($a,$b,$c,$d);
	  			for(my $I=0;$I<4;$I++){
	  				foreach my $VAL(@{$union->{$TMPARR2[$I]}}){
	  					if($VAL == $TMPARR2[0] || $VAL == $TMPARR2[1] ||$VAL == $TMPARR2[2] ||$VAL == $TMPARR2[3] ){
	  						$IMPFLAG1++;
	  					}
	  				}
	  				if($IMPFLAG1==3){$IMPFLAG2=1;}
	  				$IMPFLAG1=0;
	  			}
	  
	  			##[AtomName,ResidueIndex,prevSize]##
	  			$na = $map->{$a}->[0];
	  			$ra = $connect->[$map->{$a}->[1]];
	   			$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
	  			$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
	  			$nd = $map->{$d}->[0];$rd = $connect->[$map->{$d}->[1]];
	  			$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
	  			$sizeC=$map->{$c}->[2];$sizeD=$map->{$d}->[2];
	  
	  			($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
	  			($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
	  			($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));
	  			($id,$td) = ($sizeD+getAtomIndexInResidue($rd,$nd),getAtomBType($rd,$nd));	
	  
	        	  	if($IMPFLAG2==0){
	  				smog_quit("There is an incorrectly formed improper dihedral. Three atoms must be bonded to a central atom. Improper defined by atoms $ia-$ib-$ic-$id.\nThere may be a missing bond, or incorrectly defined improper in the .bif file.\n");
	  			}
	  
	          			($nb,$nc) =  ($map->{$b}->[1]-$map->{$c}->[1]==0)?($nb,$nc):("nb?",$nc);
	  				my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");	
	  				push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
	  			}else{
	  				smog_quit("Internal error 1.  Please inform smog team");
	  			}
			}	
	 	}
	
	 	if($resIndA != 0) {next;}
	 
	 	## WORK RESIDUE A ##
	 	foreach my $ips(@{$resAIp})
	 	{
	    		if(! (defined $ips) ) {next;}
	
			my $ia;my $ib;my $ic;my $id;
			my $ta;my $tb;my $tc;my $td;
	                my $na;my $nb;my $nc;my $nd;
			my $ra;my $rb;my $rc;my $rd;
			my $sizeA; my $sizeB;my $sizeC;my $sizeD;
			my($a,$b,$c,$d) = @{$ips};
	
	 		$a=$bondMapHashRev{"$a-$resIndA"};
	 		$b=$bondMapHashRev{"$b-$resIndA"};
	 		$c=$bondMapHashRev{"$c-$resIndA"};
	 		$d=$bondMapHashRev{"$d-$resIndA"};
			my $IMPFLAG1=0;
			my $IMPFLAG2=0;
			my @TMPARR2 = ($a,$b,$c,$d);
			for(my $I=0;$I<4;$I++){
				foreach my $VAL(@{$union->{$TMPARR2[$I]}}){
					if($VAL == $TMPARR2[0] || $VAL == $TMPARR2[1] ||$VAL == $TMPARR2[2] ||$VAL == $TMPARR2[3] ){
						$IMPFLAG1++;
					}
				}
				if($IMPFLAG1==3){$IMPFLAG2=1;}
				$IMPFLAG1=0;
			}
	
			##[AtomName,ResidueIndex,prevSize]##
			$na = $map->{$a}->[0];
			$ra = $connect->[$map->{$a}->[1]];
	 		$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
			$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
			$nd = $map->{$d}->[0];$rd = $connect->[$map->{$d}->[1]];
			$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
			$sizeC=$map->{$c}->[2];$sizeD=$map->{$d}->[2];
	
			($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
			($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
			($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));
			($id,$td) = ($sizeD+getAtomIndexInResidue($rd,$nd),getAtomBType($rd,$nd));	
	         	if($IMPFLAG2==0){
				smog_quit("There is an incorrectly formed improper dihedral. Three atoms must be bonded to a central atom. Improper defined by atoms $ia-$ib-$ic-$id.\nThere may be a missing bond, or incorrectly defined improper in the .bif file.\n");
			}
	       	
			## Adjust args for getEnergyGroup() ##
	        	($nb,$nc) =  ($map->{$b}->[1]-$map->{$c}->[1]==0)?($nb,$nc):("nb?",$nc);
			my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");	
			push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
		}
	}
}

sub connWildcardMatchAngles
{
	my($a,$b,$c) = @_;
	my $angHandle = $interactions->{"angles"};
	my $funct="";
		
	## WILD CARD MATCHING CONDITIONALS ##
	my $matchScore = 0; my $saveScore = 0; my $matchScoreCount=0; my $symmatch=0;
	foreach my $matches(keys %{$angHandle})
	{
		$matchScore = 0;
		my ($aM,$bM,$cM) = split("-",$matches);
		unless(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
			|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
			|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)){
				if($a =~ /\Q$aM\E/) {$matchScore+=2;} else {$matchScore+=1;}
				if($b =~ /\Q$bM\E/) {$matchScore+=2;} else {$matchScore+=1;}
				if($c =~ /\Q$cM\E/) {$matchScore+=2;} else {$matchScore+=1;}
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

	if(!defined $funct|| $funct eq ""){smog_quit("There is no function for bType combination $a-$b-$c. Check .b file");}
	return $funct;
}

sub connWildcardMatchBond
{
	my($typeA,$typeB) = @_;
	my $funct="";


  ## Check if atoms exists in declaration ##

		## WILD CARD MATCHING CONDITIONALS ##

	## If both bond types exists ##
	if( exists $interactions->{"bonds"}->{$typeA}->{$typeB})
	{
		$funct = $interactions->{"bonds"}->{$typeA}->{$typeB};
	}elsif (
		$typeA ne $typeB && (exists $interactions->{"bonds"}->{$typeA}->{"*"} 
                         && exists $interactions->{"bonds"}->{$typeB}->{"*"})){
		smog_quit ("Wildcard conflict in bonds $typeA-$typeB. Both $typeA-\* and $typeB-\* are defined in .b file. Can not unambiguously assign a function...");
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
		smog_quit ("Unable to unambiguously assign bond types to all bonds in a residue\n Offending btypes are $typeA $typeB. Check .b file");
		}
	}

    return $funct;
}

sub connWildcardMatchImpropers
{
    	my($a,$b,$c,$d) = @_;
	my $diheHandle = $interactions->{"impropers"};
	my $funct="";
	## WILD CARD MATCHING CONDITIONALS ##
	my $matchScore = 0; my $saveScore = 0;
	foreach my $matches(keys %{$diheHandle})
	{
		$matchScore = 0;
		my ($aM,$bM,$cM,$dM) = split("-",$matches);
		if(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
			|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
			|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)
			|| ($d !~ /\Q$dM\E/ && $dM !~ /\Q*\E/)){next;}
		if($a =~ /\Q$aM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($b =~ /\Q$bM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($c =~ /\Q$cM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($d =~ /\Q$dM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($matchScore >= $saveScore)
		{$saveScore = $matchScore;$funct = $diheHandle->{$matches};}
	}
	if(!defined $funct || $funct eq ""){smog_quit("There is no function for bType combination $a-$b-$c-$d. Check .b file");}
	return $funct;
}


sub connWildcardMatchDihes
{
    my($a,$b,$c,$d,$eG) = @_;
	my $diheHandle = $interactions->{"dihedrals"}->{"$eG"};
	my $funct="";
	## WILD CARD MATCHING CONDITIONALS ##
	my $matchScore = 0; my $saveScore = 0;;my $matchScoreCount=0; my $symmatch=0; my $Nd=0;
	my $NumOfKeys=keys %{$diheHandle};
	my @keys = keys %{$diheHandle};
	if($NumOfKeys == 1 && $keys[0] eq "*-*-*-*"){
		$funct = $diheHandle->{$keys[0]};

	}else{
		foreach my $matches(keys %{$diheHandle})
		{
			$Nd++;
			$matchScore = 0;my $saveScore = 0;
			# this step can be done once, rather than for each call.
			my ($aM,$bM,$cM,$dM) = split("-",$matches);
			if($matches eq "*-*-*-*"){
				$matchScore=4;
			}else{
	
			if(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
				|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
				|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)
				|| ($d !~ /\Q$dM\E/ && $dM !~ /\Q*\E/)){next;}
	
			if($a =~ /\Q$aM\E/) {$matchScore+=2;} else {$matchScore+=1;}
			if($b =~ /\Q$bM\E/) {$matchScore+=2;} else {$matchScore+=1;}
			if($c =~ /\Q$cM\E/) {$matchScore+=2;} else {$matchScore+=1;}
			if($d =~ /\Q$dM\E/) {$matchScore+=2;} else {$matchScore+=1;}
	
			}	
	
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
				$saveScore = $matchScore;$funct = $diheHandle->{$matches};
			}
		}
	
		if($Nd ==0){
			smog_quit ("energy group\"$eG\" is used in .bif file, or a BOND line, but it is not defined in .sif file.");
		}
		
		my $sym=0;
		if(($a eq $d and $b eq $c) || ($a eq $b and $b eq $c and $c eq $d)){
			$sym=1;
		}
		if(($symmatch ==0 && $sym == 1 && $matchScoreCount > 1)  || ($symmatch ==0 && $sym == 0 && $matchScoreCount > 0) || ($symmatch ==1 && $sym == 0 && $matchScoreCount > 0) || ($symmatch ==1 && $sym == 1 && $matchScoreCount > 0)){
	
			smog_quit ("$symmatch  $sym $matchScoreCount Multiple possible angles match $a-$b-$c-$d, and energyGroup $eG equally well. Can not determine function based on .b file.");
		}
	
		if($matchScore == 0){
			smog_quit ("Dihedral Angle between bTypes $a-$b-$c-$d and energyGroup $eG: Unable to match to a function in .b file.");
		}
	
		if(!defined $funct || $funct eq ""){smog_quit("There is no function for bType combination $a-$b-$c-$d with energyGroup=$eG. Check .b file");}
	}

	return $funct;
}



sub GenAnglesDihedrals
{
	my($connect,$chainlength) = @_;
	## $connect is a list of connected residues ##
	my @tempA;my @tempD;
	my $bonds;my $dihes; my $angles; 
	my $oneFour;
	my %union; my $connHandle;
	my $atomA = ""; my $atomB = "";
	my $i=0;my $j=0;my $mapCounter=0;
	my %bondMapHash; ##[AtomName,ResidueIndex,prevSize]##
	my %bondMapHashRev;
	my @connectList;
	my @AtomsInConnections;
	my $leftAtom;my $rightAtom;
	my $leftResidue;my $rightResidue;
	my $prevSize = 0;
	
	## Go through list of connected residues ##
	for($i = 0;$i<=$#$connect;$i++)
	{ 
		$j = $i+1;  
  		my $resABonds = $dihedralAdjList{$connect->[$i]};
		my $resAAtoms = $residues{$connect->[$i]}->{"atoms"};
		## Atoms to mapCounter renaming ##
 		foreach my $atom(keys %{$resAAtoms})
		{
			$bondMapHashRev{"$atom-$i"}=$mapCounter;
			$bondMapHash{$mapCounter}=[$atom,$i,$prevSize];
			$mapCounter++;  
		}
		foreach my $atom(keys %{$resABonds})
    		{
     			my @tempArr = map {$bondMapHashRev{"$_-$i"}} @{$resABonds->{$atom}};
			my $atomKey = $bondMapHashRev{"$atom-$i"}; 
     			$union{$atomKey} = \@tempArr;
		}

		# if this is a single residue chain, then don't try to connect it to the next residue
        	if($#$connect == 0) {last;}
		## Start of chain no inter residue connection ##
    		#  but setup leftAtom and leftResidue sizes #
    		if($i == 0) 
		{
			if(!exists $connections{$residues{$connect->[0]}->{"residueType"}}->{$residues{$connect->[1]}->{"residueType"}}){
			my $typeA=$residues{$connect->[0]}->{"residueType"};
			my $typeB=$residues{$connect->[1]}->{"residueType"};
			smog_quit("Connection not defined between residues of type $typeA and $typeB. Check .bif file.  Issue encountered when connecting residue 0 and 1 (residue index within chain, starting at 0).")
		}
		$connHandle = $connections{$residues{$connect->[0]}->{"residueType"}}->{$residues{$connect->[1]}->{"residueType"}};
	  	$leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
	  	$leftAtom = $bondMapHashRev{"$leftAtom-$i"};
	 	$prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});
      		 next;
		}
        	## $i > 0, create inter residue connection ##
		## $i-1 <--> $i
		if(!exists $connections{$residues{$connect->[$i-1]}->{"residueType"}}->{$residues{$connect->[$i]}->{"residueType"}}){
			my $typeA=$residues{$connect->[$i-1]}->{"residueType"};
			my $typeB=$residues{$connect->[$i]}->{"residueType"};
			my $ii=$i-1;
			smog_quit("Connection not defined between residues of type $typeA and $typeB. Check .bif file. Issue encountered when connecting residue $ii and $i (residue index within chain, starting at 0)")
		}
		$connHandle = $connections{$residues{$connect->[$i-1]}->{"residueType"}}->{$residues{$connect->[$i]}->{"residueType"}};
		$rightAtom = $connHandle->{"bond"}->[0]->{"atom"}->[1];
	    	$rightAtom = $bondMapHashRev{"$rightAtom-$i"};
		push(@AtomsInConnections,$leftAtom);
		push(@AtomsInConnections,$rightAtom);
		push(@connectList,$leftAtom);
		push(@connectList,$rightAtom);
    
    		## $i <--> $i+1
        	if($i == $#$connect) {last;}
		if(!exists $connections{$residues{$connect->[$i]}->{"residueType"}}->{$residues{$connect->[$i+1]}->{"residueType"}}){
			my $typeA=$residues{$connect->[$i]}->{"residueType"};
			my $typeB=$residues{$connect->[$i+1]}->{"residueType"};
			my $ii=$i+1;
			smog_quit("Connection not defined between residues of type $typeA and $typeB. Check .bif file. Issue encountered when connecting residue $i and $ii (residue index within chain, starting at 0)")
		}
        	$connHandle = $connections{$residues{$connect->[$i]}->{"residueType"}}->{$residues{$connect->[$i+1]}->{"residueType"}};
    		$leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
    		$leftAtom = $bondMapHashRev{"$leftAtom-$i"};
    		$prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});
   	}
  	## Create Inter residue connection ##
  	for($i=0;$i<scalar(@connectList)-1;$i+=2) {
   		push(@{$union{$connectList[$i]}},$connectList[$i+1]);
   		push(@{$union{$connectList[$i+1]}},$connectList[$i]);
  	}
  	my @includedatoms;
  	for(my $I=0;$I<$chainlength;$I++){
   		$includedatoms[$I]=0;
  	}
  	foreach my $atom(keys %union){
   		$includedatoms[$atom]=1;
	}
	for(my $I=0;$I<$chainlength;$I++){
		if($includedatoms[$I]==0 && $chainlength !=1){
			my $atomname=$bondMapHash{$I}->[0];
			my $residue=$bondMapHash{$I}->[1];
			$residue++;
        		smog_quit("No bonds found with atom $atomname in residue $residue. Check .bif definitions.");
   		}
	}

  	($dihes,$angles,$oneFour)=adjListTraversal(\%union);
  	return (\@connectList,$angles,$dihes,\%bondMapHash,\%bondMapHashRev,\%union,\@AtomsInConnections);
}



sub createConnection
{
	my($connect,$firstFlag,$atomA,$atomB,$counter) = @_;
	my @tempA; my @tempD;
	my $dihes; my $angles; my $oneFour;
	my %union; my $connHandle;
	my $resABonds;my $resBBonds;
	my %tempAdjList; 


  	## USES GLOBAL FLAG MISSING TO CREATE NEW MAP ##
   	## Connection via connections attribute ##
   	if(!$atomA || !$atomB){
    		$connHandle = $connections{$residues{$connect->[0]}->{"residueType"}}->{$residues{$connect->[1]}->{"residueType"}};
    		$atomA = $connHandle->{"bond"}->[0]->{"atom"}->[0];
    		$atomB = $connHandle->{"bond"}->[0]->{"atom"}->[1];
    	}
	$resABonds = $dihedralAdjList{$connect->[0]};
	$resBBonds = $dihedralAdjList{$connect->[1]};
	

	## Rename atoms in $resBBonds to avoid clashing ##
 	foreach my $atom(keys %{$resBBonds})
	{
		my @tempArr = map {"$_?"} @{$resBBonds->{$atom}};
		$tempAdjList{"$atom?"} = \@tempArr;
	}
	%union = ();
	while(my($k,$v) = each %{$resABonds}){@{$union{$k}}=@{$v};}
	while(my($k,$v) = each %tempAdjList){$union{$k}=$v;}
	
	## Connect C-N
	push(@{$union{"$atomA"}},"$atomB?");
	push(@{$union{"$atomB?"}},"$atomA");
		
	 
	($dihes,$angles,$oneFour)=adjListTraversal(\%union); 
  	 

	## REMOVE ANY ANGLES/DIHES NOT CONTAIN C-N? or N?-C and any *.-*. without '?' if firstFlag==0
	if($firstFlag==0){
		@tempA = map {$_ =~/\?/ ? ($_) : ()} @{$angles};
		@tempD = map {$_ =~ /\?/ ? ($_) : () } @{$dihes};
		@{$angles} = @tempA; @{$dihes} = @tempD;
	}
	## ADHOC BONDS ##
	## REMOVE ALL ANGLES/DIHES EXCEPT INTERDOMAIN ##
	elsif($firstFlag==-1)
	{
		@tempA = map {countQM($_) > 0 && countQM($_) < 3? ($_) : ()} @{$angles};
		@tempD = map {countQM($_) > 0 && countQM($_) < 4? ($_) : ()} @{$dihes};
		@{$angles} = @tempA; @{$dihes} = @tempD;
	}
	return ($angles,$dihes,\%union,"$atomA","$atomB");
}

sub countQM{my $in=shift;my $c = () = $in =~ /\?/g;return $c;}


sub catPDL
{
	foreach my $pdls(keys %tempPDL)
   	{
     		my @arrConvert = values %{$tempPDL{$pdls}};
	 	$resPDL{$pdls} = cat(@arrConvert);
	 	delete $tempPDL{$pdls};
   	}
}


sub parseCONTACT
{
	#lets leave this as two filename inputs in case we want to allow two sources of contacts in the future (i.e. user and shadow)
	my($fileName,$fileName2,$userProvidedMap,$CGenabled) = @_;
	my $numContacts = 0; my $garbage = 0;
	my $line = "";
	my $chain1;my $chain2; my $contact1; my $contact2; my $pdbNum1; my $pdbNum2; my $res1; my $res2;
	my $dist;
	my $x1;my $x2;my $y1;my $y2;my $z1;my $z2;
	my %resContactHash; my $skip = 0; my $COARSECONT;
	my @interiorTempPDL; #usage: push(@interiorTempPDL,[1,$contact1,$contact2,$dist]);
	#Format for this PDL has a boolean as the first argument
	#it is unused for now, but could be useful in future to use
	#as a flag to differentiate between user generated and smog generated contacts

	if(!$userProvidedMap){ #use shadow generated contact map
		## OPEN .contact FILE ##
		unless (open(CONTFILE, $fileName)) {
			smog_quit ("Internal contact file can not be read.  See shadow.log for more information.");
		}
		my $coarseFile = $fileName.".CG";
		if($CGenabled == 1) { unless (open($COARSECONT,">$coarseFile")) {
			smog_quit ("Internal contact file cannot be written.");
		} }
		while($line = <CONTFILE>)
		{
			($chain1,$contact1,$chain2,$contact2) = split(/\s+/,$line);
			if($CGenabled == 1) {
				#if we are coarse graining then we need to map AA to residue contacts
				#first see which residues they belong
				$res1 = $allAtomsBackup{$contact1}[9]; $res2 = $allAtomsBackup{$contact2}[9];
				$skip = 1; 	#skip if this res1-res2 has already been added
				if(!exists $resContactHash{"$res1,$res2"}) {
					$resContactHash{"$res1,$res2"} = 1;
					$skip = 0;
					$contact1 = $res1; $contact2 = $res2;
					print $COARSECONT "$chain1 $res1 $chain2 $res2\n";
				}
			}
			if($skip == 0) { #maybe we skip sometimes if coarse graining				
				$x1 = $allAtoms{$contact1}[6];$y1 = $allAtoms{$contact1}[7];$z1 = $allAtoms{$contact1}[8];
				$x2 = $allAtoms{$contact2}[6];$y2 = $allAtoms{$contact2}[7];$z2 = $allAtoms{$contact2}[8];
				$dist = sqrt( ($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2) * $angToNano;
				if($dist < $interactionThreshold->{"contacts"}->{"shortContacts"})
				{
				  if($main::setContacttoLimit){
				    $dist=$interactionThreshold->{"contacts"}->{"shortContacts"};
	                            print "NOTE: Contact between atoms $contact1 $contact2 exceed contacts threshold with value $dist\n";
				    print "-limitcontactlength is being used, will set distance of contact to $dist\n\n";
				  }else{
	                            smog_quit("Contact between atoms $contact1 $contact2 exceed contacts threshold with value $dist");
	 			  }
				}
				push(@interiorTempPDL,[0,$contact1,$contact2,$dist]);
				$numContacts++;
			}
		}
		if($CGenabled == 1) { close($COARSECONT); }
	} else { #read in contact from file
		print "\nNOTE: Not calculating contact map\n";
		print "Reading contacts from $fileName2\n";
		## OPEN user provided contact FILE ##
		unless (open(CONTFILE1, $fileName2)) {
			smog_quit ("Cannot read contact file '$fileName2'.");
		}
		## User contact map should be in format below and use input PDB numbering ##
		## chain1 atom1 chain2 atom2 ##
		while($line = <CONTFILE1>) {
			my ($chain1,$pdbNum1,$chain2,$pdbNum2,$dist) = split(/\s+/,$line);
			$chain1--;$chain2--; #moving to zero based numbering
			if(!exists $indexMap{"$chain1-$pdbNum1"}) { 
				$chain1++;
				smog_quit("Seems that PDB number $pdbNum1 in chain $chain1 does not exist. Check input contact map.\n");
			}
			if(!exists $indexMap{"$chain2-$pdbNum2"}) { 
				$chain2++;
				smog_quit("Seems that PDB number $pdbNum2 in chain $chain2 does not exist. Check input contact map.\n");
			}
			$contact1 = $indexMap{"$chain1-$pdbNum1"};
			$contact2 = $indexMap{"$chain2-$pdbNum2"};
			if($dist) { #check if distance provided
				if(whatAmI($dist)!=3) { #check if it is numeric
					if($dist < 0 || $dist > 1000) { #check that it is a sensible number
						$chain1++;$chain2++;
						smog_quit("Input contact map has distance for contact $chain1 $pdbNum1 $chain2 $pdbNum2 less than 0 or greater than 1000nm. Maybe something is wrong? Distance is: $dist.\n");
					} 
				} else {
					$chain1++;$chain2++;
					smog_quit("Input contact map has non-numeric distance for contact $chain1 $pdbNum1 $chain2 $pdbNum2: $dist.\n");
				}
			} else { #distace was not provided, lets calculate it from structure
				$x1 = $allAtoms{$contact1}[6];$y1 = $allAtoms{$contact1}[7];$z1 = $allAtoms{$contact1}[8];
				$x2 = $allAtoms{$contact2}[6];$y2 = $allAtoms{$contact2}[7];$z2 = $allAtoms{$contact2}[8];
				$dist = sqrt( ($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2);
			}
			$dist = $dist * $angToNano;
			if(!exists $allAtoms{$contact1}){warn("ATOM $contact1 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}
			if(!exists $allAtoms{$contact2}){warn("ATOM $contact2 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}
			if($dist < $interactionThreshold->{"contacts"}->{"shortContacts"})
			{

			  if($main::setContacttoLimit){
			    $dist=$interactionThreshold->{"contacts"}->{"shortContacts"};
                            print "CONTACT between atoms $contact1 $contact2 exceed contacts threshold with value $dist\n";
			    print "-limitcontactlength is being used, will set distance of contact to $dist\n\n";
			  }else{
                            smog_quit("CONTACT between atoms $contact1 $contact2 exceed contacts threshold with value $dist");
 			  }
		        }
			push(@interiorTempPDL,[1,$contact1,$contact2,$dist]);
			$numContacts++;
		}
	} 
	$contactPDL = pdl(@interiorTempPDL);
	return $numContacts;
  
}

sub whatAmI {
	if($_[0] =~ /^[0-9,eE]+$/) {return 1;} #integer
	if($_[0] =~ /^[0-9,.Ee+-]+$/) {return 2;} #float
	return 3; #not integer or float
}

1;
