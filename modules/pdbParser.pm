#!/usr/bin/perl -w
##############################################################################
# pdbParser.pl: parses PDB file and obtains ATOM, residue and coordinate info.
# PDB file has to comply to the standard column format for each attributes.
# Author: Mohit Raghunathan													
# Date: May 2012															 
##############################################################################
package pdbParser;

use bifParser;
use setRatios;
use strict;
use warnings;
use Data::Dumper;
use Exporter;
use PDL; ## LOAD PDL MODULE
use Carp;

## DELETE LATER ##
#use Devel::Size qw(size total_size);

## DECLEARATION TO SHAR DATA STRUCTURES ##
our @ISA = 'Exporter';
our @EXPORT = 
qw($interactionThreshold %residues $termRatios %allAtoms parseCONTACT $contactPDL parseATOM catPDL $totalAtoms returnFunction intToFunc funcToInt %connAngleFunctionals %connDiheFunctionals %connBondFunctionals %resPDL %connPDL %bondFunctionals %dihedralFunctionals %angleFunctionals setInputFileName parseBif parseSif parseBonds createBondFunctionals createDihedralAngleFunctionals parseNonBonds getContactFunctionals $contactSettings $interactions clearPDBMemory clearBifMemory parseATOMCoarse);

my @vector;
my $coorPDL;
my %results;
my %residueCount = ("ASN"=>8);
my %residuePDL=();

my %testASN = ("N"=>0,"CA"=>1,"C"=>2,"O"=>3,"CB"=>4,"CG"=>5,"OD1"=>6,"ND2"=>7); 

our %tempPDL = ();
our %resPDL;
our %connPDL;


our %connAngleFunctionals;
our %connDiheFunctionals;
our %connBondFunctionals;

my @consecResidues;


our $totalAtoms;
our $contactPDL;
our %allAtoms;

our %extContacts;


###########################
## CLEAR VARIABLE MEMORY ##
###########################
sub clearPDBMemory {
undef %tempPDL;undef %resPDL;undef %connPDL;
undef %connAngleFunctionals;undef %connDiheFunctionals;
undef %connBondFunctionals;undef $totalAtoms;
undef $contactPDL; undef %allAtoms;
}


####################################################################
# parseATOM(\@pdbLines,\%hashRef):extract x, y, and z coordinates, #
# serial number and element symbol from PDB ATOM record type. 	   #
####################################################################
sub parseATOM
{
  my ($fileName) = @_;
  ## INTERNAL VARIABLES ##
  my $counter = 0;
  my @temp; my @connResA; my @connResB; my @union;
  my $x;my $y;my $z;
  my $residue; my $interiorResidue; my $atom;my $atomSerial;
  my $resCount; my $lineEnd;
  my $i; my $putIndex=0; my $strLength;
  my $resType;
  my $angH; my $diheH;
  my $bondStrA;my $bondStrB;my $typeA;my $typeB;
  my $endFlag=0; my $headFlag=1;my $outLength;
  $totalAtoms = 0;my $nbType;my $residueType;
  my $atomCounter=0;my $singleFlag = 1;
  my $chainNumber = 0;my $linkFlag = 0;
  my $residueIndex=0;
  my $pdbResidueIndex=0; my $interiorPdbResidueIndex=0;
  my %indexMap;my $lineNumber = 0;
  my $PDB; 
 
  ## OPEN .PDB FILE ##
 unless (open($PDB, $fileName)) {
    confess "Cannot read from '$fileName'.\nProgram closing.\n";
}

  ## LOOP THROUGH EACH LINE ##
 while(my $record = <$PDB>)
 {
 
 my @impAtoms = ();
 $lineNumber++;
 ## IF END LINE ##
 if($record =~ m/^END/) 
 {last;}
 ## PARSE BOND LINES ##
 if($record =~ m/^BOND/)
 {
    chomp($record);
   
    my($trig,$atoma,$atomb,$eG) = split(/\s+/,$record);
    my $idxA = $indexMap{$atoma};
    my $idxB = $indexMap{$atomb};
    my $resA = $allAtoms{$idxA}->[5];
    my $resB = $allAtoms{$idxB}->[5];
    my $atomA = $allAtoms{$idxA}->[3];
    my $atomB = $allAtoms{$idxB}->[3];
    my $resAIdx = $allAtoms{$idxA}->[2];
    my $resBIdx = $allAtoms{$idxB}->[2];

    my $sizeA = scalar(keys %{$residues{$resA}->{"atoms"}});
    $counter++;
    my $union;
    $union=($tempPDL{$resA}->{$resAIdx})->glue(1,$tempPDL{$resB}->{$resBIdx});
    print "\nNOTE:";
    print "Generating ad-hoc BOND between atoms $atoma,$atomb\n";
    $connPDL{$counter}=$union;
    ## Check if improper directive is present ##
    if($record =~ m/IMPROPER/)
    {
     my($left,$right) = split(/IMPROPER/,$record);
     $right =~ s/^\s+|\s+$//g;
     @impAtoms = split(/\s+/,$right);
     print "IMPROPER DETECTED @impAtoms\n";
    }
    
    
    connCreateInteractionsBOND([$resA,$resB],$sizeA,$counter,$atomA,$atomB,$resAIdx,$resBIdx,$eG,\@impAtoms); 
 }

	## IF TER LINE  ##
	if($record =~ m/^TER/)
	{
	   
			$chainNumber++; ## INCREMENT CHAIN NUMBER ##
			## PREV CHAIN WAS SINGLE ## 
			if($headFlag == 1)
			{
					singleCreateInteractions($residue,$counter);
					$connPDL{$counter}=pdl(@connResA);
			}
			## ASSUME START OF PDB, headFlag = 1 ##
			$headFlag = 1;next;
	} 
	
	## ONLY WORK WITH ATOM LINES ##
	if($record =~ m/^ATOM/ || $record =~ m/^HETATM/)
	{
	
		$outLength = length($record);
	    ## OBTAIN RESIDUE NAME ##
		$residue = substr($record,17,4);
		$residue =~ s/^\s+|\s+$//g;
        $pdbResidueIndex = substr($record,22,4); 
        if(!exists $residues{$residue}){confess "\nPDB PARSE ERROR: \"$residue\" doesn't exist in .bif at line $lineNumber\n";}


		$resCount = scalar(keys(%{$residues{$residue}->{"atoms"}}));
	    seek($PDB, -$outLength, 1); # place the same line back onto the filehandle
	
	    ## Incr local residue index ## 
	    $residueIndex++;
        $lineNumber--;
	    ## Parse residue hash ##
	     my %uniqueAtom;
		for($i=0;$i<$resCount;$i++)
		{
            $lineNumber++;
			$record = <$PDB>;
			if($record !~ m/^ATOM|^HETATM/)
			{confess("PDB PARSE ERROR\n Expected ATOM or HETATM line at Line $lineNumber. Residue $residue might have been truncated at $lineNumber\n");}
			$interiorResidue = substr($record,17,4);
			$interiorResidue =~ s/^\s+|\s+$//g;
            $interiorPdbResidueIndex = substr($record,22,4);  
			## CHECK IF ALL ATOMS CONFORM TO BIF RESIDUE DECLARATION ##
			if($interiorResidue !~ /$residue/
            || $pdbResidueIndex != $interiorPdbResidueIndex)
			{
			$residueIndex--;$lineNumber--;
            confess "\nPDB PARSE ERROR\nThere must be missing atoms in residue $residue at line $lineNumber. $residue at $lineNumber must have been truncated, because the next entry is $interiorResidue\n";
			}
			
			$x = substr($record, 30, 8);
			$y = substr($record, 38, 8);
			$z = substr($record, 46, 8);
			$atom = substr($record, 12, 4);	
			#$atomSerial = substr($record,6,5);
			$atomCounter++;
			$atomSerial=$atomCounter;
			$atomSerial =~ s/^\s+|\s+$//g;	
			$atom =~ s/^\s+|\s+$//g;

                       if(exists $uniqueAtom{$atom})
			{confess("\nPDB PARSE ERROR: $atom is a duplicate in $residue at line $lineNumber\n");} 	
                	else {$uniqueAtom{$atom}=1;}
			if(!exists $residues{$residue}){confess "\nPDB PARSE ERROR: \"$residue\" doesn't exist in .bif at line $lineNumber";}
			if(!exists $residues{$residue}->{"atoms"}->{$atom}){	
                confess "\nPDB PARSE ERROR:$atom doesn't exists in .bif declaration of $residue at line $lineNumber.\n\n";}			
			
			$putIndex = $residues{$residue}->{"atoms"}->{$atom}->{"index"};
			$nbType = $residues{$residue}->{"atoms"}->{$atom}->{"nbType"};
			$residueType = $residues{$residue}->{"type"};
			##[nbType,residueType,resIndex,atom,chainNumber,resName]
			$allAtoms{$atomSerial}=[$nbType,$residueType,$residueIndex,$atom,$chainNumber,$residue];
			my $pdbIndex = substr($record,6,5);
			$pdbIndex =~ s/^\s+|\s+$//g;
			$indexMap{$atomSerial}=$pdbIndex;

			$temp[$putIndex]=[$x,$y,$z,$atomSerial];
			if($residues{$residue}->{"atomCount"} == -1){$totalAtoms++;}
		}
		
		## Add total atoms in residue ##
		if($residues{$residue}->{"atomCount"} != -1)
		{
			$totalAtoms+=($residues{$residue}->{"atomCount"});
		}
		 

		if($i != $resCount){confess "\nPDB PARSE ERROR\nTotal number of atoms of $residue at line $lineNumber doesn't match with .bif declaration\n\n";}
		
		$counter++;
        ## HEAD RESIDUE WAIT FOR NEXT RESIDUE TO CONNECT ##
		if($headFlag) 
		{
		   @connResA = @temp;@connResB = @temp;
		   $consecResidues[0] = $residue;
		   $consecResidues[1] = $residue;
		   $headFlag = 0;
		   $linkFlag = 1;
		}		
		## CONNECT TO N AND N-1 RESIDUE ##
		elsif(!$headFlag)
		{
		   my @testArr; my $sizeA; my $sizeB;
		   @connResA = @connResB;
	       @connResB=@temp;
		   $consecResidues[0] = $consecResidues[1];
	       $consecResidues[1] = $residue;   
           $sizeA = scalar(@connResA);
		   @union = (@connResA,@connResB);
		   ## Extract Index info ##
	       connCreateInteractions(\@consecResidues,$sizeA,$counter,$linkFlag,"","","");
		   $headFlag = 0;
		   $connPDL{$counter}=pdl(@union);
		   @union = ();
		   @connResA = @connResB;
		   $headFlag = 0;$linkFlag = 0;
		}
		#else {confess "PDB PARSE ERROR";}
		$tempPDL{$residue}->{$residueIndex}=pdl(@temp);
		@temp = ();
	}
	$record = "";
	
 }
			## ONLY SINGLE RESIDUE IN PDB ## 
			if($headFlag == 1)
			{
					singleCreateInteractions($residue,$counter);
					$connPDL{$counter}=pdl(@connResA);
			}
	close($PDB);
}



#########################################################################
# parseExternalContacts($filename)
####################################################################
sub parseExternalContacts
{
  my($file) = @_;

  


}


####################################################################
# parseATOMCoarse(\@pdbLines,\%hashRef):extract x, y, and z coordinates, #
# serial number and element symbol from PDB ATOM record type. 	   #
####################################################################
sub parseATOMCoarse
{

  my ($fileName) = @_;

  ## INTERNAL VARIABLES ##
  my $counter = 0;
  my @temp; my @connResA; my @connResB; my @union;
  my @tempBond;
  my @consecResidues;
  my $x;my $y;my $z;
  my $residue; my $interiorResidue; my $atom;my $atomSerial;
  my $resCount; my $lineEnd;
  my $i; my $putIndex=0; my $strLength;
  my $resType;
  my $angH; my $diheH;
  my $bondStrA;my $bondStrB;my $typeA;my $typeB;
  my $endFlag=0; my $headFlag=1;my $outLength;
  $totalAtoms = 0;my $nbType;my $residueType;
  my $atomCounter=0;my $singleFlag = 1;
  my $chainNumber = 0;my $linkFlag = 0;
  my $residueIndex=1;
   
  ## OPEN .PDB FILE ##
 unless (open(MYFILE, $fileName)) {
    confess "Cannot read from '$fileName'.\nProgram closing.\n";
}

  ## LOOP THROUGH EACH LINE ##
 while(my $record = <MYFILE>)
 {
	## IF TER LINE  ##
	if($record =~ m/TER|END/)
	{
			$chainNumber++; ## INCREMENT CHAIN NUMBER ##
			## PREV CHAIN WAS SINGLE ## 
			if($headFlag == 1)
			{
				confess("There is a chain with only one residue. Cannot coarse grain\n");
			}
			## CREATE INTERACTION ##
            coarseCreateInteractions(\@consecResidues,$counter);
		   	$connPDL{$counter}=pdl(@union);
			@union = ();$counter++;
            @consecResidues = ();
			if($record =~ m/END/){last;}
			else {next;}
	} 
	
	## ONLY WORK WITH ATOM LINES ##
	if($record =~ m/ATOM/ || $record =~ m/HETATM/)
	{
	
		$outLength = length($record);
	 	## OBTAIN RESIDUE NAME ##
		$residue = substr($record,17,4);
		
		$residue =~ s/^\s+|\s+$//g;
		$resCount = scalar(keys(%{$residueBackup{$residue}->{"atoms"}}));
	 	seek(MYFILE, -$outLength, 1); # place the same line back onto the filehandle
	
		for($i=0;$i<$resCount;$i++)
		{
			$record = <MYFILE>;

			$interiorResidue = substr($record,17,4);
			$interiorResidue =~ s/^\s+|\s+$//g;
			## CHECK IF ALL ATOMS CONFORM TO BIF RESIDUE DECLARATION ##
			if($interiorResidue !~ /$residue/)
			{confess "PDB PARSE ERROR\nResidue doesn't conform with coarse grain .bif:: $record";}
			$atom = substr($record, 12, 4);
			$atom =~ s/^\s+|\s+$//g;
			if(!exists $residueBackup{$residue}->{"atoms"}->{$atom})
			{confess "PDB PARSE ERROR\n$atom doesn't exists in .bif declaration of $residue\n\n";}
			
			## CHECK IF ATOM IS COARSE GRAINED ##
                        if(!exists $residues{$residue}->{"atoms"}->{$atom}){next;}
			
			$x = substr($record, 30, 8);
			$y = substr($record, 38, 8);
			$z = substr($record, 46, 8);
			$atomCounter++;
			$atomSerial=$atomCounter;
			$atomSerial =~ s/^\s+|\s+$//g;	
			
			$putIndex = $residues{$residue}->{"atoms"}->{$atom}->{"index"};
			$nbType = $residues{$residue}->{"atoms"}->{$atom}->{"nbType"};
			$residueType = $residues{$residue}->{"type"};
		    ##[nbType,residueType,resIndex,atom,chainNumber,resName]
			$allAtoms{$atomSerial}=[$nbType,$residueType,$residueIndex,$atom,$chainNumber,$residue]; ## SAVE UNIQUE NBTYPES --> obtain info from nbtype

			$temp[$putIndex]=[$x,$y,$z,$atomSerial];
            $tempBond[$putIndex]=[$x,$y,$z,$atomSerial];
			$totalAtoms++;
		}

		if($i != $resCount){confess "PDB PARSE ERROR\nTotal number of atoms of $residue doesn't match with .bif declaration\n\n";}
		## CONCAT RESIDUE ##
	  	@union = (@union,@temp);@temp=();
		push(@consecResidues,$residue);
		$headFlag = 0;
        $tempPDL{$residue}->{$residueIndex}=pdl(@tempBond);
		@tempBond = ();
		$residueIndex++;
				
	}
	$record = "";
	
 }


			## ONLY SINGLE RESIDUE IN PDB NO COARSE GRAINING ## 
			if($headFlag == 1)
			{
				confess("There is only a single residue in your system, cannot coarse grain.\n");
			}
}



##
# returnFunction: Return the fType of a specified function
sub returnFunction
{
 my($funcString) = @_;
 if(!exists $functions->{$funcString}){confess "$funcString cannot be found";}
 return $functions->{$funcString}->{"fType"};

}

##
# getAtomAbsoluteIndex: Return the index of an atom from storef indexing
sub getAtomAbsoluteIndex
{
 my($residue,$atom) = @_;
 if(!exists $residues{$residue}){confess "$residue wasn't defined in bif";}
 if(!exists $residues{$residue}->{"atoms"}->{$atom}){confess "$atom wasn't defined in $residue in the bif";}
 return $residues{$residue}->{"atoms"}->{$atom}->{"index"};

}

##
# getAtomBType: Return the bondType of an atom
sub getAtomBType
{
 my($residue,$atom) = @_;
 return $residues{$residue}->{"atoms"}->{$atom}->{"bType"};
}


##
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
		elsif(exists $residues{$residueIn}->{"rigidGroups"}->{"$atoma-$atomb"})
			{return $residues{$residueIn}->{"rigidGroups"}->{"$atoma-$atomb"};}
		else{confess("A specified energy group for $residuea:$atoma, $residueb:$atomb doesn't exists");}
	}
 	## If Bond is between two residues ##
	else
	{
		$residueTypea =$residues{$residuea}->{"type"};
		$residueTypeb =$residues{$residueb}->{"type"};
		return $connections{$residueTypea}->{$residueTypeb}->{"bond"}->[0]->{"energyGroup"};
	}

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
	 	($ia,$ta) = (getAtomAbsoluteIndex($residue,$a),getAtomBType($residue,$a));
		($ib,$tb) = (getAtomAbsoluteIndex($residue,$b),getAtomBType($residue,$b));
		($ic,$tc) = (getAtomAbsoluteIndex($residue,$c),getAtomBType($residue,$c));
   		my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
   		push(@tempArr,pdl($ia,$ib,$ic,$if));		
	}
		$connAngleFunctionals{$counter} = cat(@tempArr);
		@tempArr = ();
				
	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
 		($ia,$ta) = (getAtomAbsoluteIndex($residue,$a),getAtomBType($residue,$a));
		($ib,$tb) = (getAtomAbsoluteIndex($residue,$b),getAtomBType($residue,$b));
		($ic,$tc) = (getAtomAbsoluteIndex($residue,$c),getAtomBType($residue,$c));
		($id,$td) = (getAtomAbsoluteIndex($residue,$d),getAtomBType($residue,$d));
				
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
		($ia,$ta) = (getAtomAbsoluteIndex($residue,$ia),getAtomBType($residue,$ia));
	 ($ib,$tb) = (getAtomAbsoluteIndex($residue,$ib),getAtomBType($residue,$ib));
		($ic,$tc) = (getAtomAbsoluteIndex($residue,$ic),getAtomBType($residue,$ic));
		($id,$td) = (getAtomAbsoluteIndex($residue,$id),getAtomBType($residue,$id));
	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
	## [a,b,c,d,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
	push(@tempArr,[$ia,$ib,$ic,$id,$if,1,-1]);	
  }
		
		$connDiheFunctionals{$counter} = pdl(@tempArr);
		@tempArr = ();
		
}

sub coarseCreateInteractions {

	my ($connect,$counter) = @_;
	## $connect is a list of connected residues ##
   	my($bH,$angH,$diheH,$map) = createMultiConnections($connect);
    my @tempArr=();
	## BOND ##
    for(my $i=0;$i<scalar(@{$bH})-1;$i+=2) {	
	my $bondStrA = $bH->[$i];
    $bondStrA = $map->{$bondStrA}->[0];
    my $bondStrB = $bH->[$i+1];
    $bondStrB = $map->{$bondStrB}->[0];
    my $sizeA = $map->{$bH->[$i]}->[2];my $sizeB = $map->{$bH->[$i+1]}->[2];
	my $ra=$connect->[$map->{$bH->[$i]}->[1]];my $rb=$connect->[$map->{$bH->[$i+1]}->[1]];
	my ($ia,$ta) = ($sizeA+getAtomAbsoluteIndex($ra,$bondStrA)
		       ,getAtomBType($ra,$bondStrA));
	my ($ib,$tb) = ($sizeB+getAtomAbsoluteIndex($rb,$bondStrB)
		       ,getAtomBType($rb,$bondStrB));
    my $if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");	
	push(@tempArr,pdl($ia,$ib,$if));
	}
	$connBondFunctionals{$counter}=cat(@tempArr);
	@tempArr=();
	## ANGLES ##
	foreach my $angs(@{$angH})
	{
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
                my $na;my $nb;my $nc;
		my $ra;my $rb;my $rc;
		my $sizeA; my $sizeB;my $sizeC;
		my($a,$b,$c) = split("-",$angs);
		##[AtomName,ResidueIndex,prevSize]##
		$na = $map->{$a}->[0];$ra = $connect->[$map->{$a}->[1]];
 		$nb = $map->{$b}->[0];$rb = $connect->[$map->{$b}->[1]];
		$nc = $map->{$c}->[0];$rc = $connect->[$map->{$c}->[1]];
		$sizeA=$map->{$a}->[2];$sizeB=$map->{$b}->[2];
		$sizeC=$map->{$c}->[2];

		($ia,$ta) = ($sizeA+getAtomAbsoluteIndex($ra,$na),getAtomBType($ra,$na));
		($ib,$tb) = ($sizeB+getAtomAbsoluteIndex($rb,$nb),getAtomBType($rb,$nb));
		($ic,$tc) = ($sizeC+getAtomAbsoluteIndex($rc,$nc),getAtomBType($rc,$nc));	
        	my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
        	push(@tempArr,pdl($ia,$ib,$ic,$if));		
	}
		$connAngleFunctionals{$counter} = cat(@tempArr);
		@tempArr = ();


	## DIHEDRALS ##
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

		($ia,$ta) = ($sizeA+getAtomAbsoluteIndex($ra,$na),getAtomBType($ra,$na));
		($ib,$tb) = ($sizeB+getAtomAbsoluteIndex($rb,$nb),getAtomBType($rb,$nb));
		($ic,$tc) = ($sizeC+getAtomAbsoluteIndex($rc,$nc),getAtomBType($rc,$nc));
		($id,$td) = ($sizeD+getAtomAbsoluteIndex($rd,$nd),getAtomBType($rd,$nd));	
        	
		## Adjust args for getEnergyGroup() ##
        ($nb,$nc) =  ($map->{$b}->[1]-$map->{$c}->[1]==0)?($nb,$nc):("nb?",$nc);
		my $eG = getEnergyGroup($rb,$rc,$nb,$nc);
		my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);	
        $eG = $eGRevTable{$eG};
		## [x,y,z,func,countDihedrals,energyGroup]
		push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);

	}
	## DOESN'T SUPPORT IMPROPER ##
	$connDiheFunctionals{$counter} = pdl(@tempArr);
	@tempArr = ();

}


sub returnBondTypeFromIndex
{
  my($idx) = @_;
  my $residue = $allAtoms{$idx}->[5];
  my $atom = $allAtoms{$idx}->[3];
  if(!$residue || !$atom)
  {confess("Error finding the residue for atom $idx. Perhaps your indices are wrong?\n");}
  if(!$residues{$residue}->{"atoms"}->{$atom})
  	{confess("$atom is not part of $residue\n");}
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
 	{confess("Ad-hoc Improper Create Error: Atom $a is part of neither residue $resIDA or $resIDB\n");}
 	if(returnResidueIndexFromIndex($b)!=$resIDA && returnResidueIndexFromIndex($b)!=$resIDB)
 	{confess("Ad-hoc Improper Create Error: Atom $b is part of neither residue $resIDA or $resIDB\n");}
 	if(returnResidueIndexFromIndex($c)!=$resIDA && returnResidueIndexFromIndex($c)!=$resIDB)
 	{confess("Ad-hoc Improper Create Error: Atom $c is part of neither residue $resIDA or $resIDB\n");}
 	if(returnResidueIndexFromIndex($d)!=$resIDA && returnResidueIndexFromIndex($d)!=$resIDB)
 	{confess("Ad-hoc Improper Create Error: Atom $d is part of neither residue $resIDA or $resIDB\n");}
 	
	my($ia)= (returnResidueIndexFromIndex($a)==$resIDA ? (getAtomAbsoluteIndex($resA,$iia)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$iia)));
 	my($ib)= (returnResidueIndexFromIndex($b)==$resIDA ? (getAtomAbsoluteIndex($resA,$iib)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$iib)));
	my($ic)= (returnResidueIndexFromIndex($c)==$resIDA ? (getAtomAbsoluteIndex($resA,$iic)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$iic)));
	my($id)= (returnResidueIndexFromIndex($d)==$resIDA ? (getAtomAbsoluteIndex($resA,$iid)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$iid)));
 	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
 	## [x,y,z,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
 	push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);
  	
}




sub connCreateInteractionsBOND
{

    my($consecResiduesH,$sizeA,$counter,$atomA,$atomB,$resAIdx,$resBIdx,$bEG,$imp) = @_;
	my @consecResidues = @{$consecResiduesH};
	my $residue = $consecResidues[1];


    ## AD-HOC BONDS ##
	my($angH,$diheH,$adjList,$bondStrA,$bondStrB)=createConnection($consecResiduesH,0,$atomA,$atomB);

	## BOND ##
	my ($ia,$ta) = (getAtomAbsoluteIndex($consecResidues[0],$bondStrA)
							,getAtomBType($consecResidues[0],$bondStrA));
	my ($ib,$tb) = ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$bondStrB)
							,getAtomBType($consecResidues[1],$bondStrB));
	
	
	my $if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");
	
    $connBondFunctionals{$counter}=pdl($ia,$ib,$if);
			
	## ANGLES ##
	my @tempArr;
	foreach my $angs(@{$angH})
	{
		my($a,$b,$c) = split("-",$angs);
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
		($ia,$ta) = ($a =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
				
		($ib,$tb) = ($b =~ /(.*)\?/)
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
				
		($ic,$tc) = ($c =~ /(.*)\?/) 
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
				
			
        my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
        push(@tempArr,pdl($ia,$ib,$ic,$if));

	}
        if(@tempArr)
        {$connAngleFunctionals{$counter} = cat(@tempArr);}
        else{#warn("PDB PARSE WARN:: There are no angles between ",$consecResidues[0]," and ",$consecResidues[1]);}
        }
		@tempArr = ();
			
			
			
			
	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
		my $eG;
		
		($ia,$ta) = ($a =~ /(.*)\?/)
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
				
		($ib,$tb) = ($b =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
				
		($ic,$tc) = ($c =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
			
		($id,$td) = ($d =~ /(.*)\?/) 
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$d),getAtomBType($consecResidues[0],$d));
		
		
		#if($consecResidues[0] eq "ASN" && $consecResidues[1] eq "NAG"
		#&& $b !~ /(.*)\?/ && $c !~ /(.*)\?/)
		#{
		#if($bEG=~ m/^$/){confess "HI";}
		#confess("$ia,$ib,$ic,$id,$if,1,$bEG, $consecResidues[0],$consecResidues[1] $b $c");
		#}
		
		
		if(!$bEG || ($b =~/.*\?/ && $c =~/.*\?/)|| ($b !~/.*\?/ && $c !~/.*\?/))
		{$eG=getEnergyGroup($consecResidues[0],$consecResidues[1],$b,$c);}
		else{$eG=$bEG;}
		my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);
		
		
		
		$eG = $eGRevTable{$eG};
		
		
		
		## [x,y,z,func,countDihedrals,energyGroup]
		push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);	
		
	}
	     
	   
		## Manually add Improper dihedrals ##
		if(scalar(@{$imp})!=0){
		appendImpropersBOND($consecResidues[0],$consecResidues[1],$resAIdx,$resBIdx,$sizeA,$imp,\@tempArr);
		}
		##CALL getSetDiheCounts####
		##getSetDiheCounts(\@tempArr);
		
        if(@tempArr)
        {$connDiheFunctionals{$counter} = pdl(@tempArr);}
        else{#warn("PDB PARSE WARN:: There are no dihedrals between ",$consecResidues[0]," and ",$consecResidues[1]);
            }
				@tempArr = ();
}


sub connCreateInteractions
{
    my($consecResiduesH,$sizeA,$counter,$startFlag,$atomA,$atomB,$bEG) = @_;
	my @consecResidues = @{$consecResiduesH};
	my $residue = $consecResidues[1];


    ## AD-HOC BONDS ##
	my($angH,$diheH,$adjList,$bondStrA,$bondStrB)=createConnection($consecResiduesH,$startFlag,$atomA,$atomB);

	## BOND ##
	my ($ia,$ta) = (getAtomAbsoluteIndex($consecResidues[0],$bondStrA)
							,getAtomBType($consecResidues[0],$bondStrA));
	my ($ib,$tb) = ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$bondStrB)
							,getAtomBType($consecResidues[1],$bondStrB));
	
	
	my $if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");
	
    $connBondFunctionals{$counter}=pdl($ia,$ib,$if);
			
	## ANGLES ##
	my @tempArr;
	foreach my $angs(@{$angH})
	{
		my($a,$b,$c) = split("-",$angs);
		my $ia;my $ib;my $ic;
		my $ta;my $tb;my $tc;
		($ia,$ta) = ($a =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
				
		($ib,$tb) = ($b =~ /(.*)\?/)
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
				
		($ic,$tc) = ($c =~ /(.*)\?/) 
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
				
			
        my $if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
        push(@tempArr,pdl($ia,$ib,$ic,$if));

	}
        if(@tempArr)
        {$connAngleFunctionals{$counter} = cat(@tempArr);}
        else{#warn("PDB PARSE WARN:: There are no angles between ",$consecResidues[0]," and ",$consecResidues[1]);}
        }
		@tempArr = ();
			
			
			
			
	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);
		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
		my $eG;
		
		($ia,$ta) = ($a =~ /(.*)\?/)
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$a),getAtomBType($consecResidues[0],$a));
				
		($ib,$tb) = ($b =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$b),getAtomBType($consecResidues[0],$b));
				
		($ic,$tc) = ($c =~ /(.*)\?/ )
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$c),getAtomBType($consecResidues[0],$c));
			
		($id,$td) = ($d =~ /(.*)\?/) 
			? ($sizeA+getAtomAbsoluteIndex($consecResidues[1],$1),getAtomBType($consecResidues[1],$1))
			: (getAtomAbsoluteIndex($consecResidues[0],$d),getAtomBType($consecResidues[0],$d));
		
		
		#if($consecResidues[0] eq "ASN" && $consecResidues[1] eq "NAG"
		#&& $b !~ /(.*)\?/ && $c !~ /(.*)\?/)
		#{
		#if($bEG=~ m/^$/){confess "HI";}
		#confess("$ia,$ib,$ic,$id,$if,1,$bEG, $consecResidues[0],$consecResidues[1] $b $c");
		#}
		
		
		if(!$bEG || ($b =~/.*\?/ && $c =~/.*\?/)|| ($b !~/.*\?/ && $c !~/.*\?/))
		{$eG=getEnergyGroup($consecResidues[0],$consecResidues[1],$b,$c);}
		else{$eG=$bEG;}
		my $if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);
		
		
		
		$eG = $eGRevTable{$eG};
		
		
		
		## [x,y,z,func,countDihedrals,energyGroup]
		push(@tempArr,[$ia,$ib,$ic,$id,$if,1,$eG]);	
		
	}
	     
	   
		## Manually add Improper dihedrals ONLY FOR NONBOND Connections##
		if(!$atomA || !$atomB){
		appendImpropers($consecResidues[0],$consecResidues[1],$sizeA,\@tempArr,$startFlag);
		}
		##CALL getSetDiheCounts####
		##getSetDiheCounts(\@tempArr);
		
        if(@tempArr)
        {$connDiheFunctionals{$counter} = pdl(@tempArr);}
        else{#warn("PDB PARSE WARN:: There are no dihedrals between ",$consecResidues[0]," and ",$consecResidues[1]);
            }
				@tempArr = ();
}

sub appendImpropers
{
 my($resA,$resB,$sizeA,$tempArr,$startFlag) = @_;
 my $resAIp = $residues{"$resA"}->{"impropers"};
 my $resBIp = $residues{"$resB"}->{"impropers"};
 my @connImproper; my $connHandle;
 
 ## WORK RESIDUE B ##
 foreach my $ips(@{$resBIp})
 {
    if(! (defined $ips) ) {next;}
	my($a,$b,$c,$d) = @{$ips};
	my($ia,$ta)=($sizeA+getAtomAbsoluteIndex($resB,$a),getAtomBType($resB,$a));
	my($ib,$tb)=($sizeA+getAtomAbsoluteIndex($resB,$b),getAtomBType($resB,$b));
	my($ic,$tc)=($sizeA+getAtomAbsoluteIndex($resB,$c),getAtomBType($resB,$c));
	my($id,$td)=($sizeA+getAtomAbsoluteIndex($resB,$d),getAtomBType($resB,$d));
	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
	## [a,b,c,d,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
	push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
 }
 
 ## WORK ON INTER-RESIDUAL IMPROPERS ##
 ### CHANGE THIS, ONLY HANDLES SINGLE IMPROPERS ###
 $connHandle = $connections{$residues{$resA}->{"type"}}->{$residues{$resB}->{"type"}};
 #@connImproper = @{$connHandle->{"improper"}};
 foreach my $ips(@{$connHandle->{"improper"}})
 {
  	my @atomsHandle = @{$ips->{"atom"}}; 
 	my ($a,$b,$c,$d) = ('O','CA','C','N?');
        ($a,$b,$c,$d) = @atomsHandle;
 	my($ia,$ta)= ($a !~/(.*)\?/ ? (getAtomAbsoluteIndex($resA,$a),getAtomBType($resA,$a)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$1),getAtomBType($resB,$1)));
 	my($ib,$tb)= ($b !~/(.*)\?/ ? (getAtomAbsoluteIndex($resA,$b),getAtomBType($resA,$b)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$1),getAtomBType($resB,$1)));
 	my($ic,$tc)= ($c !~/(.*)\?/ ? (getAtomAbsoluteIndex($resA,$c),getAtomBType($resA,$c)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$1),getAtomBType($resB,$1)));
 	my($id,$td)= ($d !~/(.*)\?/ ? (getAtomAbsoluteIndex($resA,$d),getAtomBType($resA,$d)) 
	: ($sizeA+getAtomAbsoluteIndex($resB,$1),getAtomBType($resB,$1)));
 	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
 	## [x,y,z,func,countDihedrals,energyGroup] energyGroup is negative signifies improper
 	push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);
  }	



 if($startFlag == 0) {return;}
 
 ## WORK RESIDUE A ##
 foreach my $ips(@{$resAIp})
 {
    if(! (defined $ips) ) {next;}
	my($a,$b,$c,$d) = @{$ips};
	my($ia,$ta)=(getAtomAbsoluteIndex($resA,$a),getAtomBType($resA,$a));
	my($ib,$tb)=(getAtomAbsoluteIndex($resA,$b),getAtomBType($resA,$b));
	my($ic,$tc)=(getAtomAbsoluteIndex($resA,$c),getAtomBType($resA,$c));
	my($id,$td)=(getAtomAbsoluteIndex($resA,$d),getAtomBType($resA,$d));
	my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");
	## [x,y,z,func,countDihedrals,energyGroup]
	push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
 }



}

sub connWildcardMatchAngles
{
    my($a,$b,$c) = @_;
	my $angHandle = $interactions->{"angles"};
	my $funct="";
	## WILD CARD MATCHING CONDITIONALS ##
	my $matchScore = 0; my $saveScore = 0;
	foreach my $matches(keys %{$angHandle})
	{
		$matchScore = 0;
		my ($aM,$bM,$cM) = split("-",$matches);
		if(($a !~ /\Q$aM\E/ && $aM !~ /\Q*\E/)
			|| ($b !~ /\Q$bM\E/ && $bM !~ /\Q*\E/)
			|| ($c !~ /\Q$cM\E/ && $cM !~ /\Q*\E/)){next;}
		if($a =~ /\Q$aM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($b =~ /\Q$bM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($c =~ /\Q$cM\E/) {$matchScore+=2;} else {$matchScore+=1;}
		if($matchScore >= $saveScore)
		{$saveScore = $matchScore;$funct = $angHandle->{$matches};}
	}
	if(!defined $funct|| $funct eq ""){confess("\nINTERACTION GENERATE ERROR\n There is no function for bType combination $a-$b-$c. Check .b file\n");}
	return $funct;
}

sub connWildcardMatchBond
{
	my($typeA,$typeB) = @_;
	my $funct="";
    ## If both bond types exists ##
	if(exists $interactions->{"bonds"}->{$typeA}
		&& exists $interactions->{"bonds"}->{$typeA}->{$typeB})
	{$funct = $interactions->{"bonds"}->{$typeA}->{$typeB};}
		
	## If typeA exists while TypeB is a wildcard ##
	elsif (exists $interactions->{"bonds"}->{$typeA}
		&& !(exists $interactions->{"bonds"}->{$typeA}->{$typeB}))
	{
        $funct = $interactions->{"bonds"}->{$typeA}->{"*"};
        if(!defined $funct || $funct eq "")
        {$funct = $interactions->{"bonds"}->{"*"}->{"*"};}
    }
		
	## If typeA is a wildcard while TypeB exists ##
	elsif(!(exists $interactions->{"bonds"}->{$typeA}) &&
			(exists $interactions->{"bonds"}->{"*"}
            && exists $interactions->{"bonds"}->{"*"}->{$typeB})
            )
	{
        $funct = $interactions->{"bonds"}->{"*"}->{$typeB};
        if(!defined $funct || $funct eq "")
        {$funct = $interactions->{"bonds"}->{"*"}->{"*"};}

    }
		
	## If both types are wildcard ## 
	else
	{
        $funct = $interactions->{"bonds"}->{"*"}->{"*"};
    }
    if(!defined $funct || $funct eq ""){confess("\nINTERACTION GENERATE ERROR\n There is no function for bType combination $typeA-$typeB. Check .b file\n");}
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
	if(!defined $funct || $funct eq ""){confess("\nINTERACTION GENERATE ERROR\n There is no function for bType combination $a-$b-$c-$d. Check .b file\n");}
	return $funct;
}


sub connWildcardMatchDihes
{
    my($a,$b,$c,$d,$eG) = @_;
	my $diheHandle = $interactions->{"dihedrals"}->{"$eG"};
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
	if(!defined $funct || $funct eq ""){confess("\nINTERACTION GENERATE ERROR\n There is no function for bType combination $a-$b-$c-$d
	with energyGroup=$eG. Check .b file\n");}
	return $funct;
}



sub createMultiConnections
{
  my($connect) = @_;
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

    #print Dumper %union; print Dumper $resABonds;exit;
	## Start of chain no inter residue connection ##
    #  but setup leftAtom and leftResidue sizes #
    if($i == 0) 
	{
	  $connHandle 
      = $connections{$residues{$connect->[0]}->{"type"}}->{$residues{$connect->[1]}->{"type"}};
	  $leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
	  $leftAtom = $bondMapHashRev{"$leftAtom-$i"};
	  $prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});
       next;
	}
        ## $i > 0, create inter residue connection ##
	## $i-1 <--> $i
	$connHandle = $connections{$residues{$connect->[$i-1]}->{"type"}}->{$residues{$connect->[$i]}->{"type"}};
	$rightAtom = $connHandle->{"bond"}->[0]->{"atom"}->[1];
    $rightAtom = $bondMapHashRev{"$rightAtom-$i"};
	push(@connectList,$leftAtom);
	push(@connectList,$rightAtom);
    
    ## $i <--> $i+1
        if($i == $#$connect) {last;}
        $connHandle = $connections{$residues{$connect->[$i]}->{"type"}}->{$residues{$connect->[$i+1]}->{"type"}};
    $leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
    $leftAtom = $bondMapHashRev{"$leftAtom-$i"};
    $prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});

  }
   ## Create Inter residue connection ##
   for($i=0;$i<scalar(@connectList)-1;$i+=2) {
  	push(@{$union{$connectList[$i]}},$connectList[$i+1]);
	push(@{$union{$connectList[$i+1]}},$connectList[$i]);
        
   }
   ($dihes,$angles,$oneFour)=adjListTraversal(\%union);
   return (\@connectList,$angles,$dihes,\%bondMapHash);

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
    $connHandle = $connections{$residues{$connect->[0]}->{"type"}}->{$residues{$connect->[1]}->{"type"}};
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
  my($fileName,$fileName2,$noChainFlag,$noAllFlag,$coarseGraining) = @_;
  my $numContacts = 0; my $garbage = 0;
  my $line = "";
  my $type1;my $type2; my $contact1; my $contact2;
  my $dist;
  my @interiorTempPDL;
		## noChainFlag ignores contact between chains ##  

  ## OPEN .contact FILE ##
  unless (open(MYFILE, $fileName)) {
    confess "CONTACT FILE READ ERROR: Contact file $fileName doesn't exist, check shadow.log for possible error from Shadow Map program\n";
  }
  
  if($noAllFlag){print "NOTE::Ignoring contacts calculated from Shadow Map Program\n";}

  if(!$noAllFlag){
  while($line = <MYFILE>)
  {
	($contact1,$type1,$contact2,$type2,$dist) = split(/\s+/,$line);
 	## INTER CHAIN CONTACT IGNORE ##
 	if($noChainFlag && ($allAtoms{$type1}->[4] ne $allAtoms{$type2}->[4]))
 	{next;}
	## NOTE RESIDUE INDEX == CONTACT ##
	if($dist < $interactionThreshold->{"contacts"}->{"shortContacts"})
	{confess("ERROR: CONTACTS between atoms $type1 $type2 exceed contacts threshold with value $dist\n");}
	
	push(@interiorTempPDL,[0,$type1,0,$type2,$dist]);
	$numContacts++;
  }
    print "Adding additional contacts from $fileName2\n";

   ## NO DCA FILE RETURN ##
   if($fileName2 eq ""){$contactPDL = pdl(@interiorTempPDL);return $numContacts;}
  }
   ## Else Proceed ##
   ## NO DCA FILE RETURN ##
   if($fileName2 eq ""){$contactPDL = pdl(@interiorTempPDL);return $numContacts;}
  
  print "Adding additional contacts from $fileName2\n";
  ## OPEN .dca FILE ##
  unless (open(MYFILE1, $fileName2)) {
    confess "CONTACT FILE READ ERROR:Cannot read additional contact file '$fileName2'.\nProgram closing.\n";
  }
  
   ## DCA CONTACTS ARE IN FORM
   ## atom1 atom2 dist(angstrom) ##
   while($line = <MYFILE1>)
  {
	my ($contact1,$contact2,$dist,$epsilon) = split(/\s+/,$line);
 	## INTER CHAIN CONTACT IGNORE ##
 	#if($noChainFlag && ($allAtoms{$type1}->[2] ne $allAtoms{$type2}->[2])){next;}
	## NOTE RESIDUE INDEX == CONTACT ##
	## ANGSTROM TO NM CONVERSION ##
	$dist/=10;
	##[dcaContacts_boolean,c1,value,c2,dist]##
    if(!exists $allAtoms{$contact1}){warn("ATOM $contact1 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}
    if(!exists $allAtoms{$contact2}){warn("ATOM $contact2 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}
    push(@interiorTempPDL,[1,$contact1,$epsilon,$contact2,$dist]);
	$numContacts++;
  }
  
  $contactPDL = pdl(@interiorTempPDL);
  return $numContacts;
  
}


1;
