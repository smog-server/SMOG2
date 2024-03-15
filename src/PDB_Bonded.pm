#########################################################################################
#
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#                     Jeffrey Noel, Mariana Levi, Mohit Ranghunathan,
#                 Heiko Lammert, Ryan Hayes, Jose Onuchic & Paul Whitford
#
#          Copyright (c) 2015,2016,2018,2021,2022 The SMOG development team at
#                      The Center for Theoretical Biological Physics
#                       Rice University and Northeastern University
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

##############################################################################
# PDB_Bonded: A set of routines for parsing the PDB file and generating 
# bonded interaction information, residue and coordinate info.
##############################################################################
package PDB_Bonded;

use templateParser;
use setRatios;
use strict;
use warnings FATAL => 'all';
use Exporter;
use PDL; ## LOAD PDL MODULE
use Storable qw(dclone);
use smog_common;
## DECLARATION TO SHARE DATA STRUCTURES ##
our @ISA = 'Exporter';
our @EXPORT = 
qw(%eGTable $energyGroups $interactionThreshold %fTypes %residues $termRatios %allAtoms parseCONTACT $contactPDL catPDL $totalAtoms returnFunction intToFunc funcToInt %bondFunctionals %AngleData %DihedralData %BondData %resPDL %isBOND %bondPDL %dihedralFunctionals %angleFunctionals setInputFileName parseBif parseSif parseBonds createBondFunctionals createDihedralAngleFunctionals parseNonBonds getContactFunctionals $contactSettings $interactions clearPDBMemory clearBifMemory parsePDBATOMS checkPDB);

our %tempPDL = ();
our %tempPDLB = ();
# we call this "PDL" so that it is clear this lists follows tempPDL, but it isn't actually PDL
# it hold the bTypes of atoms in tempPDLB, which is saved in case we have a BOND
our %tempPDLbType = ();
our %resPDL;
our %bondPDL;
our %isBOND;
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
# keep track of BONDs, so that contacts can be removed
our %ignoreBONDed;
###########################
## CLEAR VARIABLE MEMORY ##
###########################
sub clearPDBMemory {
	undef %tempPDL;
	undef %tempPDLB;
	undef %tempPDLbType;
	undef %resPDL;
	undef %bondPDL; 
	undef %isBOND; 
	undef %AngleData;
	undef %DihedralData;
	undef %BondData;
	undef $totalAtoms;
	#undef $contactPDL; 
	undef %indexMap;
	#for coarse graining contact maps we need to know the atomNum<->resNum mapping
	%allAtomsBackup = %{ dclone (\%allAtoms) };
	undef %allAtoms; 
}

####################################################################
# checkPDB
# this routine will check to make sure the PDB conforms to all
# expected formatting requirements
####################################################################
sub checkPDB
{
	# load the PDB into memory and return a handle to the array.
	print "Checking PDB formatting\n";	
	my ($fileName,$CGenabled,$CGstage,$freecoor) = @_;
	
	## INTERNAL VARIABLES ##
	my $counter = 0;
	my $x;my $y;my $z;
	my $residue; my $atom;my $atomSerial;
	my $atomsInRes; 
	my $i; 
	my $nbType;
	my $residueType;
	my $pairType;
	my $atomCounter=0;
	my $chainNumber = 0;
	my $residueIndex=1;
	my $interiorPdbResidueIndex=0;
	my $lineNumber = 0;
	my $lastchainstart=0;
	my $endfound=0;
	my $residueSerial=0;
	my $largebase=0;
	my %duplicateatoms;
	## OPEN .PDB FILE ##

	unless (open(PDBFILE, $fileName)) {
		smog_quit ("Cannot read from '$fileName'.");
	}

	my $lastrecord="";
	# first check and make sure there is an END and there are no ATOM lines afterwards
	my @PDBDATA=();
	while(<PDBFILE>)
	{
		push(@PDBDATA,$_);
	}
	close(PDBFILE);
	my $lastline="";
	foreach my $record(@PDBDATA){
		my $lng = $record;
		chomp($lng);
		if($lng eq ""  && $endfound!=1 && $lastline ne ""){
			smog_quit("Blank line found before END.  This is generally not intentional. Blank line appears directly after:\n$lastline");		
		}
		$lastline=$record;	
		if($record =~m/^END/){
			$endfound=1;
			next;
		}
		if($lng eq "" || $record =~ m/^remark/i || $record =~ m/^HEADER/ || $record =~ m/^TITLE/){
			next;
		}

		if($record =~ m/^LARGE/){
			# large-base numbering will be used
			$largebase=1;
			&InitLargeBase;
			my $B=$smog_common::BaseN;
			print "PDB is using large numbers (expecting base $B)\n";
			next;
		}


		if($record !~ m/^BOND/ && $endfound ==1){
			smog_quit("PDB format issue: Only user-defined bonds given by BOND, or REMARK lines, may be listed after END. Offending line: \"$record\"\n");
		}
	}

        unless($endfound){smog_quit("PDB format error. END must appear at the end of the ATOM lines.")};

	$endfound=0;
	my $lastresindex="null";
        my %BONDlist;
	 ## LOOP THROUGH EACH LINE ##
	for (my $K=0;$K<=$#PDBDATA;$K++){
		my $record=$PDBDATA[$K];
 		# when reading the file, check to see if the line is a BOND, then check if TER/END
		# if neither, then read each residue and make sure they are only ATOM/HETATM
		$lineNumber++;
	
		## PARSE BOND LINES ##
		my $lng = $record;
		chomp($lng);
		$lng =~ s/\s+//g;	
		$lng =~ s/\t+//g;	
		if($record =~  m/^remark/i || $record =~ m/^LARGE/ || $lng eq ""|| $record =~ m/^HEADER/ || $record =~ m/^TITLE/){
			next;
		# make sure BOND appears after END
		}elsif($record !~ m/^BOND/ && $endfound ==1){
			smog_quit("PDB format issue: Only user-defined bonds given by BOND, or REMARK lines, may be listed after END. Offending line: \"$lng\"\n");
		}

		# CHECK IF IT IS A BOND
 		if($record =~ m/^BOND/){

			if($endfound ==0){
   				smog_quit("PDB format issue: User-defined bonds given by BOND should be listed immediately after END.");
  			}

    			chomp($record);
			unless($CGenabled == 1 && $CGstage==0){
				#only skip if this is the first round of CG
				my $resat;
				if($CGenabled==1){
					$resat="residue";
				}else{
					$resat="atom";
				} 
				my @TMP = split(/\s+/,$record);
				if(@TMP != 6){
					smog_quit("Tag BOND must have 5 arguments (chain A, $resat number A, chain B, $resat number B, energy group). Offending line:\n$record");
				}
				my($trig,$chaina,$atoma,$chainb,$atomb,$eG) = split(/\s+/,$record);
				if(exists $BONDlist{"$chaina $atoma $chainb $atomb"}){
					smog_quit("The following BOND is defined more than once in the PDB file (chainA, atomnumA, chainB, atomnumB): $chaina, $atoma, $chainb, $atomb\nIf you suppress this error, the bond (and associated angles) will appear multiple times in the force field (.top file).")
				}
				$BONDlist{"$chaina $atoma $chainb $atomb"}=1;
				$BONDlist{"$chainb $atomb $chaina $atoma"}=1;
				
				#internally, chains are indexed 0,1...
				$chaina--;
				$chainb--;
				if(!exists $indexMap{"$chaina-$atoma"}){
					$chaina++;
					smog_quit("When trying to add a user-defined BOND, could not find $resat $atoma in chain $chaina");
				}
				if(!exists $indexMap{"$chainb-$atomb"}){
					$chainb++;
					smog_quit("When trying to add a user-defined BOND, could not find $resat $atomb in chain $chainb");
				}
			}
			$counter++;
			next;
		}

		## CHECK FOR END OR TER  ##
		if($record =~ m/^TER|^END/)
		{
 			$lastresindex="null";
			$chainNumber++; ## INCREMENT CHAIN NUMBER ##
			if(!defined $atomSerial){
				smog_quit("TER or END line encountered prior to the first ATOM. This is probably a mistake in your PDB file.");
			}
			my $chainlength=$atomSerial-$lastchainstart;
			$counter++;
			$lastchainstart=$atomSerial;
			if($record =~ m/^END/){$endfound=1;}
			next;
		} 
	
		## CHECK IF ATOM/HETATM ##
		if($record =~ m/^ATOM/ || $record =~ m/^HETATM/)
		{
	        	$lineNumber--;
		 	## OBTAIN RESIDUE NAME ##
			$residue = trim(substr($record,17,4));
			if(!exists $residues{$residue}){smog_quit (" \"$residue\" doesn't exist in .bif. See line $lineNumber of PDB file.");}
			## if first iteration, save residueBackup, and use residues
			if($CGstage == 1){
				$atomsInRes = scalar(keys(%{$residueBackup{$residue}->{"atoms"}}));
			}else{
				$atomsInRes = scalar(keys(%{$residues{$residue}->{"atoms"}}));
			}
			my $atomsInBif=scalar(keys(%{$residues{$residue}->{"atoms"}}));
			if($atomsInBif != 1 && $CGstage ==1){
				smog_quit ("When using CG, each residue can only have one atom in the CG template. Check .bif definition for $residue");
			}
			my $atomsmatch=0;
			my $resname=$residue;
	        	my $resindex = trim(substr($record,22,5));
			if ($lastresindex ne "null"){
				my $diff;
				if($largebase==1){
					$diff=BaseLargetoTen($resindex)-BaseLargetoTen($lastresindex);
				}else{
					$diff=$resindex-$lastresindex;
				}

				if($diff != 1 && $resindex ne $lastresindex){
					smog_quit("Non-sequential residue numbers ($lastresindex,$resindex) appear at line $lineNumber.");
				}
			}
			$lastresindex=$resindex;
			my %uniqueAtom;
			$residueSerial++;
			for($i=0;$i<$atomsInRes;$i++)
			{
				my $record = $PDBDATA[$K];
				$K++;
	 			$lineNumber++;
				if($record =~ m/^ATOM|^HETATM/)
				{
	
		   			$residue = trim(substr($record,17,4));
		   			my $altlocator = substr($record,16,1);
					if($altlocator ne " "){
						smog_quit("Alternate location indicator found at line $lineNumber.  Alt. Loc. Indic. not supported by SMOG.");
					}
	
		        		if(!exists $residues{$residue}){smog_quit (" \"$residue\" doesn't exist in .bif. See line $lineNumber of PDB file.");}
	            			$interiorPdbResidueIndex = trim(substr($record,22,5)); 
				}
				if($resname ne $residue or $resindex ne $interiorPdbResidueIndex or  $record !~ m/^ATOM|^HETATM/){
					my $linetemp=$lineNumber-1;
					my $missingatoms="";
				        if($CGstage == 1){
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
					smog_quit("It appears that a residue in the PDB file does not contain all of the atoms defined in the .bif file.\nOffending residue: $resname (ending at line $linetemp).  Missing atoms: $missingatoms");	
					# these next lines only matter if errors are being reported as warnings
		 			$K--;
	        			$lineNumber--;
					last;
				}
	
				unless($interiorPdbResidueIndex =~ /^\d+$/ ||  $largebase ==1){;
					smog_quit ("Residue $residue$interiorPdbResidueIndex contains non integer value for the index, or an insertion code.");
				}
	
				## CHECK IF ALL ATOMS CONFORM TO BIF RESIDUE DECLARATION ##
				$atom = trim(substr($record, 12, 4));

				if(exists $uniqueAtom{$atom})
				{
					smog_quit("$atom appears twice in $residue at line $lineNumber\n");
				}
				else {
					$uniqueAtom{$atom}=1;
				}
					
				if(!exists $residues{$residue}->{"atoms"}->{$atom})
				{
					if($CGstage == 0){
						smog_quit ("\"$atom\" doesn't exist in .bif declaration of \"$residue\"");
					}else{
						# if CG is turned on, skip atom
						next;
					}
				}
				
				$atomsmatch++;
                                if(defined $freecoor){
					# Read the PDB coordinates as free-format.
					my $string=trim(substr($record, 30));
					my @coor=split(/\s+/,$string);
					$x=$coor[0];
					$y=$coor[1];
					$z=$coor[2];
				}else{
					if(substr($record,34,1) !~  m/\./ ) {
						smog_quit("X coordinate in PDB file is not properly formatted.  The decimal should be column 35. Problematic line:\n$record");
 					}
					if(substr($record,42,1) !~  m/\./ ) {
						smog_quit("Y coordinate in PDB file is not properly formatted.  The decimal should be column 43. Problematic line:\n$record");
 					}
					if(substr($record,50,1) !~  m/\./ ) {
						smog_quit("Z coordinate in PDB file is not properly formatted.  The decimal should be column 51. Problematic line:\n$record");
 					}
					$x = trim(substr($record, 30, 8));
					$y = trim(substr($record, 38, 8));
					$z = trim(substr($record, 46, 8));
					if(whatAmI($x) > 2 || whatAmI($y) > 2 || whatAmI($z) > 2){
						smog_quit("Coordinate read, but does not appear to be a number. Perhaps you are using free-formatted coordinates and should employ the -freecoor flag. Issue found at line:\n$record");
					}
				}
				$atomCounter++;
				$atomSerial=$atomCounter;

				if(!exists $duplicateatoms{"$x$y$z"}){
					$duplicateatoms{"$x$y$z"}=1;
				}else{
					smog_quit("Atom has identical coordinates as a previous atom. Error found at line:\n$record");
				}			
				my $th=$residues{$residue}->{"atoms"}->{$atom};	
				$nbType = $th->{"nbType"};
				$pairType = $th->{"pairType"};
				$residueType = $residues{$residue}->{"residueType"};
				my $pdbIndex;
				if($CGstage==1){
					$pdbIndex = $interiorPdbResidueIndex;
				}else{
					$pdbIndex = trim(substr($record,6,5));
				}
	            		my $resIDorig = trim(substr($record,22,5)); 
				if(exists $indexMap{"$chainNumber-$pdbIndex"}){
					smog_quit("Atom/Residue numbers must be unique within each chain. Offending line:\n$record");
				}
				$indexMap{"$chainNumber-$pdbIndex"}=$atomSerial;
				# the atoms are now being stored in checkPDB, not parsePDB	
				$allAtoms{$atomSerial}=[$nbType,$residueType,$residueIndex,$atom,$chainNumber,$residue,$x,$y,$z,$residueSerial,$pairType,$pdbIndex,$resIDorig];
			}
			$K--;
			if(!defined $residues{$residue}->{"atomCount"}){
				smog_quit("Residue $residue found in PDB, but it does not appear to be defined in the .bif file.  This issue is usually detected when using smog_adjustPDB. Issue found at line:\n$record");

			}
			## CONCAT RESIDUE ##
			$residueIndex++;
				
		}else{
			smog_quit("Expected ATOM or HETATM line at line $lineNumber.");
		}
		$lastrecord=$record;
		$record = "";
 	}
	print "Done checking PDB formatting\n\n";
	return (\@PDBDATA);
}



####################################################################
# parsePDBATOMS
# this routine will parse the PDB atoms and build the covalent
# geometry
####################################################################
sub parsePDBATOMS
{
	
	my ($PDBDATA,$CGenabled,$CGstage,$freecoor) = @_;
	my @PDBDATA=@{$PDBDATA};
	## INTERNAL VARIABLES ##
	my $counter = 0;
	my @temp; my @tempB;my @union;
	my @consecResidues;
	my $x;my $y;my $z;
	my $residue; my $atom;my $atomSerial;
	my $atomsInRes; 
	my $i; my $putIndex=0; 
	$totalAtoms = 0;
	my $nbType;
	my $bType;
	my $residueType; 
	my $pairType;
	my $atomCounter=0;
	my $chainNumber = 0;
	my $residueIndex=1;
	my $interiorPdbResidueIndex=0;
	my $lineNumber = 0;
	my $lastchainstart=0;
	my $residueSerial=0;
	## OPEN .PDB FILE ##
	my $atominhash;
	my $PDBresstart;
	my %bondlists;
	my %bondlistsEG;
	my $BONDNUM = 0;
	print "Organizing PDB data\n";
	if($CGstage == 1){
		$atominhash=\%residueBackup;
		$PDBresstart=22;
	}else{
		$atominhash=\%residues;
		$PDBresstart=6;
	}
	# first check and make sure there is an END and there are no ATOM lines afterwards
	my $lastrecord="";

	my $lastresindex="null";
	 ## LOOP THROUGH EACH LINE ##
        for (my $K=0;$K<=$#PDBDATA;$K++){
                my $record=$PDBDATA[$K];
		$lineNumber++;
	
		my $lng = $record;
		chomp($lng);
		$lng =~ s/\s+//g;	
		$lng =~ s/\t+//g;	
		if($record =~ m/^remark/i || $lng eq "" || $record =~ m/^LARGE/ || $record =~ m/^HEADER/ || $record =~ m/^TITLE/){
			next;
		}
 		if($record =~ m/^BOND/){

			if($CGenabled==1 && $CGstage==0){
				#we are doing the all-atom processing before CGing, so ignore the BONDs. Otherwise, use them.
				next;
			}
    			chomp($record);
   
			my @TMP = split(/\s+/,$record);
			my($trig,$chaina,$atoma,$chainb,$atomb,$eG) = split(/\s+/,$record);

			my $userchaina=$chaina;	
			my $userchainb=$chainb;	
			#internally, chain numbers begin at 0, not 1.
			$chaina--;
			$chainb--;
			# get information about the two BONDed atoms
			my $idxA = $indexMap{"$chaina-$atoma"};
			my $idxB = $indexMap{"$chainb-$atomb"};
			my $resA = $allAtoms{$idxA}->[5];
			my $resB = $allAtoms{$idxB}->[5];
			my $atomA = $allAtoms{$idxA}->[3];
			my $atomB = $allAtoms{$idxB}->[3];
			my $resAIdx = $allAtoms{$idxA}->[2];
			my $resBIdx = $allAtoms{$idxB}->[2];
			my $resAIdxOrig = $allAtoms{$idxA}->[12];
			my $resBIdxOrig = $allAtoms{$idxB}->[12];
			my $fragment;
			# keep a record of each added bond, so that any contacts associated with this pair can be pruned, later.
			smog_note("Generating user-specified bonded interaction between atoms $atomA in $resA$resAIdxOrig of chain $userchaina (residue $resAIdx in gro file) and $atomB in $resB$resBIdxOrig of chain $userchainb (residue $resBIdx in gro file).\nWill assign to energy group $eG.");


			# make a hash of res indexes that we want to include in the fragment
			# include the two nearest residues, in sequence, for each of the two sides
			# this will ensure that enough information is available to reconstruct
			# any possible dihedral that would be associated with the BOND
			my %residFragHash;
			foreach my $RI($resAIdx,$resBIdx){
				for (my $II=$RI-2;$II<=$RI+2;$II++){
					$residFragHash{$II}=1;
				}
			}
			delete($residFragHash{$resAIdx});
			$fragment=$tempPDLB{$resAIdx};
			# create the fragment and all bond lists 
			# format [x y z globalnumber]
			my @bTypesFragment=@{$tempPDLbType{$resAIdx}};
			foreach my $kk(keys %residFragHash){
				if(exists $tempPDLB{$kk}){
					$fragment=($fragment)->glue(1,$tempPDLB{$kk});
					push(@bTypesFragment,@{$tempPDLbType{$kk}});
				}
			}

			my $newindex=0;
			# %keep is a hash that maps global pdl indexes to the fragment
			my %keep;
			foreach my $kk(list(slice($fragment,4))){
				$keep{$kk}=$newindex;
				$newindex++;
			}

			# find the positions of the newly-bonded atoms in the piddle 
			my $indexA;
			my $indexB;
			my @local2global;
			$newindex=0;
			foreach my $kk(list(slice($fragment,3))){
				push(@local2global,$kk);
				if($idxA==$kk){
					$indexA=$newindex;
				}
				if($idxB==$kk){
					$indexB=$newindex;
				}
				$newindex++;
			}

			# go through fragment, find keys that are in keep (i.e. the atoms in the fragment). Then iterate through each array and convert
			# bondlistFragment and energygroupsFragment are for storing the bonds and energy groups of bonds associated with the fragment that is being used for the BOND.  
			my %bondlistFragment;
			my %energygroupsFragment;
			for my $JK(keys %bondlists){
				if(exists $keep{$JK}){
					my $newkey=$keep{$JK};
					my @eGa=@{$bondlistsEG{$JK}};
					foreach my $KL(@{$bondlists{$JK}}){
						my $eGval=shift @eGa;
						if(exists $keep{$KL}){
							my $newval=$keep{$KL};
							push(@{$bondlistFragment{$newkey}},$newval);	
							push(@{$energygroupsFragment{$newkey}},$eGval);	
						}
					}
				}
			}


			# save the coordinates of this piddle, for later use. 
			$bondPDL{$counter}=$fragment;
			$BONDNUM++;
			$isBOND{$counter}=$BONDNUM;
			connCreateInteractionsSingleBOND(\%bondlistFragment,\%energygroupsFragment,\@bTypesFragment,$counter,$chaina,$chainb,$indexA,$indexB,$atomA,$atomB,$resA,$resB,$eG,\@local2global); 

			$counter++;
			next;
		}

		## IF TER LINE  ##
		if($record =~ m/^TER|^END/)
		{
 			$lastresindex="null";
			my $chainlength=$atomSerial-$lastchainstart;
			if($chainlength==0){
				next;
			}
			$chainNumber++; ## INCREMENT CHAIN NUMBER ##
			## CREATE INTERACTION ##
			print "Building covalent geometry for chain $chainNumber\n";
        		my ($bl,$blEG)=GenerateBondedGeometry(\@consecResidues,$counter,$chainNumber,$chainlength,$lastrecord);
			# keeping track of the bondlists, so that we can use them with the user-defined BOND routines
			foreach my $I(keys %{$bl}){
				my @tmparr=();
				my $In=$I+$lastchainstart+1;
				$bondlists{$In}=\@tmparr;
				foreach my $JJJ(@{${$bl}{$I}}){
                			# convert to global numbering
					push(@{$bondlists{$In}},$JJJ+$lastchainstart+1);
				}
				my @tmparr2=();
				$bondlistsEG{$In}=\@tmparr2;
				foreach my $JJJ(@{${$blEG}{$I}}){
					push(@{$bondlistsEG{$In}},$JJJ);
				}

			}

		   	$bondPDL{$counter}=pdl(@union);
			@union = ();$counter++;
        		@consecResidues = ();
			$lastchainstart=$atomSerial;
			next;
		} 
	
		## ONLY WORK WITH ATOM LINES ##
		if($record =~ m/^ATOM/ || $record =~ m/^HETATM/)
		{
	        	$lineNumber--;
		 	## OBTAIN RESIDUE NAME ##
			$residue = trim(substr($record,17,4));
			$atomsInRes = scalar(keys(%{${$atominhash}{$residue}->{"atoms"}}));
			my $atomsInBif=scalar(keys(%{$residues{$residue}->{"atoms"}}));
			my $atomsmatch=0;
			my $resname=$residue;
	        	my $resindex = substr($record,22,5);
			$lastresindex=$resindex;
			my %uniqueAtom;
			$residueSerial++;
			for($i=0;$i<$atomsInRes;$i++)
			{
				$record = $PDBDATA[$K];
				$K++;
	 			$lineNumber++;
		   		$residue = trim(substr($record,17,4));
		   		my $altlocator = substr($record,16,1);
	
				$atom = trim(substr($record, 12, 4));
				my $ahandle=$residues{$residue}->{"atoms"}->{$atom};	
				if(!defined $ahandle){
					# this must be something we are excluding because of CGing
					next;
				}	
				$atomsmatch++;

####  this info can be taken from allAtoms.  We don't need to reparse the PDB
                                if(defined $freecoor){
					# Read the PDB coordinates as free-format.
					my $string=trim(substr($record, 30));
					my @coor=split(/\s+/,$string);
					$x=$coor[0];
					$y=$coor[1];
					$z=$coor[2];
				}else{
					$x = trim(substr($record, 30, 8));
					$y = trim(substr($record, 38, 8));
					$z = trim(substr($record, 46, 8));
				}

				$atomCounter++;
				$atomSerial=$atomCounter;
				
				$putIndex = $ahandle->{"index"};
				$nbType = $ahandle->{"nbType"};
				$bType = $ahandle->{"bType"};
				$pairType = $ahandle->{"pairType"};
				$residueType = $residues{$residue}->{"residueType"};
				$temp[$putIndex]=[$x,$y,$z,$atomSerial];
				$tempB[$putIndex]=$bType;
			}
			$K--;
			# this is added so that we know the PDL index, as well as the global gro index
			for (my $I=0;$I<=$#temp;$I++){
				push(@{$temp[$I]},$I+1+$totalAtoms);
			}
			if($residues{$residue}->{"atomCount"} == -1){
				$totalAtoms+=$#temp+1;
			}else{
				$totalAtoms+=$residues{$residue}->{"atomCount"};
			}
			## CONCAT RESIDUE ##
			push(@union,@temp);
	        	$tempPDL{$residue}->{$residueIndex}=pdl(@temp);
	        	$tempPDLB{$residueIndex}=pdl(@temp);
	        	push(@{$tempPDLbType{$residueIndex}},@tempB);

			@temp=();
			@tempB=();
			push(@consecResidues,$residue);
			$residueIndex++;
				
		}
		$lastrecord=$record;
		$record = "";
 	}
	print "Done building covalent geometry of all chains.\n\n";
}

# returnFunction: Return the fType and directive field for a specified function
sub returnFunction
{
	my($funcString) = @_;
	my $addExclusions;
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
        my %visitedList;
	my @nextround;
	my $size =keys %union;
	if($size == 0){
		return(-1,0);
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
	$found+= scalar keys %visitedList;

	my $missing=0;
	foreach my $atom(sort {$a <=> $b} keys %union){
		if(!exists $visitedList{$atom}){
			$missing++;
			print "\n!!!Unable to connect the atom at position $atom of chain $chid to the rest of the chain!!!\n";
		}
	}
	return($found,$missing);
}

sub GenerateBondedGeometry {

	my ($connect,$counter,$chid,$chainlength,$lastread) = @_;
	## $connect is a list of connected residues ##
   	my($bH,$angH,$diheH,$map,$bondMapHashRev,$union,$union2,$unionEG) = GenAnglesDihedrals($connect,$chainlength,$chid);
	my %union=%{$union};
	if($chainlength == 0){
		# no atoms in the chain can only happen if we have a TER TER, or TER END.  So, ignore
		return($union2,$unionEG);
	}elsif($chainlength != 1){
		print "Attempting to connect all atoms in chain $chid to the first atom...\n";
		my ($connected,$missed)=connectivityCheck(\%union,$chid);
		if($connected == -1){
			print "\tChain $chid has no bonds. No connections possible. May be a listing of ions.\n\n";
			# this chain has no bonds, so no need to try and connect things
			return($union2,$unionEG);
		}
		if($missed==0 && $connected == $chainlength){
			print "\tAll $connected atoms connected via covalent bonds \n"; 
		}else{
			smog_quit("We appear to have connected $connected of $chainlength atoms in chain $chid.  There is an issue connecting atoms to the rest of the chain using covalent bond definitions.\nThere is probably a missing bond definition in the .bif file.\nCheck for earlier warning messages. Last line read: $lastread ")
		}
	}elsif($chainlength == 1){
		print "Only 1 atom in chain $chid.  Will not perform connectivity checks.\n";
	}

	print "Generating bonds for chain $chid.\n";

    	my @tempArr=();
	## BOND ##
    	for(my $i=0;$i<scalar(@{$bH})-1;$i+=2) {	
		my $bH1=$map->{$bH->[$i]};
		my $bH2=$map->{$bH->[$i+1]};

		my $bondStrA = $bH1->[0];
		my $bondStrB = $bH2->[0];
		my $sizeA = $bH1->[2];
		my $sizeB = $bH2->[2];
		my $ra=$connect->[$bH1->[1]];
		my $rb=$connect->[$bH2->[1]];

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
		my($a,$b,$c) = split("-",$angs);
		my ($na,$mar,$sizeA)=@{$map->{$a}}[0..2];
		my ($nb,$mbr,$sizeB)=@{$map->{$b}}[0..2];
		my ($nc,$mcr,$sizeC)=@{$map->{$c}}[0..2];

		my $ra = $connect->[$mar];
		my $rb = $connect->[$mbr];
		my $rc = $connect->[$mcr];

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
		my($a,$b,$c,$d) = split("-",$dihes);
		##[AtomName,ResidueIndex,prevSize]##

		my ($na,$mar,$sizeA)=@{$map->{$a}}[0..2];
		my ($nb,$mbr,$sizeB)=@{$map->{$b}}[0..2];
		my ($nc,$mcr,$sizeC)=@{$map->{$c}}[0..2];
		my ($nd,$mdr,$sizeD)=@{$map->{$d}}[0..2];

		my $ra = $connect->[$mar];
		my $rb = $connect->[$mbr];
		my $rc = $connect->[$mcr];
		my $rd = $connect->[$mdr];

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

	appendImpropers($map,$connect,$bondMapHashRev,\@tempArr,\%union);

	$DihedralData{$counter} = pdl(@tempArr);
	@tempArr = ();
	print "Done generating bonded geometry of chain $chid.\n\n";
	return($union2,$unionEG);
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


sub checkconnection
{
	my ($i,$c0,$c1,$chid)=@_;
	if(!exists $connections{$residues{$c0}->{"residueType"}}->{$residues{$c1}->{"residueType"}}){
		my $typeA=$residues{$c0}->{"residueType"};
		my $typeB=$residues{$c1}->{"residueType"};
		my $ii=$i+1;
		if((!defined $residues{$c0}->{"connect"} || $residues{$c0}->{"connect"} eq "yes") && (!defined $residues{$c1}->{"connect"} || $residues{$c1}->{"connect"} eq "yes")){
			smog_quit("Connection not defined between residues of type $typeA ($c0) and $typeB ($c1). Check .bif file. Issue encountered when connecting residue $i and $ii in chain $chid (residue index within chain, starting at 1)")
		}elsif(defined $residues{$c0}->{"connect"} && $residues{$c0}->{"connect"} ne "no"){
			my $tmp=$residues{$c0}->{"connect"};
			smog_quit("residue connect value must be yes, or no.  $tmp found.")
		}elsif(defined $residues{$c1}->{"connect"} && $residues{$c1}->{"connect"} ne "no"){
			my $tmp=$residues{$c1}->{"connect"};
			smog_quit("residue connect value must be yes, or no.  $tmp found.")
		}

	}
}

sub connCreateInteractionsSingleBOND
{
    	my($union,$unioneG,$bTypesFragment,$counter,$chain1,$chain2,$index1,$index2,$atomA,$atomB,$resA,$resB,$bEG,$local2global) = @_;
	my @local2global=@{$local2global};
	my $a; my $b; my $c; my $d; my $ta; my $tb; my $tc; my $td;   
	my $if;
	my $eG;
 	## AD-HOC BONDS ##
	my @localbTypes=@{$bTypesFragment};
	
	# add the bond to the fragment and then generate all possible angles and dihedrals within the fragment.
	my($angH,$diheH)=createConnection($index1,$index2,$union,$unioneG,$bEG);
	## Add the BOND to the model##
	$ta = $localbTypes[$index1];
	$tb = $localbTypes[$index2];
	$if = funcToInt("bonds",connWildcardMatchBond($ta,$tb),"");
    	$BondData{$counter}=pdl($index1,$index2,$if);
	excludecontacts($local2global,$index1,$index2);	

#	## Add all ANGLES that contain both of the BONDed atoms##
	my @tempArr;
	foreach my $angs(@{$angH})
	{
		my($a,$b,$c) = split("-",$angs);
		if( ($a == $index1 || $b == $index1  || $c == $index1) && ($a == $index2 || $b == $index2  || $c == $index2) ) {
			$ta = $localbTypes[$a];
			$tb = $localbTypes[$b];
			$tc = $localbTypes[$c];
			$if = funcToInt("angles",connWildcardMatchAngles($ta,$tb,$tc),"");
			push(@tempArr,pdl($a,$b,$c,$if));

			excludecontacts($local2global,$a,$b,$c);	
		}
	}
        if(@tempArr)
        {
		$AngleData{$counter} = cat(@tempArr);
	}else{
		smog_quit("No angles were generated about this user-defined BOND. This is unusual. If would be good to check your structure, BOND definitions and templates, to make sure this is what you expect.",0);
        }
	# array to hold dihedrals that have the BOND as the central bond.  i.e. i-j-k-l and j-k is the bond
	@tempArr = ();
	# arrays for dihedrals to be added to chain 1 and chain 2. i.e. BOND is i-j or l-k
	my  @tempArr1=();		
	my  @tempArr2=();		

	my %eGhash;
	#unroll the energy groups for quick access.
	foreach my $atom(keys %{$unioneG}){
		my @tarra=@{${$union}{$atom}};
		my @tarreG=@{${$unioneG}{$atom}};
		foreach my $atom2(@tarra){
			my $egt = shift @tarreG;
			$eGhash{"$atom-$atom2"}=$egt;
			$eGhash{"$atom2-$atom"}=$egt;
		}
	}

	## DIHEDRALS ##
	foreach my $dihes(@{$diheH})
	{
		my($a,$b,$c,$d) = split("-",$dihes);

		$ta = $localbTypes[$a];
		$tb = $localbTypes[$b];
		$tc = $localbTypes[$c];
		$td = $localbTypes[$d];
		$eG=$eGhash{"$b-$c"};
		$if = funcToInt("dihedrals",connWildcardMatchDihes($ta,$tb,$tc,$td,$eG),$eG);
		$eG = $eGRevTable{$eG};


		if( ($b == $index1  && $c == $index2) || ( $b == $index2  &&  $c == $index1 ) ) {
			# dihedrals formed about the new bond
			push(@tempArr,[$a,$b,$c,$d,$if,1,$eG]);
			excludecontacts($local2global,$a,$b,$c,$d);	
		}
		if( ($a == $index2  && $b == $index1) || ($d == $index2  && $c == $index1)  ) {
			# dihedrals formed about existing bond in chain 1
			push(@tempArr1,[$a,$b,$c,$d,$if,1,$eG]);	
			excludecontacts($local2global,$a,$b,$c,$d);	
		}

		if( ($a == $index1  && $b == $index2) || ($d == $index1  && $c == $index2) ) {
			# dihedrals formed about existing bond in chain 2
			push(@tempArr2,[$a,$b,$c,$d,$if,1,$eG]);	
			excludecontacts($local2global,$a,$b,$c,$d);	
		}


	}

        if(@tempArr)
        {
		$DihedralData{$counter} = pdl(@tempArr);
	}else{
		smog_quit("No dihedral angles were generated about this user-defined BOND. This is unusual. If would be good to check your structure, BOND definitions and templates, to make sure this is what you expect.",0);
        }


	if(@tempArr1){
		addBONDdihToChain($counter,\@tempArr1,$chain1,$index2);
	}
	if(@tempArr2){
		addBONDdihToChain($counter,\@tempArr2,$chain2,$index1);
	}
}

sub excludecontacts 
{
	my @vals=@_;
	my $local2global=shift @vals;
	my @local2global=@{$local2global};
	for(my $I=0;$I<scalar(@vals)-1;$I++){
		for(my $J=$I+1;$J<scalar(@vals);$J++){
			my $a=$vals[$I];
			my $b=$vals[$J];
			$ignoreBONDed{"$local2global[$a] $local2global[$b]"}=1;
			$ignoreBONDed{"$local2global[$b] $local2global[$a]"}=1;
		}
	}
}

sub addBONDdihToChain
{
	my ($counter,$tempArr1,$chain1,$index2)=@_;
	my @tempArr1=@{$tempArr1};
	# find $index2 from piddle for counter and add it to the piddle for chain1 
	# globals: bondPDL, 
	#find PDL entry $bondPDL{$chain1} for chain1
	my $index2t = slice($bondPDL{$counter},":,$index2");
	#glue information about atom $a (i.e. chain2 atom) in $bondPDL{$counter}
	$bondPDL{$chain1}=$bondPDL{$chain1}->glue(1,$index2t);
	#find mapping between current BOND piddle for counter and PDL entry for chain1(column 3)

	my @listcounter=slice($bondPDL{$counter},3);
	my @listchain1=slice($bondPDL{$chain1},3);
	my $nc=0;
	my %mapto1;
	for my $h(list(@listcounter)){
		my $nc1=0;
		for my $hh(list(@listchain1)){
			if($h==$hh){
				$mapto1{$nc}=$nc1;
			}
			$nc1++;
		}
		$nc++;
	}
	#re-assign all indexes used in tempArr1 according to mapping to chain 1
	my @newArr1;
	foreach my $Y(@tempArr1){
		push(@newArr1,[$mapto1{$Y->[0]},$mapto1{$Y->[1]},$mapto1{$Y->[2]},$mapto1{$Y->[3]},$Y->[4],$Y->[5],$Y->[6]]);
	}
	#glue converted tempArr1 to $DihedralData{$chain1}
	$DihedralData{$chain1}=($DihedralData{$chain1})->glue(1,pdl(@newArr1));
}


sub appendImpropers
{
	my($map,$connect,$bondMapHashRev,$tempArr,$union) = @_;
	my $connHandle;
	#loop through the residues in the chain
	# don't loop through the last residue, since the loop always checks the current, plus next
	for(my $resIndA=0;$resIndA<=$#$connect;$resIndA++){

		my $resA=$connect->[$resIndA];
		my $resAIp = $residues{"$resA"}->{"impropers"};

	 	if($resIndA == 0) {
			DoImproperSet($map,$connect,$bondMapHashRev,$tempArr,$union,$resAIp,$resIndA,"no","nohandle");	 
		}

		# if not the final residue, then also check the next residue
		if($resIndA != $#$connect){
			my $resIndB=$resIndA+1;
			my $resB=$connect->[$resIndB];
			my $resBIp = $residues{"$resB"}->{"impropers"};
			DoImproperSet($map,$connect,$bondMapHashRev,$tempArr,$union,$resBIp,$resIndB,"no","nohandle");	 
	   		# check for connection-associated impropers
	   		$connHandle = $connections{$residues{$resA}->{"residueType"}}->{$residues{$resB}->{"residueType"}}->{"improper"};
			DoImproperSet($map,$connect,$bondMapHashRev,$tempArr,$union,$resAIp,$resIndA,$resIndB,$connHandle);	 
	 	}
	}
}

sub DoImproperSet
{
	my ($map,$connect,$bondMapHashRev,$tempArr,$union,$curres,$resInd,$resInd2,$connecthandle)=@_;
	my $imphandle;
	my $isconnection;
	# check if this is a connection, or single residue
	if($resInd2 eq "no"){
		$imphandle=$curres;
		$isconnection="no";
	}else{
		$imphandle=$connecthandle;
		$isconnection="yes";
	}

	foreach my $ips(@{$imphandle})
	{
		if(! (defined $ips) ) {next;}

		my $ia;my $ib;my $ic;my $id;
		my $ta;my $tb;my $tc;my $td;
	        my $na;my $nb;my $nc;my $nd;
		my $ra;my $rb;my $rc;my $rd;
		my $sizeA; my $sizeB;my $sizeC;my $sizeD;
		my ($a,$b,$c,$d);
		if($isconnection eq "no"){
			#just process a single residue
			($a,$b,$c,$d) = @{$ips};
			$a=$bondMapHashRev->{"$a-$resInd"};
			$b=$bondMapHashRev->{"$b-$resInd"};
			$c=$bondMapHashRev->{"$c-$resInd"};
			$d=$bondMapHashRev->{"$d-$resInd"};
		}else{
			# do a connection improper
			($a,$b,$c,$d)=splitconnectionatoms($ips->{"atom"},$bondMapHashRev,$resInd,$resInd2);
		}

		##[AtomName,ResidueIndex,prevSize]##
		$na = $map->{$a}->[0];
		$nb = $map->{$b}->[0];
		$nc = $map->{$c}->[0];
		$nd = $map->{$d}->[0];
		$ra = $connect->[$map->{$a}->[1]];
		$rb = $connect->[$map->{$b}->[1]];
		$rc = $connect->[$map->{$c}->[1]];
		$rd = $connect->[$map->{$d}->[1]];
		$sizeA=$map->{$a}->[2];
		$sizeB=$map->{$b}->[2];
		$sizeC=$map->{$c}->[2];
		$sizeD=$map->{$d}->[2];
		($ia,$ta) = ($sizeA+getAtomIndexInResidue($ra,$na),getAtomBType($ra,$na));
		($ib,$tb) = ($sizeB+getAtomIndexInResidue($rb,$nb),getAtomBType($rb,$nb));
		($ic,$tc) = ($sizeC+getAtomIndexInResidue($rc,$nc),getAtomBType($rc,$nc));
		($id,$td) = ($sizeD+getAtomIndexInResidue($rd,$nd),getAtomBType($rd,$nd));	

		if($isconnection ne "no"){
			# we only check for connections, since impropers within a 
			# dihedral are checked during template validation
			checkimproperdefined($union,$a,$b,$c,$d,$na,$nb,$nc,$nd,$resInd,$resInd2);
		}

		($nb,$nc) =  ($map->{$b}->[1]-$map->{$c}->[1]==0)?($nb,$nc):("nb?",$nc);
		my $if = funcToInt("impropers",connWildcardMatchImpropers($ta,$tb,$tc,$td),"");	
		push(@{$tempArr},[$ia,$ib,$ic,$id,$if,1,-1]);	
	}
}

sub checkimproperdefined
{
	my ($union,$a,$b,$c,$d,$na,$nb,$nc,$nd,$resInd,$resInd2)=@_;
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

 	if($IMPFLAG2==0){
		my $rr=$resInd+1;
		my $rr2=$resInd2+1;
		smog_quit("There is an incorrectly formed improper dihedral defined between residues $rr and $rr2 (index starting at 1). Three atoms must be bonded to a central atom. Improper defined by atoms $na-$nb-$nc-$nd.\nThere may be a missing bond, or incorrectly defined improper in the .bif file.\n");
	}
}

sub splitconnectionatoms
{
	# this determines if the atoms in a connection improper are in the first, or second residue
	my ($atomshandle,$bondMapHashRev,$resIndA,$resIndB)=@_;
	my @atoms = @{$atomshandle};

	for(my $I=0;$I<4;$I++){
		my $atomname;
		my $a=$atoms[$I];
		if( $a =~ s/\+$//g ){
			$atomname=$a;
			$a=$bondMapHashRev->{"$a-$resIndB"};
			if(!defined $a){smog_quit("Connection improper dihedral has unknown atom: $atomname not found in definition for residue $resIndB");}
		}else{
			$atomname=$a;
			$a=$bondMapHashRev->{"$a-$resIndA"};
			if(!defined $a){smog_quit("Connection improper dihedral has unknown atom: $atomname not found in definition for residue $resIndA");}
		}
		$atoms[$I]=$a;
	}
	return @atoms;
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
		unless(("$a" ne "$aM" && $aM ne "*")
			|| ("$b" ne "$bM" && "$bM" ne "*")
			|| ("$c" ne "$cM" && "$cM" ne "*")){
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
		smog_quit ("Multiple possible angles match $a-$b-$c equally well (score=$saveScore;$matchScoreCount). Unclear assignment of function type");
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
		# exact match
		$funct = $interactions->{"bonds"}->{$typeA}->{$typeB};
	}elsif(
           $typeA ne $typeB && (exists $interactions->{"bonds"}->{$typeA}->{"*"} 
           && exists $interactions->{"bonds"}->{$typeB}->{"*"})){
		# two partial WC matches
		smog_quit ("Wildcard conflict in bonds $typeA-$typeB. Both $typeA-\* and $typeB-\* are defined in .b file. Can not unambiguously assign a function...");
	}elsif(exists $interactions->{"bonds"}->{$typeA}->{"*"}){
	## If typeA exists while TypeB is a wildcard ##
 		$funct = $interactions->{"bonds"}->{$typeA}->{"*"};
	}elsif (exists $interactions->{"bonds"}->{$typeB}->{"*"}){
		## If typeB exists while TypeA is a wildcard ##
		$funct = $interactions->{"bonds"}->{$typeB}->{"*"};
	}
	if(!defined $funct || $funct eq ""){
		if(exists $interactions->{"bonds"}->{"*"}->{"*"}){
			$funct = $interactions->{"bonds"}->{"*"}->{"*"};
		}else{
			smog_quit ("Unable to assign bond types to all bonds in a residue\n Offending btypes are $typeA $typeB. Check .b file");
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
		if(        ("$a" ne "$aM" && "$aM" ne "*")
			|| ("$b" ne "$bM" && "$bM" ne "*")
			|| ("$c" ne "$cM" && "$cM" ne "*")
			|| ("$d" ne "$dM" && "$dM" ne "*")){next;}
		if("$a" eq "$aM") {$matchScore+=2;} else {$matchScore+=1;}
		if("$b" eq "$bM") {$matchScore+=2;} else {$matchScore+=1;}
		if("$c" eq "$cM") {$matchScore+=2;} else {$matchScore+=1;}
		if("$d" eq "$dM") {$matchScore+=2;} else {$matchScore+=1;}
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
	my $NumOfKeys=keys %{$diheHandle};
	my @keys = keys %{$diheHandle};
	if($NumOfKeys == 1 && $keys[0] eq "*-*-*-*"){
		$funct = $diheHandle->{$keys[0]};

	}else{
		my $matchScore;
		my $saveScore = 0;
		my $matchScoreCount=0; 
		my $symmatch=0; 
		my $Nd=0;
		foreach my $matches(keys %{$diheHandle})
		{
			$Nd++;
			$matchScore = 0; 
			# this step can be done once, rather than for each call.
			my ($aM,$bM,$cM,$dM) = split("-",$matches);
			if($matches eq "*-*-*-*"){
				$matchScore=4;
			}else{
	
				if(        ("$a" ne "$aM" && "$aM" ne "*")
					|| ("$b" ne "$bM" && "$bM" ne "*")
					|| ("$c" ne "$cM" && "$cM" ne "*")
					|| ("$d" ne "$dM" && "$dM" ne "*")){next;}
	
				if("$a" eq "$aM") {$matchScore+=2;} else {$matchScore+=1;}
				if("$b" eq "$bM") {$matchScore+=2;} else {$matchScore+=1;}
				if("$c" eq "$cM") {$matchScore+=2;} else {$matchScore+=1;}
				if("$d" eq "$dM") {$matchScore+=2;} else {$matchScore+=1;}
	
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
	
			smog_quit ("symm $symmatch sym $sym msc $matchScoreCount Multiple possible angles match $a-$b-$c-$d, and energyGroup $eG equally well. Can not determine function based on .b file.");
		}
	
		if($saveScore == 0){
			smog_quit ("Dihedral Angle between bTypes $a-$b-$c-$d and energyGroup $eG: Unable to match to a function in .b file.");
		}
	
		if(!defined $funct || $funct eq ""){smog_quit("There is no function for bType combination $a-$b-$c-$d with energyGroup=$eG. Check .b file");}
	}

	return $funct;
}

sub GenAnglesDihedrals
{
	my($connect,$chainlength,$chid) = @_;
	## $connect is a list of connected residues ##
	my $dihes; my $angles; 
	my $oneFour;
	my %union; 
        my %union2; # this one is specifically numbered to make BOND additions. In earlier implementations, we were limited to bonds between atoms not involved in connections.  Now, we are keeping the full connectivity in a form that is easier to work with, so that arbitrary blocks of atoms may be considered during BOND dihedral identification steps
	my %unionEG; # this carries the energy groups of each bond in the union2 lists
        my $connHandle;
	my $i=0;my $mapCounter=0;
	my %bondMapHash; ##[AtomName,ResidueIndex,prevSize]##
	my %bondMapHashRev;
	my @connectList;
	my @connectListRes;
	my $leftAtom;my $rightAtom;
	my $prevSize = 0;
	
	## Go through list of connected residues ##
	for($i = 0;$i<=$#$connect;$i++)
	{
  		my $resABonds = $dihedralAdjList{$connect->[$i]};
  		my $resABondseG = $dihedralAdjListeG{$connect->[$i]};
		my $resAAtoms = $residues{$connect->[$i]}->{"atoms"};
		## Atoms to mapCounter renaming ##
 		foreach my $atom( keys %{$resAAtoms})
		{
			$bondMapHashRev{"$atom-$i"}=$mapCounter;
			$bondMapHash{$mapCounter}=[$atom,$i,$prevSize];
			$mapCounter++;  
		}
		foreach my $atom(keys %{$resABonds})
    		{
     			my @tempArr = map {$bondMapHashRev{"$_-$i"}} @{$resABonds->{$atom}};
     			my @tempArreG = @{$resABondseG->{$atom}};
			my $atomKey = $bondMapHashRev{"$atom-$i"};


     			$union{$atomKey} = \@tempArr;

			my $atomname=$bondMapHash{$atomKey}->[0];
			my $residue=$connect->[$i];

			my $ia = $prevSize+getAtomIndexInResidue($residue,$atomname);
			my @ntempArr=();
			foreach my $JJ(@tempArr){
				$atomname=$bondMapHash{$JJ}->[0];
				my $ja = $prevSize+getAtomIndexInResidue($residue,$atomname);
				push(@ntempArr,$ja);
			}
     			$union2{$ia} = \@ntempArr;
     			$unionEG{$ia} = \@tempArreG;

		}

		# if this is a single residue chain, then don't try to connect it to the next residue
        	if($#$connect == 0) {last;}
		## Start of chain no inter residue connection ##
    		#  but setup leftAtom and leftResidue sizes #
    		if($i == 0) 
		{
			checkconnection($i+1,$connect->[0],$connect->[1],$chid);
			if(exists $connections{$residues{$connect->[0]}->{"residueType"}}->{$residues{$connect->[1]}->{"residueType"}}){
				$connHandle = $connections{$residues{$connect->[0]}->{"residueType"}}->{$residues{$connect->[1]}->{"residueType"}};
	  			$leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
				if(! defined $bondMapHashRev{"$leftAtom-$i"}){
					smog_quit("Unable to connect residues 0 and 1.  It appears that the connecting atom in residue 0 is missing. Connection expects there to be a $leftAtom atom.");
				}
	  			$leftAtom = $bondMapHashRev{"$leftAtom-$i"};
	 			$prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});		
			}
      			next;
		}
		checkconnection($i,$connect->[$i-1],$connect->[$i],$chid);
        	## $i > 0, create inter residue connection ##
		## $i-1 <--> $i
		if(exists $connections{$residues{$connect->[$i-1]}->{"residueType"}}->{$residues{$connect->[$i]}->{"residueType"}}){
			$connHandle = $connections{$residues{$connect->[$i-1]}->{"residueType"}}->{$residues{$connect->[$i]}->{"residueType"}};
			$rightAtom = $connHandle->{"bond"}->[0]->{"atom"}->[1];
			if(! defined $bondMapHashRev{"$rightAtom-$i"}){
				my $j=$i-1;
				smog_quit("Unable to connect residues $j and $i.  It appears that the connecting atom in residue $i is missing. Connection expects there to be a $rightAtom atom.");
			}
	    		$rightAtom = $bondMapHashRev{"$rightAtom-$i"};
			push(@connectList,$leftAtom);
			push(@connectList,$rightAtom);
			push(@connectListRes,$connect->[$i-1]);
			push(@connectListRes,$connect->[$i]);
		}
    
    		## $i <--> $i+1
        	if($i == $#$connect) {last;}
		checkconnection($i+1,$connect->[$i],$connect->[$i+1],$chid);
		if(exists $connections{$residues{$connect->[$i]}->{"residueType"}}->{$residues{$connect->[$i+1]}->{"residueType"}}){
        		$connHandle = $connections{$residues{$connect->[$i]}->{"residueType"}}->{$residues{$connect->[$i+1]}->{"residueType"}};
    			$leftAtom = $connHandle->{"bond"}->[0]->{"atom"}->[0];
    			$leftAtom = $bondMapHashRev{"$leftAtom-$i"};
    			$prevSize = $prevSize+scalar(keys %{$residues{$connect->[$i]}->{"atoms"}});
		}
   	}
  	## Create Inter residue connection ##
  	for($i=0;$i<scalar(@connectList)-1;$i+=2) {
   		push(@{$union{$connectList[$i]}},$connectList[$i+1]);
   		push(@{$union{$connectList[$i+1]}},$connectList[$i]);

		my @ntempArr=();
		my $JJ=$connectList[$i];
		my $residue=$connectListRes[$i];
		my $atomname=$bondMapHash{$JJ}->[0];
		my $prevSize=$bondMapHash{$JJ}->[2];
		$JJ = $prevSize+getAtomIndexInResidue($residue,$atomname);
		my $KK=$connectList[$i+1];
		my $prevSize1=$bondMapHash{$KK}->[2];
		my $atomname1=$bondMapHash{$KK}->[0];
		my $residue1=$connectListRes[$i+1];
		$KK = $prevSize1+getAtomIndexInResidue($residue1,$atomname1);
		my $eG = getEnergyGroup($residue,$residue1,$atomname,"$atomname1?");
   		push(@{$unionEG{$JJ}},$eG);
   		push(@{$unionEG{$KK}},$eG);
   		push(@{$union2{$JJ}},$KK);
   		push(@{$union2{$KK}},$JJ);

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
			if(defined $residues{$residue}->{"atoms"}->{$atomname}->{"bond"} && $residues{$residue}->{"atoms"}->{$atomname}->{"pairType"} !=0){
        			smog_quit("No bonds found with atom $atomname in residue $residue. Check .bif definitions.");
			}
   		}
	}

  	($dihes,$angles,$oneFour)=adjListTraversal(\%union);
  	return (\@connectList,$angles,$dihes,\%bondMapHash,\%bondMapHashRev,\%union,\%union2,\%unionEG);
}

sub createConnection
{
	my($atom1,$atom2,$union,$unioneG,$eG) = @_;
	my @tempA; my @tempD;
	my $dihes; my $angles; my $oneFour;
	my %union=%{$union}; 
	my %unioneG=%{$unioneG}; 
	my $resABonds;my $resBBonds;
	my %tempAdjList; 

  # simple addition of the desired bond
	push(@{$union{"$atom1"}},"$atom2");
	push(@{$union{"$atom2"}},"$atom1");
	push(@{$unioneG{"$atom1"}},"$eG");
	push(@{$unioneG{"$atom2"}},"$eG");
	 
	($dihes,$angles,$oneFour)=adjListTraversal(\%union); 
  	 
	return ($angles,$dihes);
}

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
	my($fileName,$fileName2,$userProvidedMap,$CGenabled,$saveSCMorig) = @_;
	my $numContacts = 0; 
	my $line = "";
	my $chain1;my $chain2; my $contact1; my $contact2; my $res1; my $res2;
	my $dist;
	my $x1;my $x2;my $y1;my $y2;my $z1;my $z2;
	my %resContactHash; my $skip = 1; my $COARSECONT; 
	my @interiorTempPDL; #usage: push(@interiorTempPDL,[1,$contact1,$contact2,$dist]);
	#Format for this PDL has a boolean as the first argument
	#it is unused for now, but could be useful in future to use
	#as a flag to differentiate between user generated and smog generated contacts

	if(!$userProvidedMap){ #use shadow-generated contact map
		## OPEN .contact FILE ##
		unless (open(CONTFILE, $fileName)) {
			smog_quit ("Unable to read shadow-generated contact map. This typically means there was an issue running shadow.");
		}
		my $coarseFile = $fileName.".CG";
		if($CGenabled == 1) { 
			unless (open($COARSECONT,">$coarseFile")) {
				smog_quit ("Internal contact file cannot be written.");
			} 
		}
		while($line = <CONTFILE>)
		{
			($chain1,$contact1,$chain2,$contact2) = split(/\s+/,$line);
			if(whatAmI($chain1) !=1 || whatAmI($contact1) !=1 ||whatAmI($chain2) !=1 ||whatAmI($contact2) !=1 ){
				smog_quit("Internal error: Somehow, a non-integer value was found for a chain, or atom, in the contact file.\n LINE: $line");
			}
			if($CGenabled == 1) {
				#if we are coarse graining then we need to map AA to residue contacts
				#first see which residues they belong
				$res1 = $allAtomsBackup{$contact1}[9]; $res2 = $allAtomsBackup{$contact2}[9];
				$skip = 0; 	#skip if this res1-res2 has already been added
				if(!exists $resContactHash{"$res1,$res2"}) {
					$resContactHash{"$res1,$res2"} = 1;
					$skip = 1;
					$contact1 = $res1; $contact2 = $res2;
					print $COARSECONT "$chain1 $res1 $chain2 $res2\n";
				}
			}
			if($skip == 1) { #maybe we skip sometimes if coarse graining				
				($x1,$y1,$z1) = @{$allAtoms{$contact1}}[6..8];
				($x2,$y2,$z2) = @{$allAtoms{$contact2}}[6..8];

				$dist = sqrt( ($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2) * $angToNano;
				my $mindist=$interactionThreshold->{"contacts"}->{"shortContacts"};
				if($dist < $mindist)
				{
					if($main::setContacttoLimit){
	                        		smog_note("Contact between atoms $contact1 $contact2 below threshold.\n-limitcontactlength is being used, will set distance of contact to $mindist");
						$dist=$mindist;
					}elsif($main::DeleteShortContact){
	                        		smog_note("Contact between atoms $contact1 $contact2 below threshold.\n-deleteshortcontact is being used, will exclude this contact.");
						next;
					}else{
	                            		smog_quit("Contact between atoms $contact1 $contact2 below threshold distance with value ". sprintf("%.3f",$dist));
	 			  	}
				}

				if(exists $ignoreBONDed{"$contact1 $contact2"}){
					smog_note("Calculated contact will be excluded because the atoms are associated with a bond, angle or dihedral defined by a \"BOND\". Contact is between atoms $contact1 and $contact2");
					next;
				}
				push(@interiorTempPDL,[$userProvidedMap,$contact1,$contact2,$dist]);

				$numContacts++;
			}
		}
		$contactPDL = pdl(@interiorTempPDL);
		if($CGenabled == 1) { close($COARSECONT); }
		close(CONTFILE);
		#here we will delete (or rename) the shadow output and make a new output contact file that is consistent with the input pdb. What does this mean?
		#It means that instead of using a serial numbering starting from 1 (which is the shadow output), it should have chain and atom number the same as in
		#the PDB. The conversion between serial and PDB is contained inside the allAtoms array which is indexed by serial and has the pdb number as the 11th column.
		#$allAtoms{$atomSerial}=[$nbType,$residueType,$residueIndex,$atom,$chainNumber,$residue,$x,$y,$z,$residueSerial,$pairType,$pdbIndex];
		#unlink($fileName);
		if ($saveSCMorig==1){
			`mv $fileName $fileName.ShadowOutput`;
			unlink($fileName);
		}else{
			`rm $fileName`;
		}
		
		unless (open(CONTFILE,">$fileName")) { smog_quit ("PDB-consistent contact file $fileName cannot be written."); }
		for(my $i=0;$i<$numContacts;$i++)
		{
			my $atoma = sclr(slice($contactPDL,"1:1,$i:$i"));
			my $atomb = sclr(slice($contactPDL,"2:2,$i:$i"));
			my $pdbnuma = $allAtoms{$atoma}[11];
			my $pdbnumb = $allAtoms{$atomb}[11];
			my $chaina = $allAtoms{$atoma}[4]+1;
			my $chainb = $allAtoms{$atomb}[4]+1;		
			print CONTFILE "$chaina $pdbnuma $chainb $pdbnumb\n";
		}
		close(CONTFILE);
	} else { #read in contact from file
		if($CGenabled == 1) {
			smog_quit ("User input contact map is not supported with automatic coarse graining (e.g., -CA, -CAgauss, -tCG).".
			" This is because the coarse graining options are only to enable automatic contact generation from an all-atom input pdb.".
			" If the user wishes to use their own contact map, use a .bif and .pdb that already have the correct graining.".
			" For example, for a C-alpha model, directly use the template \$SMOG2_LOCATION/templates/SBM_calpha with the -t option." ); }
		print "Not going to calculate a contact map\n";
		print "Will read contacts from $fileName2\n";
		## OPEN user provided contact FILE ##
		unless (open(CONTFILE1, $fileName2)) {
			smog_quit ("Cannot read contact file '$fileName2'.");
		}
		## User contact map should be in format below and use input PDB numbering ##
		## chain1 atom1 chain2 atom2 ##
		while($line = <CONTFILE1>) {
			my ($data,$comment)=checkcomment($line);
			if($data eq ""){next;}
			my @A = split(/\s+/,$data);
			if($#A < 3 || $#A>4){smog_quit("incorrect number of entries in contact file.  See line: $line");}
			my ($chain1,$pdbNum1,$chain2,$pdbNum2,$dist) = split(/\s+/,$line);
			if(whatAmI($chain1) !=1 || whatAmI($pdbNum1) !=1 ||whatAmI($chain2) !=1 ||whatAmI($pdbNum2) !=1 ){
				smog_quit("non-integer value given for chain, or atom, in contact file.")
			}
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

				($x1,$y1,$z1) = @{$allAtoms{$contact1}}[6..8];
				($x2,$y2,$z2) = @{$allAtoms{$contact2}}[6..8];
				$dist = sqrt( ($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2);
			}
			$dist = $dist * $angToNano;
			if(!exists $allAtoms{$contact1}){smog_quit("ATOM $contact1 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}
			if(!exists $allAtoms{$contact2}){smog_quit("ATOM $contact2 doesn't exists. Skipping contacts $contact1-$contact2\n");next;}

			my $mindist=$interactionThreshold->{"contacts"}->{"shortContacts"};
			if($dist < $mindist)
			{
				if($main::setContacttoLimit){
                        		smog_note("Contact between atoms $contact1 $contact2 below threshold.\n-limitcontactlength is being used, will set distance of contact to $mindist");
					$dist=$mindist;
				}elsif($main::DeleteShortContact){
                        		smog_note("Contact between atoms $contact1 $contact2 below threshold.\n-deleteshortcontact is being used, will exclude this contact.");
					next;
				}else{
                            		smog_quit("Contact between atoms $contact1 $contact2 below threshold distance with value $dist");
 			  	}
			}
			if(exists $ignoreBONDed{"$contact1 $contact2"}){
				smog_note("Listed contact will be excluded because the atoms are associated with a bond, angle or dihedral defined by a \"BOND\". Contact is between atoms $contact1 and $contact2");
				next;
			}

			push(@interiorTempPDL,[$userProvidedMap,$contact1,$contact2,$dist]);
			$numContacts++;
		}
		$contactPDL = pdl(@interiorTempPDL);
	} 
		
	#exit(0);
	return $numContacts;
  
}

1;
