#!/usr/bin/perl -w

# USAGE ./adjustPDB <pdb> <map file> <output file>
use strict;
use warnings;
my %map;

## VARS: ##
my $rows;
my $atomCount;
## TAIL/HEAD FLAG INITIALIZATION ##
my $isHead;
my $isTail;
my $resInd;
my $chainName;
my $resName;
my @residue;
my @pArr;
my $loopInd;
my $prevResTypeRNA;

## CHECK INPUT ##
if(scalar(@ARGV) < 3){
	confess("Usage: ./adjustPDB <PDB> <-f filename | sbmMap> <output file>\n 
For default templates use the option sbmMap, else provide a filename with the option -f \n");
}
if(!exists $ENV{"SMOG_PATH"}){
	confess("Environmental Variable SMOG_PATH not set\n");
}
my $mapOption = $ARGV[1];
if($mapOption =~ m/^sbmMap$/){
	$mapOption = "$ENV{SMOG_PATH}/tools/sbmMap";
	open(NEWPDB,">$ARGV[2]") or die "Can't open the all-atom region atoms PDB file. \n";
}
elsif($mapOption =~ m/^-f$/ && exists $ARGV[2]){
	$mapOption = $ARGV[2];
	open(NEWPDB,">$ARGV[3]") or die "Can't open the all-atom region atoms PDB file. \n";
}
else{
	confess("Usage: ./adjustPDB <PDB> <-f filename | sbmMap> <output file>\n 
			For default templates use the option sbmMap, else provide a filename with the option -f \n");
}

open(OLDPDB,"$ARGV[0]") or die "Can't open molecule PDB file.\n";

## CACHE MAP ##
#Note: Map files has to be in format <residue> <head_name> <tail_name>
open(MAP,$mapOption) || confess("Cannot open Map file\n");
while(my $line = <MAP>)
{
  ## SKIP COMMENT LINE ##
  chomp($line);
  if($line =~ /#/){next;}
  my @entries = split(/\s+/,$line);
  #map-->residue-->{head|tail} 
  $map{$entries[0]}{"head"} = $entries[1];
  $map{$entries[0]}{"tail"} = $entries[2];
}

## Returns 1 if a residue is a RNA or DNA, else returns zero. ##
sub isRNA{
	my($resName) = @_;
	my %RNAlist=("U",1,"U5",1,"G",1,"G5",1,"A",1,"A5",1,"C",1,"C5",1,"MIA",1);
	if ($RNAlist{$resName}){
			return 1;
	}
	return 0;
}


sub adjustInputFile(){
	
	# Read the molecule PDB into an array (for convinience)
	
	$rows=0;
	while(<OLDPDB>){
		my $line=$_;
		chomp($line);
		$pArr[$rows]=$line;	
		$rows++;
	}
	close OLDPDB;
	
	$isHead = 1;
	$isTail = 0;
	$resInd = 1;
	$atomCount=1;	
	my $k = 0;
	
	## Loops through all PDB rows ##
	while($k<$rows){
		
		## If end of file: ##
		if ($pArr[$k] =~ /^END/){
			print NEWPDB "END\n";
			last;
		}
		
		## If TER line: ##
		if($pArr[$k] =~ /^TER/){
			print NEWPDB "TER\n";
			$atomCount=1;
			$isHead = 1;
			$resInd = 1;
			$k++;
		}
		
		## If atom line: ##
		if ($pArr[$k] =~ /^ATOM/){

			my $resNum = substr($pArr[$k],22,5);
			$resNum =~ s/^\s+|\s+$//g;		
			my $newResNum = $resNum;
			
			
			$resName = substr($pArr[$k],17,4);
			$resName =~ s/^\s+|\s+$//g;
			
			$chainName = substr($pArr[$k],21,1);
			
			## Loop through the resdiue ##
			$loopInd = 0;
			while ($newResNum eq $resNum){
				
				## Obtain atom information ##
				my $atomName = substr($pArr[$k],12,4);
				$atomName =~ s/^\s+|\s+$//g;
				
				$residue[$loopInd]->{"atomName"} = $atomName;
				$residue[$loopInd]->{"atomIndex"} = $atomCount;
				$residue[$loopInd]->{"x"} = substr($pArr[$k],30,8);   
				$residue[$loopInd]->{"y"} = substr($pArr[$k],38,8);  
				$residue[$loopInd]->{"z"} = substr($pArr[$k],46,8);  
				
				$k++;
				
				## Check if next line is END or TER ##
				if ($pArr[$k] =~ /^END/ || $pArr[$k] =~ /^TER/){
					$newResNum = "";
					$isTail = 1;
				}	
				else{
					$newResNum = substr($pArr[$k],22,5);	#get next residue index
					$newResNum =~ s/^\s+|\s+$//g;		
					if ($newResNum eq $resNum){
						$loopInd++;
						$atomCount++;
					}
				}  	
			}
			
			## Adjust Tail\Head names ##
			if ($isHead){
				$resName = $map{$resName}{"head"};
			}
			if ($isTail){
				if (not($prevResTypeRNA)){	#check that it's not an amino-accilated tRNA
					$resName = $map{$resName}{"tail"};	
				}
			}
			
			for(my $i=0; $i<$loopInd+1; $i++){
					my $aName = $residue[$i]->{"atomName"};
					if ($aName =~ /'/){
						chop($aName);
						$aName .= "*";
					}
					if ($aName eq "OP1"){$aName="O1P"};
					if ($aName eq "OP2"){$aName="O2P"};
					if ($aName eq "OP3"){$aName="O3P"};
					if ($aName eq "CG"){$aName="CG1"};
					my $ind = $residue[$i]->{"atomIndex"};
					my $x = $residue[$i]->{"x"};
					my $y = $residue[$i]->{"y"};
					my $z = $residue[$i]->{"z"};

				    printf NEWPDB "ATOM  %5d %4s %-4s%s%4d    %8.3f%8.3f%8.3f\n",$ind,$aName,$resName,$chainName,$resInd,$x,$y,$z;
			}

			#LEAVE CURRENT RESIDUE AND UPDATE FLAGS#
			$isHead = 0;
			$isTail = 0;
			$resInd++;
			$prevResTypeRNA = isRNA($resName);
			undef @residue; 
		}
	}
	close NEWPDB;
}

adjustInputFile();
