##
## HELPS TO REALIGN RESIDUE NUMBERS AND ATOM INCIDES IN PDB
##
use strict;
use warnings;

my $inputPDB = $ARGV[0];
my $outputPDB = $ARGV[1];
open(FILE,"$inputPDB");open(OUTPUT,">$outputPDB");
my $atomCounter = 1; my $residueCounter = 1;
while(my $line = <FILE>)
{
	if($line =~ m/END/){last;}
	if($line =~ m/TER/){next;}
	my @splits = split(/\s+/,$line);
	$splits[1]=$atomCounter;
	my $atom = "ATOM";
	if($line =~ m/HETATM/){$atom="HETATM";}
 $line = sprintf "$atom  %5d %4s %-3s  %4d    %8.3f%8.3f%8.3f\n",$splits[1],$splits[2],$splits[3],$splits[5],$splits[6],$splits[7],$splits[8];
										##$serial_n,$atom_name, $res_name, $seq_n, $ins_code, @coords
	$atomCounter++;	
	print $line;
 #$splits[1]=$atomCounter;print @splits,"\n";

	


}
