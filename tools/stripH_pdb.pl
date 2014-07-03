use strict;
use warnings;

open(FILE,$ARGV[0]); ## INPUT FILE
open(OUTPUT,">ignH_$ARGV[0]"); ## OUTPUT FILE

my %excludeAtoms; my %indexMap;
my $counter = 1; my $resMod = 1;
while(my $line = <FILE>)
{
	
   	do
				{
						if($line !~ m/ATOM/) { next;}
    		my @attrs = split(/\s+/,$line);
    		my $atom = $attrs[2];
						my $index = $attrs[1];
						if($counter == 1) {$resMod = $attrs[5];}
						my $resnum = $attrs[5] % $resMod + 1;
      if($atom =~ m/^H/) 
						{
								#print "Removing $atom,$attrs[1]\n";
								$excludeAtoms{$index} = $atom;
						}
						else
						{
								$attrs[1] = $counter;
								$attrs[5] = $resnum;
								$indexMap{$index}=$counter;
								my $sendString = join("\t",@attrs);
							 print $sendString,"\n";
							
								printf OUTPUT "ATOM  %5d  %-4s%3s %s%4d    %8.3f%8.3f%8.3f\n",
 $attrs[1],$attrs[2],$attrs[3],$attrs[4],$attrs[5],$attrs[6],$attrs[7],$attrs[8];
								$counter++;						
						}
				} while(($line = <FILE>) !~ m/END/);
		}
