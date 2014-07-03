use strict;
use warnings;
use Carp;

##
# STRIP HYDROGENS OFF PDB FILE 
# usage: ./stripH.pl <pdbFile> 
# WARNING: THIS PROGRAM DOES NOT GUARANTEE CORECTNESS
# WARNING: THIS PROGRAM ONLY ASSUMES THERE ARE ATOM/TER/END LINES

my $pdb = $ARGV[0];
open(PDB,$pdb);
open(PDBOUT,">$pdb.new");

my $atomCount = 0;
while(my $line = <PDB>)
{
   if($line =~ m/TER/){print PDBOUT $line;next;}
   if($line =~ m/END/){print PDBOUT "END";last;}
   if($line !~ m/^ATOM/){confess("$pdb is not structured properly\n");}

   my $atom = substr($line,12,4);
   $atom =~ s/^\s+//; #remove leading spaces
   $atom =~ s/\s+$//; #remove trailing spacesi
   if($atom =~ m/^H/) {next;}
   $atomCount++;
   my $putString = sprintf('%5s',$atomCount);
   substr($line,6,5,$putString);
   print PDBOUT $line;




}







