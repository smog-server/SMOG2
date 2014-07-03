#!/usr/bin/perl
## ADJUST PDB TO MODIFY HEAD AND TAIL RESIDUES
# THAT CONFORMS WITH BIF DECLARATION OF RESIDUES
# USAGE ./adjustPDB <pdb> <map file>
use strict;
use warnings;
use Getopt::Long;
use Carp;

## GLOBAL VARIABLE ##
my %map;

## CHECK INPUT ##
if(scalar(@ARGV) < 2) 
{confess("Usage: ./adjustPDB <PDB> <-f filename | sbmMap>\n 
For default templates use the option sbmMap, else provide a filename with the option -f \n");}
if(!exists $ENV{"SMOG_PATH"})
{confess("Environmental Variable SMOG_PATH not set\n");}

## OUTPUT PDB FILE NAME ##
open(OUTPDB,">$ARGV[0].mpd") || confess("Cannot create output PDB\n");

my $mapOption = $ARGV[1];
if($mapOption =~ m/^sbmMap$/)
{$mapOption = "$ENV{SMOG_PATH}/tools/sbmMap";}
elsif($mapOption =~ m/^-f$/ &&
exists $ARGV[2])
{$mapOption = $ARGV[2];}
else{confess("Usage: ./adjustPDB <PDB> <-f filename | sbmMap>\n 
For default templates use the option sbmMap, else provide a filename with the option -f \n");}

## INPUT PDB FILE ##
open(PDB,"$ARGV[0]") || confess("Cannot open input PDB\n");

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

## EDIT PDB FILE ##
my $line = "";
my $resName = "";
my $resIdx = "";
my $outSTR = "";
my $oldIdx = "";my $oldRes = "";
my $startChain = 0;
while($line = <PDB>)
{
  chomp($line);

  ## NEW CHAIN MODIFY TERMINAL ##
  if($line =~ /TER/ || $line =~ /END/)
  {
   if($oldRes eq "") {confess("Your PDB is not correct\n");}

   my $tail = $map{$oldRes}{"tail"};
   if(length($tail) > length($oldRes))
    {$outSTR =~ s/$oldRes\s/$tail/g;}
   else
    {$outSTR =~ s/$oldRes/$tail/g;}
   print OUTPDB $outSTR;
   $oldIdx = "";
   $resIdx = "";$startChain = 1;
   if($line =~ /^END$/){print OUTPDB "END"; last;}
   else{print OUTPDB "TER\n";}
   next;
  }

  ## SKIP NON ATOM LINES ##
  if($line !~ /ATOM/) {print OUTPDB "$line\n";next;}
 
  $oldRes = $resName;
  $resName = substr $line,17,3;
  $oldIdx = $resIdx;
  $resIdx = substr $line,22,5;
  $resName =~ s/(^\s+|\s+$)//g;
  $resIdx =~ s/(^\s+|\s+$)//g;
  
  ## START ##
  if($oldIdx eq "")
  {
   $outSTR = "$line\n";
   $startChain = 1;
   $oldIdx = $resIdx;
   next;
  }
  ## MODIFY FIRST RESIDUE ##
  if($oldIdx ne $resIdx && $startChain)
  {
   $startChain = 0;
   my $head = $map{$oldRes}{"head"};
   
   if(length($oldRes) < length($head))
   {$outSTR =~ s/$oldRes\s/$head/g;}
   else
   {$outSTR =~ s/$oldRes/$head/g;}

   print OUTPDB $outSTR;
   $outSTR = "$line\n";
   next; 
  } 
  elsif($oldIdx ne $resIdx)
  {
    $outSTR = "$outSTR";
    print OUTPDB $outSTR;
    $outSTR = "";
  }  
 
  $outSTR = "$outSTR$line\n";
 
}
print "Your PDB file $ARGV[0] was adjusted using the map $mapOption. The file $ARGV[0].mpd was generated\n";




