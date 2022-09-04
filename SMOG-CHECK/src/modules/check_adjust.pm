package check_adjust;
use strict;
use Exporter;
use smog_common;
use check_common;
our @ISA = 'Exporter';
our @EXPORT = qw(check_adjust);

sub check_adjust
{
 my ($exec,$smogexec,$sharedir)=@_;
 my $NFAIL=0;
 my $MESSAGE="";
 my %FAIL;
 my $FAILED;
 my $FAILSUM=0;
 my $UNINIT;
 my $LINESorig=0;
 my $pdbdir="$sharedir/PDB.files";
 my $mapdir="$sharedir/mapfiles";
 my $origpdb="$pdbdir/3PTA.preadjust.pdb";
 my $newpdb="testname.pdb";
 my $tool="adjust";
 my $TESTNUM=0;
 my $gitshift=getgitver;
 if($gitshift ne ""){
  # we are using a git version, so files will be 2-lines longer
  $gitshift=2;
 }else{
  $gitshift=0;
 }

 open(ORIG,"$origpdb") or internal_error("Unable to open $origpdb");
 while(<ORIG>){
  $LINESorig++;
 }
 my @FAILLIST = ('NON-ZERO EXIT','OUTPUT NAME','FILE LENGTH','SMOG RUNS','LARGE');

 # TEST 1
 print "\tChecking smog_adjustPDB with legacy naming.\n";
 $TESTNUM++;
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("adjusted.pdb");
 `$exec -legacy -i $origpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("adjusted.pdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if ($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"adjusted.pdb") or internal_error("Unable to open adjusted.pdb");
  while(<NEW>){
   $LINESnew++;
  }
  if($LINESnew==$LINESorig+$gitshift){
   # +2 because a comment is added at the top
   # but, we are also removing 2 lines, since they are consecutive TER lines
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i adjusted.pdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }

 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts","smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 2 
 print "\tChecking smog_adjustPDB with user-specified file name (legacy).\n";
 $TESTNUM++;
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -legacy -i $origpdb -o $newpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open adjusted.pdb");
  while(<NEW>){
   $LINESnew++;
  }
  if($LINESnew==$LINESorig+$gitshift){
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 3
 print "\tChecking smog_adjustPDB with default exact matching.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.resnames.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  if($LINESnew==$LINESorig+$gitshift){
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 4
 print "\tChecking smog_adjustPDB with exact matching and alternate names.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.pdb";
 my $mapfile="$mapdir/sbmMapExact.alts";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -map $mapfile -i $origpdb -o $newpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  if($LINESnew==$LINESorig+$gitshift){
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 5
 print "\tChecking smog_adjustPDB with exact matching, alternate names and last tags.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.lf.pdb";
 my $mapfile="$mapdir/sbmMapExact.lf.alts";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -map $mapfile -i $origpdb -o $newpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $LINESorig2=2; # the 2 is because we are adding two lines when running adjust
  open(ORIG2,"$origpdb") or internal_error("Unable to open $origpdb");
  while(<ORIG2>){
   $LINESorig2++;
  }
  if($LINESnew==$LINESorig2+$gitshift){
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 6
 print "\tChecking smog_adjustPDB with exact matching, alternate names and -large format.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.pdb";
 my $mapfile="$mapdir/sbmMapExact.alts";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 removeifexists("$newpdb");
 `$exec -map $mapfile -i $origpdb -o $newpdb -large &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   my $LINE=$_;
   if ($LINE =~ m/^LARGE/){
    $FAIL{'LARGE'}=0;
   }
   $LINESnew++;
  }
  if($LINESnew==$LINESorig+1+$gitshift){
   $FAIL{"FILE LENGTH"}=0;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 7
 print "\tChecking smog_adjustPDB with default exact matching, removewater and official/dirty PDB.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/2ci2.official.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 $FAIL{"FILE LENGTH"}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb -removewater &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 8
 print "\tChecking smog_adjustPDB with default exact matching, removewater, official/dirty PDB and PDBresnum.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/2ci2.official.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 $FAIL{"FILE LENGTH"}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb -removewater -PDBresnum &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 9
 print "\tChecking smog_adjustPDB with default exact matching, removeH, PDBresnum, near-official/dirty PDB.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/1cis.nearofficial.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 $FAIL{"FILE LENGTH"}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb -removeH -PDBresnum &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 9
 print "\tChecking smog_adjustPDB with default exact matching, near-official/dirty PDB with insertions.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/6qnr-partIns.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 $FAIL{"FILE LENGTH"}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb  &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 10
 print "\tChecking smog_adjustPDB with default exact matching, near-official/dirty PDB, insertions and sorting.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/6qnr-partIns.pdb";
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 $FAIL{"FILE LENGTH"}=-1;
 removeifexists("$newpdb");
 `$exec -i $origpdb -o $newpdb -sort &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open $newpdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  print "\n";
  clearfiles(("adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }



 return ($FAILSUM, $printbuffer);

}

sub trueifexists
{
 my ($file)=@_;
 if(-e $file){
  return 0;
 }else{
  return 1;
 }
}

return 1;
