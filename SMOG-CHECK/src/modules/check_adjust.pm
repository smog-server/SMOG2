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

 my $sharerefs="share/adjustrefs";
 my %TESTED;
 my $TESTED=\%TESTED;
 open(ORIG,"$origpdb") or internal_error("Unable to open $origpdb");
 while(<ORIG>){
  $LINESorig++;
 }
 my @FAILLIST = ('NON-ZERO EXIT','OUTPUT NAME','SMOG RUNS','LARGE','IDENTICAL');

 # TEST 1
 print "\tChecking smog_adjustPDB with legacy naming.\n";
 $TESTNUM++;
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("adjusted.pdb");
 `$exec -default -legacy -i $origpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("adjusted.pdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if ($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"adjusted.pdb") or internal_error("Unable to open adjusted.pdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i adjusted.pdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","adjusted.pdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }

 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts","smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 2 
 print "\tChecking smog_adjustPDB with user-specified file name (legacy).\n";
 $TESTNUM++;
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -legacy -i $origpdb -o $newpdb &> output.$tool`;
 $FAIL{"OUTPUT NAME"}=trueifexists("$newpdb");

 $FAIL{"NON-ZERO EXIT"}=$?;
 if($FAIL{"NON-ZERO EXIT"} == 0){
  my $LINESnew=0;
  open(NEW,"$newpdb") or internal_error("Unable to open adjusted.pdb");
  while(<NEW>){
   $LINESnew++;
  }
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 3
 print "\tChecking smog_adjustPDB with default exact matching.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.resnames.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 4
 print "\tChecking smog_adjustPDB with default exact matching and altlocs.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.resnames.altlocs.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb -altloc &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 5
 print "\tChecking smog_adjustPDB with exact matching and alternate names.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.pdb";
 my $mapfile="$mapdir/sbmMapExact.alts";
 &testsperformed($TESTED,\%FAIL);
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
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 6
 print "\tChecking smog_adjustPDB with exact matching, alternate names and last tags.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.lf.pdb";
 my $mapfile="$mapdir/sbmMapExact.lf.alts";
 &testsperformed($TESTED,\%FAIL);
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
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 7
 print "\tChecking smog_adjustPDB with exact matching, alternate names and -large format.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.atomnames.pdb";
 my $mapfile="$mapdir/sbmMapExact.alts";
 &testsperformed($TESTED,\%FAIL);
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
  my $smogout=`$smogexec -AA -i $newpdb -dname adjusted &> smog.output`;
  $FAIL{'SMOG RUNS'}=$?;
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 8
 print "\tChecking smog_adjustPDB with default exact matching, removewater and official/dirty PDB.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/2ci2.official.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb -removewater &> output.$tool`;
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

  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 9 
 print "\tChecking smog_adjustPDB with default exact matching, removewater, official/dirty PDB and PDBresnum.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/2ci2.official.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb -removewater -PDBresnum &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }

 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 10
 print "\tChecking smog_adjustPDB with default exact matching, removeH, PDBresnum, near-official/dirty PDB.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/1cis.nearofficial.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb -removeH -PDBresnum &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }

 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 11
 print "\tChecking smog_adjustPDB with default exact matching, near-official/dirty PDB with insertions.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/6qnr-partIns.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb  &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }

 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 12
 print "\tChecking smog_adjustPDB with default exact matching, near-official/dirty PDB, insertions and sorting.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/6qnr-partIns.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `$exec -default -i $origpdb -o $newpdb -sort &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }

 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 13
 print "\tChecking smog_adjustPDB with default exact matching and interactive selection.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.resnames.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `echo D | $exec -i $origpdb -o $newpdb &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }

 # TEST 14
 print "\tChecking smog_adjustPDB with default exact matching and Interactive selection.\n";
 $TESTNUM++;
 my $origpdb="$pdbdir/mangled.resnames.pdb";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 $FAIL{'LARGE'}=-1;
 removeifexists("$newpdb");
 `echo "I\n0\n" | $exec -i $origpdb -o $newpdb &> output.$tool`;
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
  my $compsummary;
  ($FAIL{'IDENTICAL'},$compsummary)=compareFiles("$sharerefs/adjusted.ref.$TESTNUM.pdb","$newpdb");
  if($FAIL{'IDENTICAL'} != 0){
   open(TMP,">output.check.$tool") or internal_error("Unable to open output file output.$tool");
   print TMP "$compsummary\n";
   close(TMP);
  }
 }
 my ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed($TESTNUM,("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.check.$tool","adjusted.pdb","$newpdb","output.$tool","adjusted.gro","adjusted.top","adjusted.ndx","adjusted.contacts" ,"smog.output"));
 }



 $FAILSUM+=checkalltested(\@FAILLIST,\%FAIL);

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


sub compareFiles
{
 my ($reffile,$newfile)=@_;
 open(REF,"$reffile") or internal_error("Unable to open reference file $reffile.");
 open(NEW,"$newfile") or return (1,"Failed to generate the output file: $newfile");

 my @REFDATA;
 while(<REF>){
  my $LINE=$_;
  chomp($LINE);
  if($LINE !~ m/^REMARK/){
   push(@REFDATA,$LINE);
  }
 }

 my @NEWDATA;
 while(<NEW>){
  my $LINE=$_;
  chomp($LINE);
  if($LINE !~ m/^REMARK/){
   push(@NEWDATA,$LINE);
  }
 }

 my $length=$#REFDATA;
 if($#NEWDATA == 0 or $length == 0){
  return (1,"Generated file has no content: $newfile.");
 }

 if($length != $#NEWDATA){
  return (1,"Wrong number of lines in $newfile.");
 }

 my $match=-1;
 my $mismatchdata="";
 for(my $I=0;$I<=$length;$I++){
  if($REFDATA[$I] eq $NEWDATA[$I] and $REFDATA[$I] ne ""){
   $match++;
  }else{
   $mismatchdata .= "\nInconsistent data in ref and generated file:\nRef:$REFDATA[$I]\nGen:$NEWDATA[$I]\n";
  }
 }
 if($length == $match){
  return(0,"");
 }else{
  return(1,"$mismatchdata");
 }
}

return 1;
