package check_extract;
use strict;
use Exporter;
use smog_common;
use check_common;
our @ISA = 'Exporter';
our @EXPORT = qw(check_extract);

# in supported_directives:0 = not supported. 1 = required once. 2=required (may appear more than once). 3=optional once.
our %supported_directives = ( 'defaults' => '1',
        'atomtypes' => '1',
        'nonbond_params' => '0',
        'moleculetype' => '1',
        'atoms' => '1',
        'bonds' => '1',
        'angles' => '1',
        'dihedrals' => '1',
        'pairs' => '1',
        'exclusions' => '1',
        'system' => '1',
        'molecules' => '1',
        'position_restraints' => '3'
        );



sub check_extract
{
 my ($exec,$pdbdir,$CHECKGMX,$GMXVER,$GMXPATH,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA)=@_;
 my $NFAIL=0;
 my $MESSAGE="";
 my %FAIL;
 my $FAILED;
 my $FAILSUM=0;
 my $FATAL;
 my $UNINIT;
 my $printbuffer;
 my $tool="extract";
 my @FAILLIST = ('NON-ZERO EXIT','EXTRA MAP FILE GENERATED','GMX COMPATIBLE');

 %FAIL=resettests(\%FAIL,\@FAILLIST);

# generate an AA model RNA 
 `smog2 -i $pdbdir/tRNA.pdb -AA -dname AA.tmp > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_extract.");
 }else{
  clearfiles("output.smog");
 }
  print "\tChecking smog_extract with all-atom model: no restraints\n";
  for(my $group=0;$group<3;$group++){

   %FAIL=resettests(\%FAIL,\@FAILLIST);
   print "\tChecking with index group $group\n";
   `echo $group | $exec -f AA.tmp.top -g AA.tmp.gro -n $pdbdir/sample.AA.ndx  &> output.$tool`;

   $FAIL{"NON-ZERO EXIT"}=$?;
   $FAIL{"EXTRA MAP FILE GENERATED"} = checkrestraintfile(1,"restrained.map");
   $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","extracted","no","extracted","noG96","no");
   my $string=loadfile("extracted.top");
   my ($DATA,$DIRLIST)=checkdirectives($string);

   ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
   $FAILSUM += $FAILED;
   if($FAILED !=0){
    savefailed("AA.nores.$group",("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
    print "$printbuffer\n";
   }else{
    clearfiles(("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
   }
  } 
  clearfiles(("AA.tmp.top","AA.tmp.gro","AA.tmp.ndx","AA.tmp.contacts"));

  print "\tChecking smog_extract with all-atom model: no restraints: non-standard fields\n";
  for(my $group=0;$group<2;$group++){

   %FAIL=resettests(\%FAIL,\@FAILLIST);
   print "\tChecking with index group $group\n";
   `echo $group | $exec -f $pdbdir/large.top -g $pdbdir/large.gro -n $pdbdir/large.ndx  &> output.$tool`;
   $FAIL{"NON-ZERO EXIT"}=$?; 
   $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","extracted","no","extracted","noG96","no");
   $FAIL{"EXTRA MAP FILE GENERATED"} = checkrestraintfile(1,"restrained.map");
   my $string=loadfile("extracted.top");
   my ($DATA,$DIRLIST)=checkdirectives($string);

   ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
   $FAILSUM += $FAILED;
   if($FAILED !=0){
    savefailed("AA.nores.nonstandard.$group",("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
    print "$printbuffer\n";
   }else{
    clearfiles(("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
   }
  } 

  print "\tChecking smog_extract with all-atom model: no restraints: non-standard fields: ndxorder on\n";
  for(my $group=0;$group<2;$group++){

   %FAIL=resettests(\%FAIL,\@FAILLIST);
   print "\tChecking with index group $group\n";
   `echo $group | $exec -f $pdbdir/large.top -g $pdbdir/large.gro -n $pdbdir/large.ndx -ndxorder &> output.$tool`;

   $FAIL{"NON-ZERO EXIT"}=$?;
   $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","extracted","no","extracted","noG96","yes");
   $FAIL{"EXTRA MAP FILE GENERATED"} = checkrestraintfile(1,"restrained.map");
   my $string=loadfile("extracted.top");
   my ($DATA,$DIRLIST)=checkdirectives($string);

   ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
   $FAILSUM += $FAILED;
   if($FAILED !=0){
    savefailed("AA.nores.nonstandard.$group",("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
    print "$printbuffer\n";
   }else{
    clearfiles(("output.$tool","extracted.top","extracted.gro","atomindex.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
   }
  } 

  print "\tChecking smog_extract with all-atom model: restraints: non-standard fields\n";
  for(my $group=0;$group<2;$group++){

   %FAIL=resettests(\%FAIL,\@FAILLIST);
   print "\tChecking with index group $group\n";
   `echo $group | $exec -f $pdbdir/large.top -g $pdbdir/large.gro -n $pdbdir/large.ndx -restraints 100 &> output.$tool`;
   $FAIL{"NON-ZERO EXIT"}=$?;
   $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","extracted","no","extracted","noG96","yes");
   $FAIL{"EXTRA MAP FILE GENERATED"} = checkrestraintfile(0,"restrained.map");
   my $string=loadfile("extracted.top");
   my ($DATA,$DIRLIST)=checkdirectives($string);

   ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
   $FAILSUM += $FAILED;
   if($FAILED !=0){
    savefailed("AA.res.nonstandard.$group",("output.$tool","extracted.top","extracted.gro","atomindex.map","restrained.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
    print "$printbuffer\n";
   }else{
    clearfiles(("output.$tool","extracted.top","extracted.gro","atomindex.map","restrained.map","topol.tpr","extracted.box.gro","extracted.editconf","extracted.grompp","extracted.out.mdp"));
   }
  } 

 return ($FAILSUM, $printbuffer);

}

return 1;
