package check_scale;
use strict;
use Exporter;
use smog_common;
use check_common;
our @ISA = 'Exporter';
our @EXPORT = qw(check_scale);

sub check_scale
{
 my ($exec,$pdbdir,$CHECKGMX,$GMXVER,$GMXPATH,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA)=@_;
 my $NFAIL=0;
 my $MESSAGE="";
 my %FAIL;
 my $FAILED;
 my $FAILSUM=0;
 my $tool="scale";
 my $printbuffer="";
 my @FAILLIST = ('NON-ZERO EXIT','UNCHANGED DIRECTIVES','N DIHEDRALS','N SCALED DIHEDRALS','N CONTACTS','N SCALED CONTACTS','GMX COMPATIBLE');
 %FAIL=resettests(\%FAIL,\@FAILLIST);

# generate an AA model RNA 
 `smog2 -i $pdbdir/tRNA.pdb -AA -dname AA.tmp > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_ions.");
 }else{
  clearfiles("output.smog");
 }

 print "\tChecking smog_scale-energies with all-atom model: rescaling terms\n";

 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $indexfile="share/PDB.files/sample.AA.ndx";
 my $grpsel="$pdbdir/in.groups";
 my $RC=1.5;
 my $RD=1.2;
 `$exec -f AA.tmp.top -n $indexfile -rc $RC -rd $RD < $grpsel &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","smog.rescaled","no","AA.tmp","noG96","no");
 my ($samedirs,$dihlength,$dihmatch,$conlength,$conmatch)=comparetopsrescale("AA.tmp.top","smog.rescaled.top",$indexfile,$grpsel,$RC,$RD);
 $FAIL{"UNCHANGED DIRECTIVES"}=$samedirs;
 $FAIL{"N DIHEDRALS"}=$dihlength;
 $FAIL{"N SCALED DIHEDRALS"}=$dihmatch;
 $FAIL{"N CONTACTS"}=$conlength;
 $FAIL{"N SCALED CONTACTS"}=$conmatch;

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top"){
   `cp $file tmp`;
  }
  savefailed(1,("output.$tool","smog.rescaled.top","smog.rescaled.top","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","topol.tpr","smog.rescaled.box.gro","smog.rescaled.editconf","smog.rescaled.grompp","smog.rescaled.out.mdp"));
  print "$printbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","smog.rescaled.top","topol.tpr","smog.rescaled.box.gro","smog.rescaled.editconf","smog.rescaled.grompp","smog.rescaled.out.mdp"));
 }

 print "\tChecking smog_scale-energies with all-atom model: removing terms\n";

 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $indexfile="share/PDB.files/sample.AA.ndx";
 my $grpsel="$pdbdir/in.groups";
 my $outfile="test";
 my $RC=0;
 my $RD=0;
 `$exec -f AA.tmp.top -of "$outfile.top" -n $indexfile -rc $RC -rd $RD < $grpsel &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $FAIL{"GMX COMPATIBLE"}=runGMX("AA",$CHECKGMX,"no",$GMXEDITCONF,$GMXPATH,"",$GMXEXEC,$GMXMDP,$GMXMDPCA,"no","$outfile","no","AA.tmp","noG96","no");
 my ($samedirs,$dihlength,$dihmatch,$conlength,$conmatch)=comparetopsrescale("AA.tmp.top","$outfile.top",$indexfile,$grpsel,$RC,$RD);
 $FAIL{"UNCHANGED DIRECTIVES"}=$samedirs;
 $FAIL{"N DIHEDRALS"}=$dihlength;
 $FAIL{"N SCALED DIHEDRALS"}=$dihmatch;
 $FAIL{"N CONTACTS"}=$conlength;
 $FAIL{"N SCALED CONTACTS"}=$conmatch;
 
 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  savefailed(2,("output.$tool","smog.rescaled.top","smog.rescaled.top","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","topol.tpr","$outfile.top","$outfile.box.gro","$outfile.editconf","$outfile.grompp","$outfile.out.mdp"));
  print "$printbuffer\n";
 }else{
  clearfiles(("output.$tool","$outfile.top","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","topol.tpr","$outfile.box.gro","$outfile.editconf","$outfile.grompp","$outfile.out.mdp"));
 }
 return ($FAILSUM, $printbuffer);

}

sub comparetopsrescale
{
 my ($old,$new,$indexFile,$grpsel,$RC,$RD)=@_;
 # read in original and new top files
 my $string=loadfile("$old");
 my ($DATA,$DIRLIST)=checkdirectives($string);
 my %DIRLIST=%{$DIRLIST};
 my @DATA=@{$DATA};
 $string=loadfile("$new");
 my ($DATA2,$DIRLIST2)=checkdirectives($string);
 my %DIRLIST2=%{$DIRLIST2};
 my @DATA2=@{$DATA2};
 my $samedirs=0; 
 my $dihlength=1;
 my $dihnum=1;
 my $dscaled=0; 
 my $sameatoms=0; 
 my $conlength=1;
 my $connum=1;
 my $cscaled=0; 
 my $consameatoms=0; 

 # check that unchanged directives remain unchanged
 foreach my $DIR("defaults","atomtypes","moleculetype","atoms","bonds","angles","molecules")
 {
  $samedirs++;
  if($DATA[$DIRLIST{"$DIR"}] eq $DATA2[$DIRLIST2{"$DIR"}]){
   $samedirs--;
  }else{
   print "issue: directive $DIR changed, but it shouldn\'t\n";
  }
 }
 if($RC!=0){
  $samedirs++;
  if($DATA[$DIRLIST{"exclusions"}] eq $DATA2[$DIRLIST2{"exclusions"}]){
   $samedirs--;
  }else{
   print "issue: directive exclusions changed, but it shouldn\'t\n";
  }
 }

 # read in the ndx file
 print "\t";
 my ($Ngrps,$grpnms,$groupnames,$atomgroup) = readindexfile($indexFile);
 
 my @grpnms=@{$grpnms};
 my %groupnames=%{$groupnames};
 my %atomgroup=%{$atomgroup};

 # read in grp selection list
 open(GRPSEL,"$grpsel") or smog_quit("unable to open $grpsel");
 my $DGROUP=<GRPSEL>;
 chomp($DGROUP);
 $DGROUP=$grpnms[$DGROUP];
 my $CGROUP1=<GRPSEL>;
 chomp($CGROUP1);
 $CGROUP1=$grpnms[$CGROUP1];
 my $CGROUP2=<GRPSEL>;
 chomp($CGROUP2);
 $CGROUP2=$grpnms[$CGROUP2];

 # check dihedrals
 my @D1 = split(/\n/,$DATA[$DIRLIST{"dihedrals"}]);
 my @D2 = split(/\n/,$DATA2[$DIRLIST2{"dihedrals"}]);
 if($RD != 0){
  if($#D1 == $#D2){
   $dihlength=0;
  }
  $dihnum=$#D1+1;
  # rescaling dihedrals
  my $commentshift1=0;
  my $commentshift2=0;
  for(my $I=1;$I<=$#D1;$I++){
   until(! exists $D1[$I+$commentshift1] or  hascontent($D1[$I+$commentshift1])){
    $commentshift1++;
    $dihnum--;
   }
   until(! exists $D2[$I+$commentshift2] or  hascontent($D2[$I+$commentshift2])){
    $commentshift2++;
   }
   my @A1=split(/\s+/,$D1[$I+$commentshift1]); 
   my @A2=split(/\s+/,$D2[$I+$commentshift2]);
    if($A1[0]==$A2[0] && $A1[1]==$A2[1] && $A1[2]==$A2[2] && $A1[3]==$A2[3]){
     $sameatoms++;
     my $rescale=0;
     for(my $J=0;$J<4;$J++){
      if(exists $atomgroup{$DGROUP}{$A1[$J]}){
       $rescale++; 
      }     
     }
     my $resc;
     if($rescale==4 && $A1[4] != 2){
      $resc=$RD;
     }else{
      $resc=1.0
     }
     if(abs($A1[6]*$resc-$A2[6])<0.001){
      $dscaled++;
     }else{
      print "issue: dihedral not scaled properly. Should be scaled by $resc.\nold $D1[$I+$commentshift1]\nnew $D2[$I+$commentshift2]\n";
     } 
    }else{
     print "issue: atom numbers don\'t match:\nold $D1[$I+$commentshift1]\nnew $D2[$I+$commentshift2]\n";
    } 
  }
 }else{
  # this means we should check if the dihedrals are removed.
  $dihlength=-1;
  $dihnum=$#D1+1;
  # rescaling dihedrals
  my $I2=1;
  my $commentshift1=0;
  my $commentshift2=0;
  for(my $I=1;$I<=$#D1;$I++){
   until(! exists $D1[$I+$commentshift1] or  hascontent($D1[$I+$commentshift1])){
    $commentshift1++;
    $dihnum--;
   }
   my @A1=split(/\s+/,$D1[$I+$commentshift1]);

   my $remove=0;
   for(my $J=0;$J<4;$J++){
    if(exists $atomgroup{$DGROUP}{$A1[$J]}){
     $remove++; 
    }     
   }
   if($remove==4 && $A1[4] == 1){
    # this one should not be in the new file
    $dihnum--;
    next;
   }
   until(! exists $D2[$I2+$commentshift2] or  hascontent($D2[$I2+$commentshift2])){
    $commentshift2++;
   }
   my @A2=split(/\s+/,$D2[$I2+$commentshift2]);

   if($A1[0]==$A2[0] && $A1[1]==$A2[1] && $A1[2]==$A2[2] && $A1[3]==$A2[3]){
    $sameatoms++;
   }else{
    print "issue: atom numbers don\'t match:\nold $D1[$I+$commentshift1]\nnew $D2[$I2+$commentshift2]\n";
   } 
   if(abs($A1[6]-$A2[6])<0.001){
    $dscaled++;
   }else{
    print "issue: dihedral not set properly.\n$I, $commentshift1, $I2, $commentshift2\nold $D1[$I+$commentshift1]\nnew $D2[$I2+$commentshift2]\n";
   } 
   $I2++;
  }
 }

 # check contacts
 my @C1 = split(/\n/,$DATA[$DIRLIST{"pairs"}]);
 my @C2 = split(/\n/,$DATA2[$DIRLIST2{"pairs"}]);
 if($RC != 0){
  if($#C1 == $#C2){
   $conlength=0;
  }
  $connum=$#C1+1;
  # rescaling contacts
  my $commentshift1=0;
  my $commentshift2=0;
  for(my $I=1;$I<=$#C1;$I++){
   until(! exists $C1[$I+$commentshift1] or  hascontent($C1[$I+$commentshift1])){
    $commentshift1++;
    $connum--;
   }
   until(! exists $C2[$I+$commentshift2] or  hascontent($C2[$I+$commentshift2])){
    $commentshift2++;
   }
   my @A1=split(/\s+/,$C1[$I+$commentshift1]); 
   my @A2=split(/\s+/,$C2[$I+$commentshift2]);

    if($A1[0]==$A2[0] && $A1[1]==$A2[1]){
     $consameatoms++;
     my $rescale=0;
      if(exists $atomgroup{$CGROUP1}{$A1[0]} && exists $atomgroup{$CGROUP2}{$A1[1]}){
       $rescale=1; 
      }elsif(exists $atomgroup{$CGROUP2}{$A1[0]} && exists $atomgroup{$CGROUP1}{$A1[1]}){
       $rescale=1; 
      }
     my $resc;
     if($rescale==1 && $A1[2] == 1){
      $resc=$RC;
     }else{
      $resc=1.0
     }
     if(abs($A1[3]*$resc-$A2[3])<0.001 && abs($A1[4]*$resc-$A2[4])<0.001){
      $cscaled++;
     }else{
      print "issue: contact not scaled properly. Should be scaled by $resc.\nold $C1[$I]\nnew $C2[$I]\n";
     } 
    }else{
     print "issue: atom numbers don\'t match:\nold $C1[$I]\nnew $C2[$I]\n";
    } 
  }
 }else{
  # this means we should check if the contacts are removed.
  $conlength=-1;
  $connum=$#C1+1;
  # rescaling contacts
  my $I2=1;
#  for(my $I=0;$I<=$#C1;$I++){
#   my @A1=split(/\s+/,$C1[$I]); 
#   my @A2=split(/\s+/,$C2[$I2]);

  my $commentshift1=0;
  my $commentshift2=0;
  for(my $I=1;$I<=$#C1;$I++){
   until(! exists $C1[$I+$commentshift1] or  hascontent($C1[$I+$commentshift1])){
    $commentshift1++;
    $connum--;
   }
   my @A1=split(/\s+/,$C1[$I+$commentshift1]); 

   my $remove=0;
   if(exists $atomgroup{$CGROUP1}{$A1[0]} && exists $atomgroup{$CGROUP2}{$A1[1]}){
    $remove=1; 
   }elsif(exists $atomgroup{$CGROUP2}{$A1[0]} && exists $atomgroup{$CGROUP1}{$A1[1]}){
    $remove=1; 
   }
   if($remove==1){
    $connum--;
    next;
   }

   until(! exists $C2[$I2+$commentshift2] or  hascontent($C2[$I2+$commentshift2])){
    $commentshift2++;
   }
   my @A2=split(/\s+/,$C2[$I2+$commentshift2]);

   if($A1[0]==$A2[0] && $A1[1]==$A2[1]){
    $consameatoms++;
    if(abs($A1[3]-$A2[3])<0.001 && abs($A1[4]-$A2[4])<0.001){
     $cscaled++;
    }else{
     print "issue: contact not set properly.\nold $C1[$I+$commentshift1]\nnew $C2[$I2+$commentshift2]\n";
    } 
   }else{
    print "issue: atom numbers don\'t match:\nold $C1[$I+$commentshift1]\nnew $C2[$I2+$commentshift2]\n";
   } 
  $I2++;
  }
 }
 return ($samedirs,$dihlength,abs($dihnum-$dscaled),$conlength,abs($connum-$cscaled));
}

return 1;
