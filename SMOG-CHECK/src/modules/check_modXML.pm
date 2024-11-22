package check_modXML;
use strict;
use Exporter;
use smog_common;
use check_common;
use OpenSMOG;
our @ISA = 'Exporter';
our @EXPORT = qw(check_modXML);

sub check_modXML
{
 my ($exec,$pdbdir)=@_;
 my $NFAIL=0;
 my $MESSAGE="";
 my %FAIL;
 my $FAILED;
 my $FAILSUM=0;
 my $tool="modifyXML";
 my $printbuffer="";
 my $testnum=0;
 my @FAILLIST = ('NON-ZERO EXIT','XML TREE','CONSTANTS','CONSTANTS EXIST','CONTACTS EXIST','DIHEDRALS EXIST','PARAM LISTS','EXPRESSION','INTERACTION COUNT: CONTACTS','INTERACTION VALUES: CONTACTS','INTERACTION COUNT: DIHEDRALS','INTERACTION VALUES: DIHEDRALS','NONBOND EXIST','EXPRESSION: NONBOND','NONBOND PARAM VALUES','XML OpenSMOGversion');
 my %TESTED;
 my $TESTED=\%TESTED;
# generate an AA model RNA 
 `smog2 -i $pdbdir/tRNA.pdb -AA -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }
 $testnum++;
 print "\tChecking interactive call: test $testnum\n";
 my $tmpbuffer="";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $indexfile="share/PDB.files/xml.test.ndx";
 my $settings="share/PDB.files/xmlsettings.1.in";
 my ($settings,$conhash,$dihhash)=processsettings($settings);
 `echo "$settings" | $exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml  &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml"));
 }
 
 $testnum++;
 print "\tChecking interactive call: test $testnum\n";
 my $tmpbuffer="";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $settings="share/PDB.files/xmlsettings.2.in";
 my ($settings,$conhash,$dihhash)=processsettings($settings);
 `echo "$settings" | $exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml  &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml","AA.tmp.contacts","AA.tmp.gro","AA.tmp.ndx","AA.tmp.top","AA.tmp.xml"));
 }

 $testnum++;
 print "\tChecking interactive call: test $testnum\n";
# generate an AA model protein 
 `smog2 -i $pdbdir/1AKEapo_v2.ion.pdb -t share/templates/Ion-test -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }

 my $tmpbuffer="";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $settings="share/PDB.files/xmlsettings.3.in";
 my ($settings,$conhash,$dihhash)=processsettings($settings);
 `echo "$settings" | $exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml  &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml","AA.tmp.contacts","AA.tmp.gro","AA.tmp.ndx","AA.tmp.top","AA.tmp.xml"));
 }

 $testnum++;
 print "\tChecking interactive call: test $testnum\n";
# generate an AA model protein 
 `smog2 -i $pdbdir/1AKEapo_v2.ion.pdb -t share/templates/Ion-test -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }

 my $tmpbuffer="";
 my $indexfile="share/PDB.files/xml.test.small.ndx";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $settings="share/PDB.files/xmlsettings.4.in";
 my ($settings,$conhash,$dihhash)=processsettings($settings);
 `echo "$settings" | $exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml  &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml","AA.tmp.contacts","AA.tmp.gro","AA.tmp.ndx","AA.tmp.top","AA.tmp.xml"));
 }

 $testnum++;
 print "\tChecking command-line call - contacts: test $testnum\n";
# generate an AA model protein 
 `smog2 -i $pdbdir/1AKEapo_v2.ion.pdb -t share/templates/Ion-test -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }

 my $tmpbuffer="";
 my $indexfile="share/PDB.files/xml.test.ndx";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $settings="share/PDB.files/xmlsettings.5.in";
 my ($settings,$conhash,$dihhash)=processsettingscl($settings);
 `$exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml $settings &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml","AA.tmp.contacts","AA.tmp.gro","AA.tmp.ndx","AA.tmp.top","AA.tmp.xml"));
 }

 $testnum++;
 print "\tChecking command-line call - dihedrals: test $testnum\n";
# generate an AA model protein 
 `smog2 -i $pdbdir/1AKEapo_v2.ion.pdb -t share/templates/Ion-test -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }

 my $tmpbuffer="";
 my $indexfile="share/PDB.files/xml.test.ndx";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $settings="share/PDB.files/xmlsettings.6.in";
 my ($settings,$conhash,$dihhash)=processsettingscl($settings);
 `$exec -OpenSMOG AA.tmp.xml -n $indexfile -OpenSMOGout AA.tmp.out.xml $settings &> output.$tool`;
 $FAIL{"NON-ZERO EXIT"}=$?;
 $tmpbuffer .= compareXMLsmodify(\%FAIL,"AA.tmp.xml","AA.tmp.out.xml",$indexfile,$conhash,$dihhash);

 &testsperformed($TESTED,\%FAIL);

 ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
 $FAILSUM += $FAILED;
 if($FAILED !=0){
  `mkdir tmp`;
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `cp $file tmp`;
  }
  savefailed($testnum,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml","AA.tmp.contacts","AA.tmp.gro","AA.tmp.ndx","AA.tmp.top","AA.tmp.xml"));
 }

 $FAILSUM+=checkalltested(\@FAILLIST,$TESTED);

 return ($FAILSUM, $printbuffer);

}

sub processsettings{
 my ($settings)=@_;
 # define some changes we want to introduce. This will be moved to a file for reading, later
 my %conhash;
 my %dihhash;
 open(INSET,"$settings") or smog_quit("can not open $settings");
 $settings ="";
 while(<INSET>){
  my %parhash;
  my @grparr;
  my $LINE=$_;
  chomp($LINE);
  $LINE =~ s/^\s+//g;
  $LINE =~ s/\s+$//g;
  my @A=split(/\s+/,$LINE);
  my $last="";
  my $name="";
  my $type="";
  for(my $I=0;$I<$#A;$I++){
   if($A[$I] eq "contacts" || $A[$I] eq "dihedrals"){
    $last=$A[$I];
    $type=$A[$I];
    $settings .= "Y\n";
    if($A[$I] eq "contacts"){
     $I++;
    }else{
     $I++;
     $settings .= "  $A[$I]  \n";
    }
    $name=$A[$I];
   }elsif($A[$I] eq "grp"){
    $last=$A[$I];
    $I++;
    $settings .= "  $A[$I]  \n";
    push(@grparr,$A[$I]);
   }elsif($A[$I] eq "modremove"){
    $last=$A[$I];
    $I++;
    $settings .= "  $A[$I]  \n";
    if($A[$I] =~ m/R/){
     $conhash{$name}->{"modremove"}=0;
    }
   }elsif($A[$I] eq "param"){
    if($last eq "param"){
     $settings .= "Y\n";
    }
    $last=$A[$I];
    $I++;
    $settings .= "  $A[$I]  \n";
    $I++;
    $settings .= "  $A[$I]  \n";
    $parhash{$A[$I-1]}=$A[$I];
   }
  }
  $settings .= "N\n";
  if($type eq "contacts"){
   $conhash{$name}->{"parms"}=\%parhash;
   $conhash{$name}->{"groups"}=\@grparr;
  }elsif($type eq "dihedrals"){ 
   $dihhash{$name}->{"parms"}=\%parhash;
   $dihhash{$name}->{"groups"}=\@grparr;
  }
 }
 $settings .= "N\n";
 return ($settings,\%conhash,\%dihhash);
}

sub processsettingscl{
 my ($settings)=@_;
 # define some changes we want to introduce. This will be moved to a file for reading, later
 my %conhash;
 my %dihhash;
 open(INSET,"$settings") or smog_quit("can not open $settings");
 $settings ="";
 while(<INSET>){
  my %parhash;
  my @grparr;
  my $LINE=$_;
  chomp($LINE);
  $LINE =~ s/^\s+//g;
  $LINE =~ s/\s+$//g;
  my @A=split(/\s+/,$LINE);
  my $last="";
  my $name="";
  my $type="";
  my $grpnum=1;
  for(my $I=0;$I<$#A;$I++){
   if($A[$I] eq "contacts" || $A[$I] eq "dihedrals"){
    $last=$A[$I];
    $type=$A[$I];
    $settings .= "-inter $A[$I] ";
    $I++;
    $settings .= "-type $A[$I] ";
    $name=$A[$I];
   }elsif($A[$I] eq "grp"){
    $last=$A[$I];
    $I++;
    $settings .= "-grp$grpnum $A[$I] ";
    $grpnum++;
    $I++;
    push(@grparr,$A[$I]);
   }elsif($A[$I] eq "param"){
    $last=$A[$I];
    $I++;
    $settings .= "-param $A[$I] ";
    $I++;
    $settings .= "-modby \"$A[$I]\" ";
    $parhash{$A[$I-1]}=$A[$I];
   }
  }
  $settings .= "\n";
  if($type eq "contacts"){
   $conhash{$name}->{"parms"}=\%parhash;
   $conhash{$name}->{"groups"}=\@grparr;
  }elsif($type eq "dihedrals"){ 
   $dihhash{$name}->{"parms"}=\%parhash;
   $dihhash{$name}->{"groups"}=\@grparr;
  }
 }
 return ($settings,\%conhash,\%dihhash);
}



sub compareXMLsmodify
{
 my ($fail,$old,$new,$indexFile,$conhash,$dihhash)=@_;
 # read in original and new top files
 my $xmlold=readOpenSMOGxml($old);
 my $xmlnew=readOpenSMOGxml($new);
 my $pbuffer="";
 $pbuffer .= checkheadparams($fail,$xmlold,$xmlnew);
 # read in the ndx file
 print "\t";
 my ($Ngrps,$grpnms,$groupnames,$atomgroup) = readindexfile($indexFile);
 
 my @grpnms=@{$grpnms};
 my %atomgroup=%{$atomgroup};
 $pbuffer .=checkconstants($fail,$xmlold,$xmlnew);
 $pbuffer .=checknonbonded($fail,$xmlold,$xmlnew);

 $pbuffer.=checkcontacts($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$conhash);
 $pbuffer.=checkdihedrals($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$dihhash);

 return $pbuffer;
}

sub checkconstants{ 
 # checks that the constants are unchanged
 my ($fail,$xmlold,$xmlnew)=@_;
 my $printbuffer="";
 if((defined $xmlold->{'constants'} &&  defined $xmlnew->{'constants'}) || (! defined $xmlold->{'constants'} && ! defined $xmlnew->{'constants'})){
  ${$fail}{'CONSTANTS EXIST'}=0;
 }else{
  $printbuffer .= "\"constants\" only found in one of the xml files\n";
 }
 if(defined $xmlold->{'constants'} &&  defined $xmlnew->{'constants'}){
  my $c=0;
  my $m=0;
  my %hold=%{$xmlold->{'constants'}};
  my %hnew=%{$xmlnew->{'constants'}};
  if(scalar keys %hnew == scalar keys %hold){
   foreach my $I(keys %hnew){
    $c++;
    if(defined $hold{$I}){
     if($hold{$I}->{"value"} eq $hnew{$I}->{"value"}){
      $m++;
     }else{
      $printbuffer .= "Parameter $I has different values in the xml files. $hold{$I} $hnew{$I}\n";
     }
    }else{
     $printbuffer .= "Parameter $I not found in new XML file.\n";
    }
   }
  }
  if($c==$m && $c != 0){
   ${$fail}{'CONSTANTS'}=0;
  }else{
   $printbuffer .= "Not all \"constants\" had same values\n";
  }
 }else{
  ${$fail}{'CONSTANTS'}=-1;
 }
 return $printbuffer;
} 

sub checkdihedrals{ 
 # checks that the dihedrals were updated properly
 my ($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$dihhash)=@_;
 my $printbuffer="";
 if((defined $xmlold->{'dihedrals'} &&  defined $xmlnew->{'dihedrals'}) || (! defined $xmlold->{'dihedrals'} && ! defined $xmlnew->{'dihedrals'})){
  ${$fail}{'DIHEDRALS EXIST'}=0;
 }else{
  $printbuffer .= "\"dihedrals\" present in only one xml.\n";
 }

 my $intc=0;
 my $intm=0;
 my $inttc=0;
 my $inttm=0;
 if(defined $xmlold->{'dihedrals'} &&  defined $xmlnew->{'dihedrals'}){
  # set up the index groups and parameter changes
  #...
  my $c=0;
  my $m=0;
  my %hold=%{$xmlold->{'dihedrals'}};
  my %hnew=%{$xmlnew->{'dihedrals'}};
  foreach my $type(keys %{$hold{'dihedrals_type'}}){
   # temp hashes that will store the selected atom groups.
   my %cg0;
   my %cg1;
   my %params;
   if(defined $dihhash->{$type}){
    # found a dihedral type that was supposed to be changed.
    # get the parameters;
    %params=%{$dihhash->{$type}->{parms}};
    my $g0=${$dihhash->{$type}->{groups}}[0];
    %cg0=%{$atomgroup->{${$grpnms}[$g0]}};
   }
   my $interactionhashold=$hold{'dihedrals_type'}->{$type};
   my $interactionhashnew=$hnew{'dihedrals_type'}->{$type};
   my %compnew;
   my %compold;
   # save new xml data
   foreach my $K(@{$interactionhashnew->{'interaction'}}){
    # each entry is a hash
    my %hash=%{$K};
    my $string="";
    foreach my $key(sort keys %hash){
     my $value;
     if($key !~ m/^[ijkl]$/){
      #$value=sprintf("%.2e", $hash{$key});
      $value=truncaten($hash{$key});
     }else{
      $value=$hash{$key};
     }
     $string .= "$key $value ";
    }
    $compnew{$string}=0;
   }

   # deal with original xml and see if it needs to be modified
   if(defined $dihhash->{$type}){
    # update values to the expected numbers
    foreach my $K(@{$interactionhashold->{'interaction'}}){
     # each entry is a hash
     my %hash=%{$K};
     my $string="";

     my $i=$hash{"i"};
     my $j=$hash{"j"};
     my $k=$hash{"k"};
     my $l=$hash{"l"};
     my $mod=1;
     if(defined $cg0{$i} && defined $cg0{$j} && defined $cg0{$k} && defined $cg0{$l}){
      # change something about this dihedral 
      foreach my $key(sort keys %hash){
       my $value;
       if (defined $params{$key}){
        # this is parameter to update
        $value=eval("$hash{$key}$params{$key}");
        #$value=sprintf("%.2e", $value);
        $value=truncaten($value);
       }else{
        $value=$hash{$key};
        if($key !~ m/^[ijkl]$/){
         #$value=sprintf("%.2e", $value);
         $value=truncaten($value);
        }
       }
       $string .= "$key $value ";
      }
     }else{
      # atoms not in groups, just save
      foreach my $key(sort keys %hash){
       my $value;
       if($key !~ m/^[ijkl]$/){
        #$value=sprintf("%.2e", $hash{$key});
        $value=truncaten($hash{$key});
       }else{
        $value=$hash{$key};
       }
       $string .= "$key $value ";
      }
     }
     $compold{$string}=0;
    } 
   }else{
    # not changing this set, so just copy

    foreach my $K(@{$interactionhashold->{'interaction'}}){
     # each entry is a hash
     my %hash=%{$K};
     my $string="";
     foreach my $key(sort keys %hash){
      my $value;
      if($key !~ m/^[ijkl]$/){
       #$value=sprintf("%.2e", $hash{$key});
       $value=truncaten($hash{$key});
      }else{
       $value=$hash{$key};
      }
      $string .= "$key $value ";
     }
     $compold{$string}=0;
    }
   } 
   $intc++;
   if(scalar keys %compold == scalar keys %compnew){
    $intm++;
    foreach my $I(keys %compold){
     $inttc++;
     if(defined $compnew{$I}){
      $inttm++;
     }else{
      $printbuffer .= "Interaction key \"$I\" expected, but not found in new xml\n";
     }
    }
    foreach my $I(keys %compnew){
     $inttc++;
     if(defined $compold{$I}){
      $inttm++;
     }else{
      $printbuffer .= "Interaction key \"$I\" found in new xml, but not expected\n";
     }
    }
 
   }else{
    my $nold=scalar keys %compold;
    my $nnew=scalar keys %compnew;
    $printbuffer .= "Different number of dihedral interactions in XML files (old, $nold; new, $nnew).\n";
   }
  }
 }else{
  ${$fail}{'DIHEDRALS EXIST'}=-1;
 }

 if($intc==$intm && $intc > 0){
  ${$fail}{'INTERACTION COUNT: DIHEDRALS'}=0;
 }
 if($inttc==$inttm && $inttc > 0){
  ${$fail}{'INTERACTION VALUES: DIHEDRALS'}=0;
 }else{
  $printbuffer .= "Not all dihedral interactions matched.\n";
 }
 return "$printbuffer";
} 

sub checkcontacts{ 
 # checks that the contacts in XML files are updated correctly.
 my ($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$conhash)=@_;
 my $printbuffer="";
 if((defined $xmlold->{'contacts'} &&  defined $xmlnew->{'contacts'}) || (! defined $xmlold->{'contacts'} && ! defined $xmlnew->{'contacts'})){
  ${$fail}{'CONTACTS EXIST'}=0;
 }else{
  $printbuffer .= "\"contacts\" present in only one xml.\n";
 }

 my $intc=0;
 my $intm=0;
 my $inttc=0;
 my $inttm=0;
 if(defined $xmlold->{'contacts'} &&  defined $xmlnew->{'contacts'}){
  # set up the index groups and parameter changes
  #...
  my $c=0;
  my $m=0;
  my %hold=%{$xmlold->{'contacts'}};
  my %hnew=%{$xmlnew->{'contacts'}};
  foreach my $type(keys %{$hold{'contacts_type'}}){
   # temp hashes that will store the selected atom groups.
   my %cg0;
   my %cg1;
   my %params;
   if(defined $conhash->{$type}){
    # found a contact type that was supposed to be changed.
    # get the parameters;
    %params=%{$conhash->{$type}->{parms}};
    my $g0=${$conhash->{$type}->{groups}}[0];
    my $g1=${$conhash->{$type}->{groups}}[1];
    %cg0=%{$atomgroup->{${$grpnms}[$g0]}};
    %cg1=%{$atomgroup->{${$grpnms}[$g1]}};
   }
   my $interactionhashold=$hold{'contacts_type'}->{$type};
   my $interactionhashnew=$hnew{'contacts_type'}->{$type};
   my %compnew;
   my %compold;
   # save new xml data
   foreach my $K(@{$interactionhashnew->{'interaction'}}){
    # each entry is a hash
    my %hash=%{$K};
    my $string="";
    foreach my $key(sort keys %hash){
     my $value;
     if($key !~ m/^[ij]$/){
      #$value=sprintf("%.2e", $hash{$key});
      $value=truncaten($hash{$key});
     }else{
      $value=$hash{$key};
     }
     $string .= "$key $value ";
    }
    $compnew{$string}=0;
   }

   # deal with original xml and see if it needs to be modified
   if(defined $conhash->{$type}){
    # update values to the expected numbers
    foreach my $K(@{$interactionhashold->{'interaction'}}){
     # each entry is a hash
     my %hash=%{$K};
     my $string="";

     my $i=$hash{"i"};
     my $j=$hash{"j"};
     my $mod=1;
     if((defined $cg0{$i} && defined $cg1{$j} ) || ( defined $cg1{$i} && defined $cg0{$j})){
      if(defined $conhash->{$type}->{"modremove"}){
       # since it should be removed, don't include it in %compold
       next;
      }else{
       # change something about this contact
       foreach my $key(sort keys %hash){
        my $value;
        if (defined $params{$key}){
         # this is parameter to update
         $value=eval("$hash{$key}$params{$key}");
         $value=truncaten($value);
        }else{
         $value=$hash{$key};
         if($key !~ m/^[ij]$/){
          $value=truncaten($value);
         }
        }
        $string .= "$key $value ";
       }
      }
     }else{
      # atoms not in groups, just save
      foreach my $key(sort keys %hash){
       my $value;
       if($key !~ m/^[ij]$/){
        $value=truncaten($hash{$key});
       }else{
        $value=$hash{$key};
       }
       $string .= "$key $value ";
      }
     }
     $compold{$string}=0;
    } 
   }else{
    # not changing this set, so just copy

    foreach my $K(@{$interactionhashold->{'interaction'}}){
     # each entry is a hash
     my %hash=%{$K};
     my $string="";
     foreach my $key(sort keys %hash){
      my $value;
      if($key !~ m/^[ij]$/){
       $value=truncaten($hash{$key});
      }else{
       $value=$hash{$key};
      }
      $string .= "$key $value ";
     }
     $compold{$string}=0;
    }
   } 

   $intc++;
   if(scalar keys %compold == scalar keys %compnew){
    $intm++;
    foreach my $I(keys %compold){
     $inttc++;
     if(defined $compnew{$I}){
      $inttm++;
     }else{
      $printbuffer .= "Interaction key \"$I\" expected, but not found in new xml\n";
     }
    }
    foreach my $I(keys %compnew){
     $inttc++;
     if(defined $compold{$I}){
      $inttm++;
     }else{
      $printbuffer .= "Interaction key \"$I\" found in new xml, but not expected\n";
     }
    }
 
   }else{
    my $nold=scalar keys %compold;
    my $nnew=scalar keys %compnew;
    $printbuffer .= "Different number of contact interactions in XML files (old, $nold; new, $nnew).\n";
   }
  }
 }else{
  ${$fail}{'CONTACTS EXIST'}=-1;
 }

 if($intc==$intm && $intc > 0){
  ${$fail}{'INTERACTION COUNT: CONTACTS'}=0;
 }
 if($inttc==$inttm && $inttc > 0){
  ${$fail}{'INTERACTION VALUES: CONTACTS'}=0;
 }else{
  $printbuffer .= "Not all contact interactions matched.\n";
 }
 return "$printbuffer";
} 

sub checkheadparams{ 
 # checks that the header is the same.
 my ($fail,$xmlold,$xmlnew)=@_;
 my $printmessage="";
 my $chead=0;
 my $mhead=0;
 my $cparam=0;
 my $mparam=0;
 my $cexp=0;
 my $mexp=0;
 foreach my $I(sort keys %{$xmlold}){
  $chead++;
  if(defined ${$xmlnew}{$I}){
   $mhead++;
  }
  if($I eq "OpenSMOGversion"){
   if($xmlnew->{$I} eq $xmlold->{$I}){
    ${$fail}{'XML OpenSMOGversion'}=0;
   }
   next;
  }
  foreach my $J(sort keys %{$xmlold->{$I}}){
   $chead++;
   if(defined ${$xmlnew->{$I}}{$J}){
    $mhead++;
   }
   foreach my $K(sort keys %{$xmlold->{$I}->{$J}}){
    $chead++;
    if(defined ${$xmlnew->{$I}->{$J}}{$K}){
     $mhead++;
    }
    if($I ne "constants" and $I ne "nonbond"){
     my @oldparams=@{$xmlold->{$I}->{$J}->{$K}->{'parameter'}};
     my @newparams=@{$xmlnew->{$I}->{$J}->{$K}->{'parameter'}};
     $cexp++;
     if($xmlnew->{$I}->{$J}->{$K}->{'expression'}->{'expr'} eq $xmlold->{$I}->{$J}->{$K}->{'expression'}->{'expr'}){
      $mexp++;
     }else{
      $printmessage .= "Expressions not the same for new and old XML files.\n";
     }
     $cparam++;
     if($#newparams == $#oldparams){
      $mparam++;
     }
     for (my $L=0;$L<=$#newparams;$L++){
      $cparam++;
      if($oldparams[$L] eq $newparams[$L]){
       $mparam++;
      }
     }
    }
   }
  }
 }
 ${$fail}{'XML TREE'}=abs($chead-$mhead);
 ${$fail}{'PARAM LISTS'}=abs($cparam-$mparam);
 ${$fail}{'EXPRESSION'}=abs($cexp-$mexp);
 return $printmessage;
}

sub checknonbonded{ 
 # checks that nonbonded sections are the same
 my ($fail,$xmlold,$xmlnew)=@_;
 my $printbuffer="";
 if((defined $xmlold->{'nonbond'} &&  defined $xmlnew->{'nonbond'}) || (! defined $xmlold->{'nonbond'} && ! defined $xmlnew->{'nonbond'})){
  ${$fail}{'NONBOND EXIST'}=0;
 }else{
  $printbuffer .= "\"nonbond\" only found in one of the xml files\n";
 }
 if(defined $xmlold->{'nonbond'} &&  defined $xmlnew->{'nonbond'}){
  my $c=0;
  my $m=0;
  my %hold=%{$xmlold->{'nonbond'}->{'nonbond_bytype'}};
  my %hnew=%{$xmlnew->{'nonbond'}->{'nonbond_bytype'}};
  if ($hnew{expression}->{'expr'} eq $hold{expression}->{'expr'}){
   ${$fail}{'EXPRESSION: NONBOND'}=0;
  }
  my @pold=@{$hold{'nonbond_param'}};
  my @pnew=@{$hnew{'nonbond_param'}};
  my $nbcount=0;
  my $nbmatch=0;
  for(my $m=0;$m<$#pold;$m++){
   my %thold=%{$pold[$m]};
   my %thnew=%{$pnew[$m]};
   foreach my $I(keys %thnew){
    $nbcount++;
    if($thnew{$I} eq $thold{$I}){
     $nbmatch++;
    }
   }
   foreach my $I(keys %thold){
    $nbcount++;
    if($thnew{$I} eq $thold{$I}){
     $nbmatch++;
    }
   }
  }
  if($nbcount==0){
   $printbuffer .= "Did not find any nonbonded parameters\n";   
  }else{
   ${$fail}{'NONBOND PARAM VALUES'}=abs($nbcount-$nbmatch);
  }
 }else{
  ${$fail}{'EXPRESSION: NONBOND'}=-1;
  ${$fail}{'NONBOND PARAM VALUES'}=-1;
 }
 return $printbuffer;
} 

sub truncaten
{
  # a likely-to-be unnecessarily complicated way of truncated scientific notion numbers
  my ($value)=@_;
  $value=sprintf("%.5e", $value);
  $value =~ s/e/ /g;
  my ($value, $places)=split(/ /,$value);
  my $factor = 10**5;
  my $v=int($value * $factor) / $factor;
  $v="${v}e$places";
  return $v;
}

return 1;
