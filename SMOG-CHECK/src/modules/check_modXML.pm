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
 my @FAILLIST = ('NON-ZERO EXIT','XML TREE','CONSTANTS','CONSTANTS EXIST','CONTACTS EXIST','PARAM LISTS','EXPRESSION','INTERACTION COUNT: CONTACTS','INTERACTION VALUES: CONTACTS');
 my %TESTED;
 my $TESTED=\%TESTED;
# generate an AA model RNA 
 `smog2 -i $pdbdir/tRNA.pdb -AA -dname AA.tmp -OpenSMOG > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_modifyXML.");
 }else{
  clearfiles("output.smog");
 }

 print "\tChecking smog_modifyXML: test 1\n";
 my $tmpbuffer="";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $indexfile="share/PDB.files/sample.AA.ndx";
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
  savefailed(1,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml"));
 }

 print "\tChecking smog_modifyXML: test 2\n";
 my $tmpbuffer="";
 &testsperformed($TESTED,\%FAIL);
 %FAIL=resettests(\%FAIL,\@FAILLIST);
 my $indexfile="share/PDB.files/sample.AA.ndx";
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
  savefailed(1,("output.$tool","AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top","AA.tmp.out.xml"));
  print "$printbuffer\nAdditional Messages\n$tmpbuffer\n";
  foreach my $file("AA.tmp.contacts" , "AA.tmp.gro","AA.tmp.ndx", "AA.tmp.top", "AA.tmp.xml"){
   `mv tmp/$file .`;
  }
  `rmdir tmp`;
 }else{
  clearfiles(("output.$tool","AA.tmp.out.xml"));
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
    $I++;
    $settings .= "  $A[$I]  \n";
    $name=$A[$I];
   }elsif($A[$I] eq "grp"){
    $last=$A[$I];
    $I++;
    $settings .= "  $A[$I]  \n";
    push(@grparr,$A[$I]);
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

 $pbuffer.=checkcontacts($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$conhash);

 # make sure elements are identical for all non-contact and non-dihedral listings

 # go through all contact and dihedral subtypes
	# if the subtype was not supposed to be modified, make sure the two sets are identical
		# for this, make temp hashes that will hold all interactions, with unique keys - this will probably be "i-j-k...-param1-param2..."  The param order will be based on ordered names.  Compare in both directions to ensure they are identical
	# if the subtype was supposed to be changed, then save to the hash with projected/modified values.  Then, compare this against what was made by modifyXML.  The two should be identical.  
 return $pbuffer;
}

sub checkconstants{ 
 # checks that the overall structure of the XML files is the same.
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
     if($hold{$I} eq $hnew{$I}){
      $m++;
     }else{
      $printbuffer .= "Parameters $I has different values in the xml files.\n";
     }
    }else{
     $printbuffer .= "Parameters $I not found in new XML file.\n";
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

sub checkcontacts{ 
 # checks that the overall structure of the XML files is the same.
 my ($fail,$xmlold,$xmlnew,$atomgroup,$grpnms,$conhash)=@_;
 my $printbuffer="";
 if((defined $xmlold->{'contacts'} &&  defined $xmlnew->{'contacts'}) || (! defined $xmlold->{'constants'} && ! defined $xmlnew->{'constants'})){
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
      $value=sprintf("%.4e", $hash{$key});
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
      # change something about this contact
      foreach my $key(sort keys %hash){
       my $value;
       if (defined $params{$key}){
        # this is parameter to update
        $value=eval("$hash{$key}$params{$key}");
        $value=sprintf("%.4e", $value);
       }else{
        $value=$hash{$key};
       }
       $string .= "$key $value ";
      }
     }else{
      # atoms not in groups, just save
      foreach my $key(sort keys %hash){
       my $value;
       if($key !~ m/^[ij]$/){
        $value=sprintf("%.4e", $hash{$key});
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
       $value=sprintf("%.4e", $hash{$key});
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
 # checks that the overall structure of the XML files is the same.
 my ($fail,$xmlold,$xmlnew)=@_;
 my $printmessage="";
 my $chead=0;
 my $mhead=0;
 my $cparam=0;
 my $mparam=0;
 foreach my $I(sort keys %{$xmlold}){
  $chead++;
  if(defined ${$xmlnew}{$I}){
   $mhead++;
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
    my @oldparams=@{$xmlold->{$I}->{$J}->{$K}->{'parameter'}};
    my @newparams=@{$xmlnew->{$I}->{$J}->{$K}->{'parameter'}};
    if($xmlnew->{$I}->{$J}->{$K}->{'expression'}->{'expr'} eq $xmlold->{$I}->{$J}->{$K}->{'expression'}->{'expr'}){
     ${$fail}{'EXPRESSION'}=0;
    }else{
     $printmessage .= "Expression not the same for new and old XML files.\n";
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
 ${$fail}{'XML TREE'}=abs($chead-$mhead);
 ${$fail}{'PARAM LISTS'}=abs($cparam-$mparam);
 return $printmessage;
}

return 1;
