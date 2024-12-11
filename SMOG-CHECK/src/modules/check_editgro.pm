package check_editgro;
use strict;
use Exporter;
use smog_common;
use check_common;
our @ISA = 'Exporter';
our @EXPORT = qw(check_editgro);

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



sub check_editgro
{
 my ($exec,$pdbdir)=@_;
 my $NFAIL=0;
 my $MESSAGE="";
 my %FAIL;
 my $FAILED;
 my $FAILSUM=0;
 my $FATAL;
 my $UNINIT;
 my $printbuffer;
 my $tool="editgro";
 my @FAILLIST = ('NON-ZERO EXIT','IDENTICAL OUTPUT');
 my %TESTED;
 my $TESTED=\%TESTED;
 %FAIL=resettests(\%FAIL,\@FAILLIST);

# generate an AA model RNA 
 `smog2 -i $pdbdir/2FP4-GDP.pdb -AA -dname AA.tmp > output.smog`;
 unless($? == 0){
  internal_error("SMOG 2 crashed.  Fix SMOG 2 before testing smog_editgro.");
 }else{
  clearfiles("output.smog");
 }
 open(SETIN,"$pdbdir/../settings/egsettings.in") or smog_quit("unable to find editgro settings file");
 my $ind=0;
 while(<SETIN>){
  my $command=$_;
  chomp($command);
  $ind++;
  print "\tChecking smog_editgro: test $ind\n";
  `$exec $command -g AA.tmp.gro -og AA.tmp.eg.gro  &> output.$tool`;
  $FAIL{"NON-ZERO EXIT"}=$?;
  if($FAIL{"NON-ZERO EXIT"} == 0){
   $FAIL{"IDENTICAL OUTPUT"}=filediff("AA.tmp.eg.gro","$pdbdir/../editgrorefs/eg.$ind.gro",1);
  }
  ($FAILED,$printbuffer)=failsum(\%FAIL,\@FAILLIST);
  $FAILSUM += $FAILED;
  if($FAILED !=0){
   savefailed($ind,("output.$tool","AA.tmp.eg.gro"));
   print "$printbuffer\n";
  }else{
   clearfiles(("output.$tool","AA.tmp.eg.gro"));
  }
 } 

 $FAILSUM+=checkalltested(\@FAILLIST,\%FAIL);

 return ($FAILSUM, $printbuffer);

}

return 1;
