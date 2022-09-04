use strict;
use warnings FATAL => 'all';
use smog_common;
use check_common;
use check_adjust;
use check_table;
use check_ions;
use check_scale;
use check_extract;

# This is the main script that runs SMOG2 tools and then checks to see if the generated files are correct.
print <<EOT;
*****************************************************************************************
                                   smog-tools-check                                   

     smog-tools-check is part of the SMOG 2 distribution, available at smog-server.org     

                  This tool will check your installation of SMOG 2 tools.

                       See the SMOG manual for usage guidelines.

            For questions regarding this script, contact info\@smog-server.org              
*****************************************************************************************
EOT

&checkForModules;
 
my $EXEC_SMOG=$ENV{'smog_exec'};
my $SMOGDIR=$ENV{'SMOG_PATH'};
my $EXEC_ADJUST=$ENV{'exec_adjust'};
my $EXEC_IONS=$ENV{'exec_ions'};
my $EXEC_EXTRACT=$ENV{'exec_extract'};
my $EXEC_SCALE=$ENV{'exec_scale'};
my $EXEC_TABLE=$ENV{'exec_table'};

our $TOLERANCE=$ENV{'TOLERANCE'};
our $MAXTHR=1.0+$TOLERANCE;
our $MINTHR=1.0-$TOLERANCE;
our $PRECISION=$ENV{'PRECISION'};

our $PDB_DIR="share/PDB.files";

quit_init();

my ($CHECKGMX,$CHECKGMXGAUSSIAN,$GMXVER,$GMXPATH,$GMXPATHGAUSSIAN,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA)=initgmxparams($SMOGDIR);
if($CHECKGMXGAUSSIAN eq "yes"){
 smog_quit("Testing gaussian potentials not supported with smog-tools-check. Please set CHECKGMXGAUSSIAN to no.");
}

my $FAILED;
my $message;
my $FAILSUM=0;
#things to check
#extract: make sure the energetics are correct.  compare to original
#	make sure the restraints are on the right atoms
#	ensure no restraints when off

my $tested=0;
my %checkthese;
my $testall=0;
if(@ARGV>0){
 $testall=1;
 foreach my $name(@ARGV){
  $checkthese{$name}=0;
 }
}
if(defined $checkthese{"ions"} || @ARGV==0){
 print "\nTesting smog_ions\n";
 ($FAILED,$message)=check_ions($EXEC_IONS,$PDB_DIR,$CHECKGMX,$GMXVER,$GMXPATH,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA);
 if($FAILED >0){$FAILSUM++};
 $tested++;
}
if(defined $checkthese{"extract"} || @ARGV==0){
 print "\nTesting smog_extract\n";
 ($FAILED,$message)=check_extract($EXEC_EXTRACT,$PDB_DIR,$CHECKGMX,$GMXVER,$GMXPATH,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA);
 if($FAILED >0){$FAILSUM++};
 $tested++;
}
if(defined $checkthese{"tablegen"} || @ARGV==0){
 print "\nTesting smog_tablegen\n";
 ($FAILED,$message)=check_table($EXEC_TABLE,$PDB_DIR);
 if($FAILED >0){$FAILSUM++};
 $tested++;
}
if(defined $checkthese{"adjustPDB"} || @ARGV==0){
 print "\nTesting smog_adjustPDB\n";
 ($FAILED,$message)=check_adjust($EXEC_ADJUST,$EXEC_SMOG,"share");
 if($FAILED >0){$FAILSUM++};
 $tested++;
}
if(defined $checkthese{"scale-energies"} || @ARGV==0){
 print "\nTesting smog_scale-energies\n";
 ($FAILED,$message)=check_scale($EXEC_SCALE,$PDB_DIR,$CHECKGMX,$GMXVER,$GMXPATH,$GMXEXEC,$GMXEDITCONF,$GMXMDP,$GMXMDPCA);
 if($FAILED >0){$FAILSUM++};
 $tested++;
}
if($FAILSUM>0){
 print "\n\nSOME TESTS FAILED.  SEE EARLIER MESSAGES\n\n";	
 exit (1);
}elsif($tested==0){
 print "\n\nNo tests performed... Possible options: adjustPDB, tablegen, scale-energies, extract, ions\n\n"; 
 exit (1);
}elsif($testall ==0){
 print "\n\nPassed all SMOG tool checks!\n\n"; 
 exit (0);
}else{
 print "\n\nPassed the following checks: "; 
 foreach my $test(keys %checkthese){
  print "$test ";
 }
 print "\n\n";
 exit (0);
}

