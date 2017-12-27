package smog_common;
use strict;
use Exporter;

#####################
# Init error vars   #
#####################
our $maxwarn;
our $warncount;
our @ISA = 'Exporter';
our @EXPORT = qw($warncount $maxwarn quit_init smog_quit warnsummary warninfo checkForModules);


#####################
# Error routiness   #
# ###################

sub quit_init
{
	$maxwarn=0;
	$warncount=0;
}

sub smog_quit
{
	my ($LINE)=@_;
		$warncount++;
	if($maxwarn > $warncount || $maxwarn ==-1){
		warn("\nWARNING $warncount : $LINE\n");
	}elsif($maxwarn == $warncount && $maxwarn>0){
		print "\nWARNING $warncount : $LINE\n";
		warn("\n\nFATAL: REACHED USER-DEFINED MAXIMUM NUMBER OF WARNINGS. QUITTING.\n\n");
		exit;
	}else{
		print "\n\nFATAL ERROR:  $LINE\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG2 manual, or you can email us at info\@smog-server.org. \n\nNOTE: For diagnostic purposes, you can try to ignore the error with the -warn flag.\nHowever, it is not recommended that output obtained with this flag be used for an actual simulation.\n";
		exit;
	}
}

sub warninfo
{
	if($maxwarn > 0 || $maxwarn ==-1 ){
		print "\n\n-warn $maxwarn selected. Errors will be reported as warnings.\nBe cautious, ignoring errors can lead to unpredictable results. \n ONLY use this option if you are sure the error is harmless.\n\n";
	}
}

sub warnsummary
{
	if ($warncount == 1){
		print "\n\n NOTE: There was $warncount warning. It is recommended that you read all warnings carefully.\n\n"; 
	}elsif ($warncount > 1){
		print "\n\n NOTE: There were $warncount warnings. It is recommended that you read all warnings carefully.\n\n"; 
	}
}

sub checkForModules {
	my $checkPackage; my $sum=0;
	$checkPackage=`echo \$perl4smog | wc | awk '{print \$3}'`;
	if($checkPackage < 2) { print "\nSMOG 2 failed to launch.\n\nEnvironment variable perl4smog not set, maybe you need to edit the configure.smog2 script and run it with \"source configure.smog2\"\n"; $sum++;}else{
		$checkPackage=`\$perl4smog -e "use XML::Simple" 2>&1 | wc -l | awk '{print \$1}'`;
		if($checkPackage > 0) { print "Perl module XML::Simple not installed!\n"; $sum++;}
		$checkPackage=`\$perl4smog -e "use XML::Validator::Schema" 2>&1 | wc -l | awk '{print \$1}'`;
		if($checkPackage > 0) { print "Perl module XML::Validator::Schema not installed!\n"; $sum++;}
		$checkPackage=`\$perl4smog -e "use Exporter" 2>&1 | wc -l | awk '{print \$1}'`;
		if($checkPackage > 0) { print "Perl module Exporter not installed!\n"; $sum++;}
		$checkPackage=`\$perl4smog -e "use String::Util" 2>&1 | wc -l | awk '{print \$1}'`;
		if($checkPackage > 0) { print "Perl module String::Util not installed!\n"; $sum++;}
		$checkPackage=`\$perl4smog -e "use PDL" 2>&1 | wc -l | awk '{print \$1}'`;
		if($checkPackage > 0) { print "Perl Data Language not installed!\n"; $sum++;}
	}
	$checkPackage=`which java | wc -l | awk '{print \$1}'`;
	if($checkPackage < 1) { print "Java might not be installed. This package assumes Java 1.7 or greater is in the path as 'java'.\n"; $sum++;}
	if($sum > 0) { print "Need above packages before smog-check (and smog2) can run. Some hints may be in the SMOG 2 manual.\n"; exit(1); }
}



1;
