package smog_quit;
use strict;
use Exporter;
use smog_quit;

#####################
# Init error vars   #
#####################
our $maxwarn=0;
our $warncount=0;
our @ISA = 'Exporter';
our @EXPORT = qw($warncount $maxwarn smog_quit warnsummary warninfo);


#####################
# Error routiness   #
# ###################

sub smog_quit
{
	my ($LINE)=@_;
		$warncount++;
	if($maxwarn > $warncount || $maxwarn ==-1){
		warn("\nWARNING $warncount : $LINE\n");
	}elsif($maxwarn == $warncount && $maxwarn>0){
		print "\nWARNING $warncount : $LINE\n";
		warn("\n\nREACHED $maxwarn WARNINGS. QUITTING.\n\n");
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
	if ($warncount != 0){
		print "\n\n NOTE: There were $warncount warnings. It is recommended that you read all warnings carefully.\n\n"; 
	}
}


1;
