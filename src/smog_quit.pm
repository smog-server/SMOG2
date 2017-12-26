package smog_quit;
#####################
# Error call        #
# ##################

sub smog_quit
{
	my ($LINE)=@_;
	if($main::maxwarn > $main::warncount || $main::maxwarn ==-1){
		$main::warncount++;
		warn("\nWARNING $main::warncount : $LINE\n");
	}else{
		print "\n\nFATAL ERROR:  $LINE\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG2 manual, or you can email us at info\@smog-server.org. \n\nNOTE: For diagnostic purposes, you can try to ignore the error with the -warn flag.\nHowever, it is not recommended that top files generated with this flag be used for an actual simulation.\n";
		exit;
	}
}

1;
