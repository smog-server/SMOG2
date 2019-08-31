package smog_common;
use strict;
use Exporter;

#####################
# Init error vars   #
#####################
our $maxwarn;
our $warncount;
our @ISA = 'Exporter';
our @EXPORT = qw($warncount $maxwarn quit_init smog_quit warnsummary warninfo checkForModules checkcomment hascontent loadfile checkdirectives %supported_directives checkforinclude readindexfile printdashed printcenter checksuffix checkalreadyexists);
our %supported_directives;

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
	my ($LINE,$warn)=@_;
		$warncount++;
	if($maxwarn >= $warncount || $maxwarn ==-1 || defined $warn){
		warn("\nWARNING $warncount : $LINE\n\n");
	}elsif($maxwarn < $warncount && $maxwarn>0){
		print "\nWARNING $warncount : $LINE\n\n";
		warn("\n\nEXCEEDED USER-DEFINED MAXIMUM NUMBER OF WARNINGS. QUITTING.\n\n");
		exit(1);
	}else{
		print "\n\nFATAL ERROR:  $LINE\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG2 manual, or you can email us at info\@smog-server.org. \n\nNOTE: For diagnostic purposes, you can try to ignore the error with the -warn flag.\nHowever, it is not recommended that output obtained with this flag be used for an actual simulation.\n";
		exit(1);
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

################
# module check
################


sub checkForModules {
	my $checkPackage; my $sum=0;
	$checkPackage=`echo \$perl4smog | wc | awk '{print \$3}'`;
	if($checkPackage < 2) { print "\nFailed to launch.\n\nEnvironment variable perl4smog not set, maybe you need to edit the configure.smog2 script and run it with \"source configure.smog2\"\n"; $sum++;}else{
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
	if($sum > 0) { smog_quit("Need above packages before smog2, smog-check and smog tools can run. Some hints may be in the SMOG 2 manual."); }
}

##################
# PARSING TOOL
##################
sub checkcomment
{
	my ($LINE) = @_;
	my $comment;
	if($LINE =~ m/(;.*)/){
		$comment=$1;
	}else{
		$comment="";
	}
	# remove comments
	$LINE =~ s/;.*$//g; 
	$LINE =~ s/\t/ /g; 
	$LINE =~ s/^\s+|\s+$//g;
	$LINE =~ s/ +/ /g;
	if( $LINE =~ m/^#/ ){
		print "preprocessing line found in .top file\n";
	}elsif( $LINE =~ m/[#!\^\$]/ ){
		smog_quit("Special characters not recognized\n  Offending line: $LINE\n");
	}
	return ($LINE,$comment);
}

sub hascontent
{
	my ($LINE) = @_;
	# remove comments
	$LINE =~ s/;.*$//g; 
	# remove spaces and tabs
	$LINE =~ s/\s|\t//g;
	if( $LINE =~ m/[#!\^\$]/ ){
		smog_quit("Special characters not recognized in input file\n  Offending line: $LINE\n");
	}
	if($LINE eq ""){
		return 0;
	}else{
		return 1;
	}
}

## reading routines
sub loadfile
{
	my ($file)=@_;
	open(FILE, "$file") or smog_quit("Can not open $file"); 
	my $string = "";
	while (<FILE>){
		my $LINE = $_;
		chomp($LINE);
		# remove blank lines
		unless($LINE =~ m/^[\s+|\t+]$|^$/ ){ 
			if($LINE =~ m/\[\S/){
				$LINE =~ s/\[/\[ /g;
 			}
			if($LINE =~ m/\S\]/){
				$LINE =~ s/\]/ \]/g;
 			}
			$string .= "$LINE\n";
		}
	}
	close(FILE);
	return $string;
}

sub checkdirectives
{
	my ($string) = @_;
# split the top file and check that only supported directives are included.
	my %DIRLIST;
	my @DATA=split(/\n\s+\[|\n\[|^\s+\[|^\[/,$string);
	for (my $I=1;$I<=$#DATA;$I++){
		my $string1 = $DATA[$I];
		open my($fh), "<", \$string1 or smog_quit("internal error 1") ; # reading from the data in $string
		my $first_line = <$fh>; 
		close $fh;
		my ($line,$comment)=checkcomment($first_line);
	        $line =~ s/\]/ \]/g;
	        $line =~ s/\t+/ /g;
	        $line =~ s/\s+/ /g;
		my @B=split(/ /,$line);
		my $DIR=$B[0];
		chomp($DIR);
		if(!defined $B[1]){
			smog_quit("Format error near directive \"$DIR\"");
		}
		chomp($B[1]);
		if($B[1] ne "]"){
			smog_quit("Format error near directive \"$DIR\"");
		}
		if(!exists $supported_directives{$DIR}){
			smog_quit("Directive \"$DIR\" not recognized.");
		}elsif(exists $DIRLIST{$DIR}){
			smog_quit("Directive \"$DIR\" defined more than once. Currently not supported.");
		}else{
			$DIRLIST{$DIR}=$I;
		}
	}

	for my $keys(keys %supported_directives){
		if($supported_directives{$keys} == 1 && !exists $DIRLIST{$keys}){
			smog_quit("Directive \"$keys\" not found in input .top file.");
		}
		if($supported_directives{$keys} == 0 && exists $DIRLIST{$keys}){
			smog_quit("Directive \"$keys\" not supported by script.");
		}
	
	}

	return (\@DATA,\%DIRLIST);
}

sub checkforinclude
{
	my ($line,$data,$handle)=@_;
	if($data =~ m/^#/){
		smog_quit("#include file listed in .top file.  Can not process associated file (not supported). If you would like smog-tools to simply copy this include line, then use the -warn option.  Offending line:\n$line\n\n");
		print $handle "$line\n";
		return 1;
	}
	return 0;
}

sub readindexfile
{
	my ($indexFile)=@_;
	my $groupname;
	my @grpnms;
	my %groupnames;
	my $Ngrps=0;
	my %atomgroup;
	my $groupindex;
	open(ATOMLIST,"$indexFile") or smog_quit("Can\'t open $indexFile.");
	print "Reading index file $indexFile\n";
	while(<ATOMLIST>){
		my $LINEt=$_;
		chomp($LINEt);
		# remove comments first
	        my ($LINE,$comment)=checkcomment($LINEt);
		# in case we have a group directive w/o spaces
		$LINE =~ s/\[/\[ /g;
		$LINE =~ s/\]/ \]/g;
		$LINE =~ s/^\s+|\s+$//g;
		$LINE =~ s/\s+|\t+/ /g;
		my @A=split(/\s+/,$LINE);
		if($#A == -1){
			# blank line
			next;
		}
		if($A[0] eq "[" and $A[2] eq "]"){
			# must be a new group
			$groupindex=0;
			$groupname=$A[1];
			$grpnms[$Ngrps]=$groupname;
			if(exists $A[3]){
				smog_quit("Group name declarations must not have trailing characters.  $A[3] appears after $A[0] $A[1] $A[2]\n;")	
			}
			if(exists $groupnames{$groupname}){
				smog_quit("Group name $groupname declared more than once.");
			}
			$groupnames{$groupname}=$Ngrps;
			$Ngrps++;
	
			next; # this is a new group, so go to next line.
		}
	
		for(my $I=0;$I<=$#A;$I++){
			unless($A[$I] =~ m/^\d+$/){
				smog_quit("Non-numerical value for atom number in index file: $A[$I]");
			}
			if(exists $atomgroup{$groupname}{$A[$I]}){
				smog_quit("Duplicate atom $A[$I] in group $groupname");
			}else{
				$atomgroup{$groupname}{$A[$I]}=$groupindex;
				$groupindex++;
			}
		}
	}
	close(ATOMLIST);
	return ($Ngrps,\@grpnms,\%groupnames,\%atomgroup);
}

sub checksuffix
{
	my ($name,$suf)=@_;
	my ($ext)= $name =~ /(\.[^.]+)$/;
	if(!defined $ext || $ext ne "$suf"){
	        $name .=  "$suf";
	}
	return $name;
}

sub printdashed
{
	my ($headw)=@_;
	for(my $I=0;$I<$headw;$I++){
  		print "*";
	}
	print "\n";
}

sub checkalreadyexists
{
	my ($filen)=@_;
	# check if a file exists and back it up, if so.
	my $maxbu=10;
	my ($ext) = $filen =~ /(\.[^.]+)$/;
	if($filen ne "" && -e $filen){
		for(my $bu=1;$bu<=$maxbu;$bu++){
			my $buname="$filen.bu$bu";
			if( ! -e $buname){	
			print "$filen already exists.  Backing up to $buname\n";
			system("mv $filen $buname");
			last;
			}
			if($bu == $maxbu){
		 	smog_quit ("Already backed up $maxbu copies of $filen."); 
			}
		}
	}
}

sub printcenter
{
	my ($headw,$text)=@_;
	my @textarray=split(/\n/,$text);

	foreach my $line(@textarray)
	{
        	$line =~ s/^\s+|\s+$//g;
		my $L=length($line);
		my $pad=int(($headw-$L)/2);
		if($pad < 0)
		{
			print "$line\n";
		}else{
			my $string=(" " x ($pad));
			print "$string$line\n";
		}		

	}
}


1;
