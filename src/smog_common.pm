#########################################################################################
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#            Jeffrey Noel, Mariana Levi, Antonio Oliveira, Vinícius Contessoto,
#             Mohit Raghunathan, Joyce Yang, Prasad Bandarkar, Udayan Mohanty,
#                          Ailun Wang, Heiko Lammert, Ryan Hayes,
#                               Jose Onuchic & Paul Whitford
#
#            Copyright (c) 2015,2016,2018,2021, The SMOG development team at
#                        Rice University and Northeastern University
#
#          SMOG 2, Shadow and OpenSMOG are available at http://smog-server.org
#
#          You can direct questions to info@smog-server.org, or the smog-users forum,
#          which you can find at https://mailman.rice.edu/mailman/listinfo/smog-users
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#########################################################################################

#############################
# smog_common has many routines that are used by smog2 and smog-check
#############################

package smog_common;
use strict;
use warnings FATAL => 'all';
use Exporter;
use XML::Simple qw(:strict);
use OpenSMOG;
#####################
# Init error vars   #
#####################
our $VERSION="2.5beta";
our $maxwarn;
our $warncount;
our $allwarncount;
our $notecount;
our @convarray;
our %reverthash;
our $BaseN;
our @ISA = 'Exporter';
our @EXPORT = qw($allwarncount $warncount $maxwarn note_init smog_note quit_init smog_quit warnsummary warninfo checkForModules checkcomment hascontent loadfile checkdirectives %supported_directives checkforinclude readindexfile printdashed printcenter checksuffix checkalreadyexists InitLargeBase BaseTentoLarge BaseLargetoTen printhostdate whatAmI trim evalsub validateXML checkPotentialFunction GetCustomParms $VERSION);
our %supported_directives;
#####################
# Error routiness   #
# ###################

sub quit_init
{
	$maxwarn=0;
	$warncount=0;
	$allwarncount=0;
}

sub smog_quit
{
	my ($LINE,$warn)=@_;
	$allwarncount++;
	if(defined $warn){
		#if $warn is defined, it means we that this call is never treated as fatal
		print "\nWARNING $allwarncount (non-fatal warning) : $LINE\n\n";
	}elsif($maxwarn > $warncount || $maxwarn ==-1 ){
		$warncount++;
		print "\nWARNING $allwarncount (suppressed fatal error no. $warncount): $LINE\n\n";
	}elsif($maxwarn <= $warncount && $maxwarn>0){
		print "\nWARNING $allwarncount : $LINE\n\n";
		print "\n\nEXCEEDED USER-DEFINED MAXIMUM NUMBER OF WARNINGS. QUITTING.\n\n";
		exit(1);
	}else{
		print "\n\nFATAL ERROR:  $LINE\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG2 manual, or you can email us at info\@smog-server.org. \n\nNOTE: For diagnostic purposes, you can try to ignore the error with the -warn flag.\nHowever, it is not recommended that output obtained with this flag be used for an actual simulation.\n";
		exit(1);
	}
}

sub note_init
{
	$notecount=0;
}

sub smog_note
{
	my ($LINE)=@_;
	$notecount++;
	print "\nNOTE $notecount: $LINE\n\n";
}

sub warninfo
{
	if($maxwarn > 0 || $maxwarn ==-1 ){
		print "\n\n-warn $maxwarn selected. Errors will be reported as warnings.\nBe cautious, ignoring errors can lead to unpredictable results. \n ONLY use this option if you are sure the error is harmless.\n\n";
	}
}

sub warnsummary
{

	if ($notecount == 1){
		print "\n\nTHERE WAS $notecount NOTE. Even if you know what you are doing, it is always worth checking any notes.\n\n"; 
	}elsif ($notecount > 1){
		print "\n\nTHERE WERE $notecount NOTES. Even if you know what you are doing, it is always worth checking any notes.\n\n"; 
	}

	if ($allwarncount == 1){
		print "\n\n NOTE: There was $allwarncount warning. It is recommended that you read all warnings carefully.\n\n"; 
	}elsif ($allwarncount > 1){
		print "\n\n NOTE: There were $allwarncount warnings. It is recommended that you read all warnings carefully.\n\n"; 
	}
}

################
# module check
################

sub checkForModules {
	my $checkPackage; my $sum=0;
	$checkPackage=`which java `;
	if(length($checkPackage) < 1) { 
		smog_quit("Java not found. SMOG 2 requires Java 1.7 or greater to be in the path as 'java'.\n")
	}
}

##################
# PARSING TOOL
##################
sub checkcomment
{
	my ($LINE) = @_;
	chomp($LINE);
	my $comment;
	if($LINE =~ m/(;.*)/){
		$comment=$1;
	}else{
		$comment="";
	}
	# remove comments
	$LINE =~ s/;.*$//g; 
	$LINE =~ s/\t/ /g; 
	$LINE = trim($LINE);
	$LINE =~ s/ +/ /g;
	if( $LINE =~ m/[#!\^\$]/ ){
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
		# add the \n since we just stripped them using split
		$DATA[$I] .="\n";
		my $string1 = $DATA[$I];
		open my($fh), "<", \$string1 or smog_quit("internal error 1") ; # reading from the data in $string
		my $first_line = <$fh>; 
		close $fh;
		my ($line,$comment)=checkcomment($first_line);
	        $line =~ s/\]/ \]/g;
	        #$line =~ s/\t+/ /g;
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
			smog_quit("Directive \"$keys\" not found in input .top file.",0);
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
		$LINE = trim($LINE);
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

		if(!defined $groupname){
			smog_quit("Index file doesn\'t have any group names listed.");
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
        	$line = trim($line);
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

sub InitLargeBase {
        # this is only done once.
        @convarray = (0..9,"a".."z","A".."Z");
	$BaseN=$#convarray+1;
        my $J=0;
        foreach my $I(@convarray){
                $reverthash{$I}=$J;
                $J++;
        }
}

sub BaseTentoLarge {
        my ($base10,$digits)=@_;
        my $baselarge="";

        my $val=$base10;
        for (my $I=0;$I<$digits;$I++){
                my $rem=$val % $BaseN;
                $val=int($val/$BaseN);
                $rem=$convarray[$rem];
                $baselarge = $rem . $baselarge;
        }

	if($val != 0) {

		my $N=scalar(@convarray)-1;
		$N=$convarray[$N];
		my $maxV="";
        	for (my $I=0;$I<$digits;$I++){
			$maxV .= "$N";
		}
		my $max10=BaseLargetoTen($maxV);
		smog_quit("Using base-$BaseN numbering for atoms/residues and hit maximum defined value of $maxV ($max10)");

	}

        return $baselarge;
}

sub BaseLargetoTen {
        my ($baselarge)=@_;
	# remove any space
	$baselarge =~ s/\s+//g;
	my $length=length($baselarge);

        my $base10=0;
	my $power=0;
        for (my $I=$length-1;$I>=0;$I--){
		my $char=substr($baselarge,$I,1);
		$base10+=$reverthash{$char}*$BaseN**$power;
		$power++;
        }

        return $base10;
}

sub printhostdate {

	my $date=`date`;
	my $hostname=`hostname`;
	chomp($date);
	chomp($hostname);
	my $string = "; date: $date\n; hostname: $hostname\n";
	return $string;
}

sub whatAmI {
	if($_[0] =~ /^[+-]?[0-9]+[0-9,eE+-]*$/) {return 1;} #integer
	# there is certainly a more compact way of writing thie regex.  oh well, I'll come back to it...
	if($_[0] =~ /^[+-]?[0-9]*\.[0-9]*[eE]?$/ ) {return 2;} #float
	return 3; #not integer or float
}

#removes white space from the beginning and from the end of input string 
#essentially String::Util trim()
sub trim {
	my $string = $_[0];
	$string =~ s/^\s+|\s+$//g;
	return $string;
}

sub evalsub
{
	my ($expression,$value)=@_;
	$expression =~ s/\?/$value/g;
	$expression = eval($expression);
	return $expression;
}

sub validateXML
{
	my ($file,$type) = @_;

	my $validator = XML::Validator::Schema->new(file => "$ENV{SMOG_PATH}/share/schemas/$type.xsd");
	my $parser = XML::SAX::ParserFactory->parser(Handler => $validator);
	eval { $parser->parse_uri($file) };
 smog_quit("Failed at validating $file: $@ \nThis is due to an XML formatting issue.  The most common issue is that an element is missing a tagline. For example, something like this may appear in your file,
  <child>
    <subchild>.....
  </child>

whereas the following would be appropriate:
  <child>
    <subchild>.....

    </subchild>
  </child>") if $@;

}

sub checkPotentialFunction{
	my ($hd,$parms,$OSref)=@_;
	my $func=$hd->{"OpenSMOGpotential"};
	$func =~ s/\s+//g;
	if($func =~ m/^(.*?)?([^\s+a-zA-Z0-9_\;\,\=\/\-\+\*\^\(\)])(.*)?$/ ){
		# regex explained: ^ start, (.*?)? non-greedy and optional any character, (not any letter, number, or allowed operator), (.*)? optional any character, $ end.
		my $pos=length($1)+1;
		smog_quit("Unsupported characters found in OpenSMOG energy function:\n$func\nEnergy functions may only have the following characters: a-z A-Z 0-9 _ ; = , / - + * ^ ( )\nCharacter $pos is a \"$2\".");
	}
	if($func =~ m/\(\)/){
		smog_quit("Empty closed parentheses found in OpenSMOG energy function:\n$func");
	}
	if($func =~ /\*\*/){
		smog_quit("Unsupported use of ** (only ^ allowed) in OpenSMOG energy function:\n$func");
	}
	my @array = split(//, $func);

	my $count=0;
	my $charn=0;
	foreach my $char(@array){
		if($char =~ m/\(/){
			$count++;
		}
		if($char =~ m/\)/){
			$count--;
		}
		if($count < 0){
			my $badstring="";
			for (my $I=$charn-3;$I<$charn+3;$I++){
				if($I>=0){
					$badstring .= $array[$I];
				}
			}
			smog_quit("Function used in templates has unbalanced parentheses. Closed parentheses found before opening. Problem found around \"$badstring\" in the following definition:\n$func\n");
		}
		$charn++;
	}
	if($count > 0){
		smog_quit("Function used in templates has unbalanced parentheses (more open than closed). Problematic definition:\n$func\n");
	}

        # since it at least has ok parentheses, let's make sure it only uses allowed parameters and functions

	if($func =~ /\*\*/){
		smog_quit("Unsupported use of \"**\" in OpenSMOG energy function:\n$func\nA ^ character must be used.");
	}
	my %ref;
	if(defined $OSref->{"constants"}->{"constant"}){
        	%ref=%{$OSref->{"constants"}->{"constant"}};
	}
	$func =~ s/\s+//g;
	my @par=split(/[\(\)\*\/\-\+\,\^]+/,$func);
	my %phash;
	for my $I(@{$parms}){
		$phash{$I}=0;
	}
	foreach my $M(@par){
		unless(whatAmI($M) < 3 || exists $OSrestrict{$M} || exists $phash{$M} || defined $ref{$M}){
			smog_quit("\"$M\" identified as a function, but it is not found in OpenSMOG energy function:\n$func");
		}
	}
}

sub GetCustomParms{
	# $data is a hashref that contains sif information imported with XMLin 
	my ($data)=@_;
	if(exists $data->{"CustomNonBonded"}){
		my @interHandle = @{$data->{"CustomNonBonded"}};
		if(defined $interHandle[0]->{"OpenSMOGparameters"} && $interHandle[0]->{"parameters"}){
			if($interHandle[0]->{"OpenSMOGparameters"} ne $interHandle[0]->{"parameters"}){
				smog_quit("When using CustomNonbonded, give either the parameters or OpenSMOGparameters child element in the .sif file.  If both are listed, they must be identical.");
			}
		}elsif(defined $interHandle[0]->{"OpenSMOGparameters"}){
			$interHandle[0]->{"parameters"}=$interHandle[0]->{"OpenSMOGparameters"};
		}elsif(!defined $interHandle[0]->{"parameters"}){
			smog_quit("In sif file, when using CustomNonBonded, you must give either the parameters of OpenSMOGparameters child element.");
		}
		# parameters required for any custom potential definition
		# set the number of required parameters
		my $parmstring=$interHandle[0]->{"parameters"};
		my $parmstringorig=$parmstring;
		$parmstring =~ s/\s+//g;
		if($parmstring =~ m/^\,|\,\,|\,$/){
			smog_quit("Incorrectly formatted parameter list given for nonbonded CustomNonBonded. Found \"$parmstringorig\"\nCheck .nb file.");
		}
		my @parmarr=split(/\,/,$parmstring);
		return (1,\@parmarr);
	}else{
		return (0,"");
	}
}


1;
