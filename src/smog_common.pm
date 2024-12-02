#########################################################################################
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#            Jeffrey Noel, Mariana Levi, Antonio Oliveira, Vinícius Contessoto,
#             Mohit Raghunathan, Joyce Yang, Prasad Bandarkar, Udayan Mohanty,
#                          Ailun Wang, Heiko Lammert, Ryan Hayes,
#                               Jose Onuchic & Paul Whitford
#
#        Copyright (c) 2015,2016,2018,2021,2022,2024 The SMOG development team at
#                      The Center for Theoretical Biological Physics
#                       Rice University and Northeastern University
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
use SMOGglobals;
#####################
# Init error vars   #
#####################
our $VERSION="2.5-beta";
our $maxwarn;
our $warncount;
our $allwarncount;
our $notecount;
our @convarray;
our %reverthash;
our $BaseN;
our %OSrestrict;
our @ISA = 'Exporter';
our @EXPORT = qw($allwarncount $warncount $maxwarn note_init smog_note quit_init smog_quit warnsummary warninfo checkForModules checkcomment hascontent loadfile checkdirectives %supported_directives checkforinclude readindexfile printdashed printcenter checksuffix checkalreadyexists InitLargeBase BaseTentoLarge BaseLargetoTen printhostdate whatAmI trim evalsub  validateXML checkPotentialFunction GetCustomParms getgitver selectgroup listgroups getXYZfromLine selectTemplates findIonDefs findMapFile findSif $VERSION );
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
		print "\n\nFATAL ERROR:  $LINE\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG 2 manual, or you can email us at info\@smog-server.org. \n\nNOTE: For diagnostic purposes, you can try to ignore the error with the -warn flag.\nHowever, it is not recommended that output obtained with this flag be used for an actual simulation.\n";
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
		print "\n\nNOTE: There was $allwarncount warning. It is recommended that you read all warnings carefully.\n\n"; 
	}elsif ($allwarncount > 1){
		print "\n\nNOTE: There were $allwarncount warnings. It is recommended that you read all warnings carefully.\n\n"; 
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

sub getXYZfromLine
{
	my ($line,$freecoor)=@_;
	my $x;
	my $y;
	my $z;
	if(defined $freecoor){
		# Read the PDB coordinates as free-format.
		my $string=trim(substr($line, 30));
		my @coor=split(/\s+/,$string);
		$x=$coor[0];
		$y=$coor[1];
		$z=$coor[2];
		if(whatAmI($x) > 2 || whatAmI($y) > 2 || whatAmI($z) > 2){
			smog_quit("Coordinate read, but does not appear to be a number. Since you are using free-formatted coordinates (with the -freecoor flag), perhaps you are missing a delimiter. Issue found at line:\n$line");
		}
	
	}else{
		# this is also done in adjustPDB. It is also here, in case someone skips that step and makes a mistake
		if(substr($line,34,1) !~  m/\./ ) {
			smog_quit("X coordinate in PDB file is not properly formatted.  The decimal should be column 35. Problematic line:\n$line");
		}
		if(substr($line,42,1) !~  m/\./ ) {
			smog_quit("Y coordinate in PDB file is not properly formatted.  The decimal should be column 43. Problematic line:\n$line");
		}
		if(substr($line,50,1) !~  m/\./ ) {
			smog_quit("Z coordinate in PDB file is not properly formatted.  The decimal should be column 51. Problematic line:\n$line");
		}
		$x = trim(substr($line, 30, 8));
		$y = trim(substr($line, 38, 8));
		$z = trim(substr($line, 46, 8));
		if(whatAmI($x) > 2 || whatAmI($y) > 2 || whatAmI($z) > 2){
			smog_quit("Coordinate read, but does not appear to be a number. Perhaps you are using free-formatted coordinates and should employ the -freecoor flag. Issue found at line:\n$line");
		}
	}

	return ($x,$y,$z);
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
		smog_quit("#include file listed in .top file.  Can not process associated file (not supported). If you would like smog-tools to copy and save this include line, then use the -warn option.  Offending line:\n$line\n\n");
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
	my $groupindex=-1;
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

			if($groupindex == 0){
				smog_quit("Index file has an empty group ($groupname). This is probably a mistake.");
			}
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
	$expression =~ s/\?/\($value\)/g;
	$expression = eval($expression);
	return $expression;
}

sub getgitver
{
	if(defined $ENV{SMOG_COMMIT}){
		return $ENV{SMOG_COMMIT};
	}else{
		return "";
	}
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
	my ($fullfunc,$func,$parms,$OSref,$c)=@_;
	chomp($func);
	if($func =~ m/^\s+$/){
		return;
	}
	$func =~ s/\;(\s+)?$//g;
	my %phash;
 	my %definedexpressions;
	for my $I(@{$parms}){
		$phash{$I}=0;
	}
	my %ref;
	if(defined $OSref->{"constants"}->{"constant"}){
        	%ref=%{$OSref->{"constants"}->{"constant"}};
	}

	$func =~ s/\s+//g;

	# get all elements of the function, including sub-definitions
	my @funclist;
	if($func =~ m/\;/){
		@funclist=split(/\;/,$func);
	} else {
		$funclist[0]=$func;
	}

	if ($c == 0){
		#my $func=$hd->{"OpenSMOGpotential"};
		if($func =~ m/^(.*?)?([^\s+a-zA-Z0-9_\;\.\,\=\/\-\+\*\^\(\)])(.*)?$/ ){
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
		if($func =~ /\+\+/){
			smog_quit("Unsupported use of ++ in OpenSMOG energy function:\n$func");
		}
		if($func =~ /\-\-/){
			smog_quit("Unsupported use of -- in OpenSMOG energy function:\n$func");
		}
		if($func =~ /\/\//){
			smog_quit("Unsupported use of // in OpenSMOG energy function:\n$func");
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
				my $pos = (" " x $charn);
				$pos .= "^";
				smog_quit("OpenSMOG function used in templates has unbalanced parentheses. Closed parentheses found before opening. Problem found at:\n$func\n$pos");
			}
			$charn++;
		}
		if($count > 0){
			smog_quit("OpenSMOG function used in templates has unbalanced parentheses (more open than closed). Problematic definition:\n$func\n");
		}
		# get the names of subsequently-defined expressions
		if ($funclist[0] =~ m/=/) {
			smog_quit("OpenSMOG energy function can not have an equal sign in the first element. Functions:\n\t$func\n\tProblem detected with: $funclist[0]");
		}
		my %tmpfncs;
		for (my $I=1;$I<=$#funclist;$I++){
			if($funclist[$I] !~ m/^(\w+)=/) {
				smog_quit("When defining subfunctions for an OpenSMOG energy function, each must begin with a string and then an equal sign. Issue detected in element $I ($funclist[$I]) of the following energy function:\n\t$func");
			}
			if (exists $OSrestrict{$1}) {
				smog_quit("When defining subfunctions for an OpenSMOG energy function, it appears that you are trying to redefine a default OpenMM function. Function \"$1\" encountered in the following energy function:\n\t$func");
			}
			if (exists $phash{$1}) {
				smog_quit("When defining subfunctions for an OpenSMOG energy function, it appears that you are trying to redefine a parameter. \"$1=\" encountered in the following energy function:\n\t$func");
			}
			if (exists $ref{$1}) {
				smog_quit("When defining subfunctions for an OpenSMOG energy function, it appears that you are trying to redefine a constant. \"$1=\" encountered in the following energy function:\n\t$func");
			}
			if (exists $tmpfncs{$1}) {
				smog_quit("When defining subfunctions for an OpenSMOG energy function, it appears that you are trying to define \"$1\" more than once. Issue with energy function:\n\t$func");
			}
			$tmpfncs{$1}=0;
		}
	}

	for (my $I=1;$I<=$#funclist;$I++){
		$funclist[$I] =~ m/^(\w+)=/;
		$definedexpressions{$1}=0;
	}
        # since it at least has ok parentheses, let's make sure it only uses allowed parameters and functions
	$funclist[0] =~ s/^(\w+)=//g;
	my @par=split(/[\(\)\*\/\-\+\,\^]+/,$funclist[0]);
	foreach my $M(@par){
		unless(whatAmI($M) < 3 || exists $OSrestrict{$M} || exists $phash{$M} || defined $ref{$M} || exists $definedexpressions{$M} || $M eq ""){
			smog_quit("OpenSMOG energy function appears to be expressed using the quantity \"$M\". However, this does not appear to correspond to a default OpenMM function, parameter, constant, or subsequently-defined expression. Issue with energy function:\n$fullfunc",1);
		}
	}

	if($func =~ m/\;/){
		$func =~ s/(.*\;)?//;
		$c++;
		checkPotentialFunction($fullfunc,$func,$parms,$OSref,$c);
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

sub listgroups
{
	my ($Ngrps,$grpnms)=@_;
	my @grpnms=@{$grpnms};
	print "index\t:\tgroup name\n";
	for(my $I=0;$I<$Ngrps;$I++){
		print "$I\t:\t$grpnms[$I]\n";
	}
}

sub selectgroup
{
	my ($maxN)=@_;
	my $tmp=<STDIN>;
	chomp($tmp);
	$tmp =~ s/^\s+//g;
	$tmp =~ s/\s+$//g;
	unless($tmp =~ m/^\d+$/){
		smog_quit("$tmp is an invalid selection");
	}
	if($tmp <0 or $tmp >=$maxN){
		smog_quit("selection must be positive and less than or equal to $maxN");
	}
	return $tmp;
}

sub selectTemplates
{
	my ($checkforfiles)=@_;
	my $filename;
	# if an arg is given, then check if ion defs are provided.  Otherwise, just list all templates
	my $SMOGLIST;
	if(defined $ENV{'SMOG_FFDIR'}){
		$SMOGLIST=$ENV{'SMOG_FFDIR'} ;
	}else{
		$SMOGLIST=$ENV{'SMOG_PATH'} . "/share/templates/";
	}
	print "Will check for templates located in $SMOGLIST\n";
        print "Please select a force field from the list below:\n";
	print "FF Number - name : description\n";
	open(FLIST,"$SMOGLIST/ff.info") or smog_quit("Unable to open force field library file $SMOGLIST/ff.info");
	my @FLISTA;
	my $FNUM=0;
	my $line;
	my @filenames;
	while(<FLIST>){
		$line = $_;
		chomp($line);
		$line =~ s/^\s+//g;
		$line =~ s/\#.*$//g;
		unless($line =~ m/^#/ || $line eq ""){
			my @A=split(/:/,$line);
			my $dc=1;
			my $dn;
			if(defined $checkforfiles){
				if($checkforfiles eq "ions"){
					($dc,$dn)=findIonDefs("$SMOGLIST/$A[0]");
				}elsif($checkforfiles eq "adjustPDB"){
					($dc,$dn)=findMapFile("$SMOGLIST/$A[0]");
					$filenames[$FNUM]=$dn;
				}else{
					smog_quit("Internal Error: filefile");
				}
			}
			if($dc>1){
				smog_quit("Issue with template library. Found more than one set of ion definitions in $A[0]");
			}elsif($dc==1){
				$FLISTA[$FNUM]=$A[0];
				$FLISTA[$FNUM] =~ s/^\s+//g;
				$FLISTA[$FNUM] =~ s/\s+$//g;
				print "$FNUM - $line\n";

				$FNUM++;
			}
		}
	}
	$FNUM--;
	if($FNUM ==-1){
		if(defined $checkforfiles){
			if($checkforfiles eq "ions"){
				smog_quit("No templates in the library have ions.def files");
			}elsif($checkforfiles eq "adjustPDB"){
				smog_quit("No templates in the library have a mapping file");
			}	
		}else{
			smog_quit("No template directories found in the library.");
		}
	}
	my $FFN;
	my $FFv=-1;
	until($FFv == 1){
		print "Please select a force field, by FF Number:";
		$FFN=<STDIN>;
		chomp($FFN);
		$FFv=whatAmI($FFN);
		if($FFv==1){
			if($FFN<0 || $FFN > $FNUM){
				print "Invalid force field number. Please provide a value between 0 and $FNUM\n";
				$FFv=-1;
			}
		}else{
			print "Invalid force field selection. Must provide an integer between 0 and $FNUM.\n";
		}
	}
	return ($SMOGLIST,$FLISTA[$FFN],$filenames[$FFN]);
}

sub findIonDefs
{
	my ($folderName)=@_;
	my $defsexists=0;
	my $defsfile;
	my @B=split(/\s+/,$folderName);
	$folderName=$B[0];
        opendir(my $folder,$folderName) or smog_quit("Unable to open template directory listed in ff.info: $folderName");
	while(my $file = readdir($folder)){
		if($file =~ m/\.ions\.def$/ || $file =~ m/^ions\.def$/) {
			$defsexists++;
			$defsfile = "$folderName/$file";
			next;
		}
	}
	if($defsexists >1){
		smog_quit ("Found multiple ion definition files in directory $folderName");
	}
	return ($defsexists,$defsfile);

}

sub findMapFile
{
	my ($folderName)=@_;
	my $defsexists=0;
	my $defsfile;
	my @B=split(/\s+/,$folderName);
	$folderName=$B[0];
        opendir(my $folder,$folderName) or smog_quit("Unable to open template directory listed in ff.info: $folderName");
	while(my $file = readdir($folder)){
		if($file =~ m/\.map$/) {
			$defsexists++;
			$defsfile = "$folderName/$file";
			next;
		}
	}
	if($defsexists >1){
		smog_quit ("Found multiple mapping file in directory $folderName");
	}
	return ($defsexists,$defsfile);

}

sub findSif
{
	my ($folderName)=@_;
	my $sifexists=0;
	my $siffile;
        opendir(my $folder,$folderName) or smog_quit("Unable to open template directory : $folderName");
	while(my $file = readdir($folder)){
		if($file =~ m/\.sif$/) {
			$sifexists++;
			$siffile = "$folderName/$file";
			next;
		}
	}
	if($sifexists ==0){
		smog_quit ("No sif file found in directory $folderName");
	}
	if($sifexists >1){
		smog_quit ("Found multiple sif files in directory $folderName");
	}
	return $siffile;

}




1;
