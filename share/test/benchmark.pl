#!/usr/bin/env perl 

#template top is first argument
#top to compare is second argument
$template_top=$ARGV[0];
$check_top=$ARGV[1];

#list of headings in the .top. Will break if there is a heading that doesn't exist
@list=("defaults","atomtypes","moleculetype","atoms","bonds","angles","dihedrals","pairs","exclusions","system","molecules");
@contents=(); #dummy array to contain .top
$keywordIndex;
$differenceThreshold=0.002; #fractional difference between two floats in order to trigger an error
$EXACT_MATCH=0;
#read in template
&readtop($template_top);
@template_top=@contents;
@contents=();
#read in top to be checked;
&readtop($check_top);
@check_top=@contents;
@contents=();

$i=0;
#check the total number of lines agrees under each heading
foreach $i (0 .. $#list) {
	if( scalar @{ $check_top[$i]} != scalar @{ $template_top[$i]} ) {
		print "ERROR.\nDifferent number of entries under directive --> [ $list[$i] ]\n";
		print "Your generated file $check_top is different from the template $template_top .\n";
		exit(1); 
	}
}
#since number of lines agrees, I can assume the loops will not die
#check absolute agreement
if(EXACT_MATCH > 0) {
	foreach $i (0 .. $#list) {
		foreach $j (0 .. scalar @{$check_top[$i]}-1) {
			#print $i," ",$j," ";
			foreach $k (0 .. scalar @{$check_top[$i][$j]}-1) {
				#print $check_top[$i][$j][$k]," ";
				if($check_top[$i][$j][$k] ne $template_top[$i][$j][$k]) {
					&complain($i,$j,$k,$template_top[$i][$j][$k],$check_top[$i][$j][$k]);
				}
			}
			#print "\n";
		}
	}
} else { #allow $differenceThreshold variation in integers and floats

	foreach $i (0 .. $#list) {
		foreach $j (0 .. scalar @{$template_top[$i]}-1) {
			foreach $k (0 .. scalar @{$template_top[$i][$j]}-1) {
				if(&whatAmI($template_top[$i][$j][$k]) == 1) { #integer, must be exact
					if($check_top[$i][$j][$k] ne $template_top[$i][$j][$k]) {
						&complain($i,$j,$k,$template_top[$i][$j][$k],$check_top[$i][$j][$k]);
					}
				} elsif(&whatAmI($template_top[$i][$j][$k]) == 2) { #float, within differenceThreshold
					if($check_top[$i][$j][$k] >= 0) {
						if($check_top[$i][$j][$k]*(1+$differenceThreshold) < $template_top[$i][$j][$k] ||
						 	$check_top[$i][$j][$k]*(1-$differenceThreshold) > $template_top[$i][$j][$k]) {
								 &complain($i,$j,$k,$template_top[$i][$j][$k],$check_top[$i][$j][$k]);
						}
					} else {
						if($check_top[$i][$j][$k]*(1+$differenceThreshold) > $template_top[$i][$j][$k] ||
						 	$check_top[$i][$j][$k]*(1-$differenceThreshold) < $template_top[$i][$j][$k]) {
								 &complain($i,$j,$k,$template_top[$i][$j][$k],$check_top[$i][$j][$k]);
						}
					}
				} else { #string, must be exact
					if($check_top[$i][$j][$k] ne $template_top[$i][$j][$k]) {
						&complain($i,$j,$k,$template_top[$i][$j][$k],$check_top[$i][$j][$k]);
					}
				}
			}
		}
	}
}
sub complain {
	print "ERROR.\n";
	print "\nDifference at: Directive [ ",$list[$_[0]]," ]\n";
	print "Line: ",($_[1]+1)," Token: ",($_[2]+1)," Type: ",&whatAmItext($_[3]),"\n";
	print "Should be: ",$_[3],"\n";
	print "But is instead: ",$_[4],"\n";
}

sub whatAmI {
	if($_[0] =~ /^[+-]?[0-9]*[0-9,eE+-]*$/) {return 1;} #integer
	if($_[0] =~ /^[+-]?[0-9]*\.[0-9]*[0-9,eE+-]*$/) {return 2;} #float
	return 3; #not integer or float
}

sub whatAmItext {
	if($_[0] =~ /^[+-]?[0-9]*[0-9,eE+-]*$/) {return "integer";} #integer
	if($_[0] =~ /^[+-]?[0-9]*\.[0-9]*[0-9,eE+-]*$/) {return "float";} #float
	return "string"; #not integer or float
}

sub readtop {
	$topname=$_[0];
	open(TOP,$topname) or die " $topname can not be opened...";
	print "*******************\n";
	print "Reading .top: ",$topname,"\n";
	#loop through topfile
	$LINE=<TOP>;
	chomp($LINE);
	@tokens=split(' ',$LINE);
	while($tokens[1] ne "molecules") { #assumes [ molecules ] is last and doesn't check if correct
		if($tokens[0] ne "[") {
			$LINE=<TOP>;
			chomp($LINE);
			@tokens=split(' ',$LINE);
		} elsif( @tokens > 0 && &tokenIsKeyword($tokens[1]) ){
			print "Found heading keyword: ",$list[$keywordIndex]," ";
			$lineNum=0;
			$LINE=<TOP>;
			chomp($LINE);
			@tokens=split(' ',$LINE);
			until($tokens[0] eq "[") {
				if($LINE !~ /^;/ && @tokens != 0) { #ignore comments and blank lines
					push @{$contents[$keywordIndex][$lineNum]}, @tokens;
					$lineNum++;
				}
				$LINE=<TOP>;
				chomp($LINE);
				@tokens=split(' ',$LINE);
			}
			print "\t",scalar @{ $contents[$keywordIndex] }," entries.\n";
		} else {  die "ERROR. Unknown directive: ",$LINE," If correct please add to \@list in the code.";exit(1) }
	}
	close(TOP);		
}

sub tokenIsKeyword {
	$value = 0;$count = 0;
	foreach $a (@list){
		if($_[0] eq $a) {$value = 1;$keywordIndex=$count;}
		$count++;
	}
	return $value;
}
