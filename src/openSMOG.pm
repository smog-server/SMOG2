package openSMOG;
use strict;
use warnings FATAL => 'all';
use Exporter;
use XML::Simple qw(:strict);

#####################
# Init error vars   #
#####################
our @ISA = 'Exporter';
our @EXPORT = qw(OShashAddFunction openSMOGfunctionExists addOShash readopenSMOGxml openSMOGwriteXML);

########## openSMOG routines
sub OShashAddFunction{
	my ($OSref,$type,$name,$expr,$params)=@_;
	if($type ne "contacts"){
		smog_quit("openSMOG currently only supports modified contact potentials.  Issue processing $name");
	}
	my $ref=\%{$OSref->{$type}->[0]->{$type . "_type"}->[0]->{$name}->[0]};
	$ref->{expression}->[0]->{"expr"}=$expr;

	my @params=@{$params};
	foreach my $en(@params){
		push(@{$ref->{parameters}},"$en");
	}
}

sub openSMOGfunctionExists{
	my ($OSref,$type,$name)=@_;
	if(exists $OSref->{$type}->[0]->{$type . "_type"}->[0]->{$name}){
		return 1;
	}else{
		return 0;
	}
}

sub addOShash{
	my ($OSref,$stuff)=@_;
	my @stuff=@{$stuff};
	my $ref=\%{$OSref->{$stuff[2]}->[0]->{$stuff[2] . "_type"}->[0]->{$stuff[3]}->[0]};
# @stuff is the array that contains the following information: i, j, interaction type, function name, @parameters (in the order found in $OSref->{$type}->[0]->{$name}->[0]->{parameters} array)
	my %tmphash;
	$tmphash{"i"}=$stuff[0];
	$tmphash{"j"}=$stuff[1];
	my @newstuff;
	push(@newstuff,"i");
	push(@newstuff,"$stuff[0]");
	push(@newstuff,"j");
	push(@newstuff,"$stuff[1]");
	my $parn=4;
	foreach my $param(@{$ref->{parameters}}){
		push(@newstuff,"$param");
		push(@newstuff,"$stuff[$parn]");
		$tmphash{$param}=$stuff[$parn];
		$parn++;
	}
	push(@{$ref->{interaction}},\@newstuff);

}

sub readopenSMOGxml {
	my ($XMLin)=@_;
	if(-f $XMLin){
		my $xml = new XML::Simple;
		my $data = $xml->XMLin($XMLin,KeyAttr=>{contacts_type=>"name"},ForceArray=>1);
		return $data;
	}else{
		return 1;
	}
}

sub openSMOGwriteXML{
	my ($OSref,$openSMOGxml)=@_;
        # the arg is a handle to the hash holding all information to be written.
	# Only load the module if we are writing an openSMOG file
 	my $checkPackage=`\$perl4smog -e "use XML::LibXML" 2>&1`;
        if(length($checkPackage) > 0) { smog_quit("Perl module XML::LibXML not installed. Since you are using openSMOG, we can not continue...")}
        # this was a workaround to a cryptic shared variable error in perl
	use if 0==0 , "XML::LibXML";
	my $space=" ";
	my $ones="$space";
	my $twos="$space$space";
	my $threes="$space$space$space";
	# this is a very limited XML writer that is made specifically for openSMOG-formatted contact hashes
	# we will make a more versatile version later.
	my $size = keys %{$OSref};
	if($size != 0){
		my $xmlout="<openSMOGforces>\n";
		my $handle0=$OSref;

		foreach my $type(keys %{$handle0}){
			$xmlout .= "$space<$type>\n";
			my $handle1=$handle0->{$type}->[0];
			foreach my $subtype(keys %{$handle1}){
			   	my $handle2=$handle1->{$subtype}->[0];
			   	foreach my $name(keys %{$handle2}){
			   		my $handle3=$handle2->{$name}->[0];
			   	     	$xmlout .= "$twos<$subtype name=\"$name\">\n";
			   	     	my $expr=$handle3->{expression}->[0]->{"expr"};
			   	     	$xmlout .= "$threes<expression expr=\"$expr\" />\n";
			   	     	foreach my $param(@{$handle3->{parameters}}){
			   	     		$xmlout .= "$threes<parameter>$param</parameter>\n";
			   	     	}
			   	     	foreach my $param(@{$handle3->{interaction}}){
			   	     		$xmlout .="$threes<interaction ";
			   	     		my @tmparr=@{$param};
			   	     		for(my $I=0;$I<$#tmparr;$I++){
			   	     			my $fmt;
			   	     			# write integers as integers.  Everything else as scientific notation
			   	     			if($tmparr[$I+1] =~ m/^[0-9]*$/){
			   	     				$fmt="%i";
			   	     			}else{
			   	     				$fmt="%7.5e";
			   	     			}
			   	     			$tmparr[$I+1]=sprintf("$fmt",$tmparr[$I+1]);
			   	     			$xmlout .="$tmparr[$I]=\"$tmparr[$I+1]\" ";
			   	     			$I++;
			   	     		}
			   	     		$xmlout .="/>\n";
			   	     	}

			   	     	$xmlout .= "$twos</$subtype>\n";
                           	 }
			}

			$xmlout .= "$ones</$type>\n";
		}	
		$xmlout.="</openSMOGforces>\n";
	open(XMLOO,">$openSMOGxml") or smog_quit("Unable to open $openSMOGxml for writing");
	print XMLOO $xmlout;
	close(XMLOO);

	# check that the written file aligned with the intended schema.
	my $doc = XML::LibXML->new->parse_file($openSMOGxml);
	my $xmlschema = XML::LibXML::Schema->new( location => "$ENV{SMOG_PATH}/share/schemas/openSMOG.xsd", no_network => 1 );
	eval { $xmlschema->validate( $doc ); };
 smog_quit("Failed to validate $openSMOGxml against schema file $ENV{SMOG_PATH}/share/schemas/openSMOG.xsd . $@ \nThis is due to an XML formatting issue, which is probably a bug in the code.  Please send a bug report to info\@smog-server.org") if $@;

	return "\t$openSMOGxml\n";
	}
	return "";
}

1;
