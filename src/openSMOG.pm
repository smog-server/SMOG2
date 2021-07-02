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
	my $ref=\%{$OSref->{$type}->{$type . "_type"}->{$name}};
	$ref->{expression}->{"expr"}=$expr;

	my @params=@{$params};
	foreach my $en(@params){
		push(@{$ref->{parameter}},"$en");
	}
}

sub openSMOGfunctionExists{
	my ($OSref,$type,$name)=@_;
	if(exists $OSref->{$type}->{$type . "_type"}->{$name}){
		return 1;
	}else{
		return 0;
	}
}

sub addOShash{
	my ($OSref,$stuff)=@_;
	my @stuff=@{$stuff};
	my $ref=\%{$OSref->{$stuff[2]}->{$stuff[2] . "_type"}->{$stuff[3]}};
# @stuff is the array that contains the following information: i, j, interaction type, function name, @parameter (in the order found in $OSref->{$type}->{$name}->{parameters} array)
	my %tmphash;
	$tmphash{"i"}=$stuff[0];
	$tmphash{"j"}=$stuff[1];
	my @newstuff;
	push(@newstuff,"i");
	push(@newstuff,"$stuff[0]");
	push(@newstuff,"j");
	push(@newstuff,"$stuff[1]");
	my $parn=4;
	foreach my $param(@{$ref->{parameter}}){
		push(@newstuff,"$param");
		push(@newstuff,"$stuff[$parn]");
		$tmphash{$param}=$stuff[$parn];
		$parn++;
	}
	#push(@{$ref->{interaction}},\@newstuff);
	push(@{$ref->{interaction}},\%tmphash);

}

sub readopenSMOGxml {
	my ($XMLin)=@_;
	if(-f $XMLin){
		my $xml = new XML::Simple;
		my $data = $xml->XMLin($XMLin,KeyAttr=>{contacts_type=>"name"},ForceArray=>["contacts_type","parameter","interaction"]);
		return $data;
	}else{
		return 1;
	}
}

sub openSMOGwriteXML{
	my ($OSref,$openSMOGxml)=@_;
        # OSref is a handle to the hash holding all information to be written.
        # $openSMOGxml is the output file name
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
	# we will make a more versatile version later. We will probably replace this with an XMLlib call, later. But, in order to format it exactly the way we want, directly writing it is ok.
	my $size = keys %{$OSref};
	if($size != 0){
		my $xmlout="<openSMOGforces>\n";
		my $handle0=$OSref;

		foreach my $type(sort keys %{$handle0}){
			$xmlout .= "$space<$type>\n";
			my $handle1=$handle0->{$type};
			foreach my $subtype(sort keys %{$handle1}){
			   	my $handle2=$handle1->{$subtype};
			   	foreach my $name(sort keys %{$handle2}){
			   		my $handle3=$handle2->{"$name"};
			   	     	$xmlout .= "$twos<$subtype name=\"$name\">\n";
			   	     	my $expr=$handle3->{expression}->{"expr"};
			   	     	$xmlout .= "$threes<expression expr=\"$expr\"/>\n";
					my @paramlist=@{$handle3->{parameter}};
			   	     	foreach my $param(@paramlist){
			   	     		$xmlout .= "$threes<parameter>$param</parameter>\n";
			   	     	}

			   	     	foreach my $param(@{$handle3->{interaction}}){
			   	     		$xmlout .="$threes<interaction";
			   	     		my %tmphash=%{$param};
						foreach my $key("i","j",@paramlist){
		   	     				my $fmt;
		   	     				# write integers as integers.  Everything else as scientific notation
		   	     				if($tmphash{$key} =~ m/^[0-9]*$/){
		   	     					$fmt="%i";
		   	     				}else{
		   	     					$fmt="%7.5e";
		   	     				}
		   	     				my $val=sprintf("$fmt",$tmphash{$key});
		   	     				$xmlout .=" $key=\"$val\"";
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

		# check that the written file aligns with the intended schema.
		my $doc = XML::LibXML->new->parse_file($openSMOGxml);
		my $xmlschema = XML::LibXML::Schema->new( location => "$ENV{SMOG_PATH}/share/schemas/openSMOG.xsd", no_network => 1 );
		eval { $xmlschema->validate( $doc ); };
		smog_quit("Failed to validate $openSMOGxml against schema file $ENV{SMOG_PATH}/share/schemas/openSMOG.xsd . $@ \nThis is due to an XML formatting issue, which is probably a bug in the code.  Please send a bug report to info\@smog-server.org") if $@;

		return "\t$openSMOGxml\n";
	}
	return "";
}

1;
