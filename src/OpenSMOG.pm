#########################################################################################
#                          Structure-based Model (SMOG) software
#    This package is the product of contributions from a number of people, including:
#            Jeffrey Noel, Mariana Levi, Antonio Oliveira, Vinícius Contessoto,
#             Mohit Raghunathan, Joyce Yang, Prasad Bandarkar, Udayan Mohanty,
#                          Ailun Wang, Heiko Lammert, Ryan Hayes,
#                               Jose Onuchic & Paul Whitford
#
#          Copyright (c) 2015,2016,2018,2021,2022 The SMOG development team at
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
########################
#OpenSMOG package provides all routines used for OpenSMOG support
########################

package OpenSMOG;
use strict;
use warnings FATAL => 'all';
use SMOGglobals;
use smog_common;
use Exporter;
use XML::Simple qw(:strict);
use XML::LibXML;
our @ISA = 'Exporter';
our @EXPORT = qw(OShashAddFunction OShashAddConstants OShashAddNBFunction OpenSMOGfunctionExists checkOpenSMOGparam AddInteractionOShash AddDihedralOShash AddNonbondOShash AddSettingsOShash readOpenSMOGxml OpenSMOGwriteXML OpenSMOGextractXML newOpenSMOGfunction OpenSMOGAddNBsettoXML %fTypes %fTypesArgNum %OSrestrict);
our %fTypes;
our %fTypesArgNum;
our $OpenSMOG;
our %OpenSMOGatoms2restrain;
# make a list of names that can't be used as parameters in OpenSMOG.
our %OSrestrict;
our %NBtypespresent;
########## OpenSMOG routines
sub OShashAddFunction{
	my ($OSref,$type,$name,$expr,$params)=@_;
	if($type ne "contacts" and $type ne "dihedrals"){
		smog_quit("OpenSMOG currently only supports modified contact and dihedral potentials through \"functions\" declarations. Nonbonded custom potentials may be defined in the .nb file. Issue processing $name");
	}
	my $ref=\%{$OSref->{$type}->{$type . "_type"}->{$name}};
	$ref->{expression}->{"expr"}=$expr;

	my @params=@{$params};
	foreach my $en(@params){
		push(@{$ref->{parameter}},"$en");
	}
}

sub OShashAddNBFunction{
	my ($OSref,$interactions)=@_;
	my $expr = $interactions->{"CustomNonBonded"}->{"OpenSMOGpotential"};
	my $params=$interactions->{"CustomNonBonded"}->{"parameters"};
	my $cutoff=$interactions->{"CustomNonBonded"}->{"cutoff"};
	my $ref=\%{$OSref->{nonbond}->{nonbond_bytype}};
	$ref->{expression}->{"expr"}=$expr;
	$ref->{cutoff}=$cutoff;
	my @params=@{$params};
	foreach my $en(@params){
		push(@{$ref->{parameter}},"$en");
	}
}

sub OpenSMOGfunctionExists{
	my ($OSref,$type,$name)=@_;
	if(exists $OSref->{$type}->{$type . "_type"}->{$name}){
		return 1;
	}else{
		return 0;
	}
}

sub AddInteractionOShash{
	my ($inttype,$OSref,$stuff)=@_;
	my @stuff=@{$stuff};
	my $nameindex;
	if($inttype eq "contact"){
		$nameindex=2;
	}elsif ($inttype eq "dihedral"){
		$nameindex=4;
	}else{
		smog_quit("Internal error: OS Interactions hash issue 1");
	}

	my $ref=\%{$OSref->{$stuff[$nameindex]}->{$stuff[$nameindex] . "_type"}->{$stuff[$nameindex+1]}};
# @stuff is the array that contains the following information: i, j, interaction type, function name, @parameter (in the order found in $OSref->{$type}->{$name}->{parameters} array)
	my %tmphash;
	$tmphash{"i"}=$stuff[0];
	$tmphash{"j"}=$stuff[1];
	if ($inttype eq "dihedral"){
		$tmphash{"k"}=$stuff[2];
		$tmphash{"l"}=$stuff[3];
	}
	my $parn=$nameindex+2;
	foreach my $param(@{$ref->{parameter}}){
		$tmphash{$param}=$stuff[$parn];
		$parn++;
	}
	push(@{$ref->{interaction}},\%tmphash);
}

sub AddNonbondOShash{
	my ($OSref,$stuff)=@_;
	my @stuff=@{$stuff};
	my $ref=\%{$OSref->{nonbond}->{nonbond_bytype}};
# @stuff is the array that contains the following information: type1,type2, @parameter (in the order found in $OSref->{$type}->{$name}->{parameters} array)
	my %tmphash;
	$tmphash{"type1"}=$stuff[0];
	$tmphash{"type2"}=$stuff[1];
	my $parn=2;
	foreach my $param(@{$ref->{parameter}}){
		$tmphash{$param}=$stuff[$parn];
		$parn++;
	}
	push(@{$ref->{nonbond_param}},\%tmphash);

}

sub OShashAddConstants{
	my ($data,$OSref)=@_;

        if(defined $data->{"OpenSMOGsettings"}->[0]->{"constants"}){

		my $ref=\%{$OSref->{"constants"}};
		my $get=$data->{"OpenSMOGsettings"}->[0]->{"constants"}->[0]->{"constant"};
		for my $i (keys %{$get}){
			$ref->{"constant"}->{$i}->{"value"}=$get->{$i}->{"value"};
			$OSrestrict{$i}=0;
		}
	}
}

sub readOpenSMOGxml {
	my ($XMLin)=@_;
	if(-f $XMLin){
		my $xml = new XML::Simple;
		my $data = $xml->XMLin($XMLin,KeyAttr=>{contacts_type=>"name",dihedrals_type=>"name",constant=>"name"},ForceArray=>["contacts_type","dihedrals_type","constant","parameter","interaction","nonbond_param"]);
		return $data;
	}else{
		return 1;
	}
}

sub OpenSMOGwriteXML{
	my ($OSref,$OpenSMOGxml,$comments,$ionnm)=@_;
        # OSref is a handle to the hash holding all information to be written.
        # $OpenSMOGxml is the output file name
	my $space=" ";
	my $ones="$space";
	my $twos="$space$space";
	my $threes="$space$space$space";
	# this is a very limited XML writer that is made specifically for OpenSMOG-formatted contact hashes
	# we will make a more versatile version later. We will probably replace this with an XMLlib call, later. But, in order to format it exactly the way we want, directly writing it is ok.
	my $size = keys %{$OSref};
	if($size != 0){
		my $xmlout="<!--    OpenSMOG xml file generated by SMOG2.\n$comments -->\n";
		$xmlout .= "<OpenSMOGforces>\n";

		my $handle0=$OSref;

		foreach my $type(sort keys %{$handle0}){
			if($type eq "contacts" or $type eq "dihedrals"){
				$xmlout .= OpenSMOGwriteXMLinteractions($type,$handle0,$type,$space);
			}elsif($type eq "constants"){
				$xmlout .= OpenSMOGwriteXMLconstants($handle0,$type,$space);
			}elsif($type eq "nonbond"){
                                $xmlout .= OpenSMOGwriteXMLnonbond($handle0,$type,$space,$ionnm);
			}elsif($type eq "OpenSMOGminVersion"){
				next;
			}else{
				smog_quit("Internal Error: When writing OpenSMOG XML file, type $type not supported.");
			}
		}	
		$xmlout.="</OpenSMOGforces>\n";
		open(XMLOO,">$OpenSMOGxml") or smog_quit("Unable to open $OpenSMOGxml for writing");
		print XMLOO $xmlout;
		close(XMLOO);

		# check that the written file aligns with the intended schema.
		my $doc = XML::LibXML->new->parse_file($OpenSMOGxml);
		my $xmlschema = XML::LibXML::Schema->new( location => "$ENV{SMOG_PATH}/share/schemas/OpenSMOG.xsd", no_network => 1 );
		eval { $xmlschema->validate( $doc ); };
		smog_quit("Failed to validate $OpenSMOGxml against schema file $ENV{SMOG_PATH}/share/schemas/OpenSMOG.xsd . $@ \nThis is due to an XML formatting issue, which is probably a bug in the code.  Please send a bug report to info\@smog-server.org") if $@;

		return "\t$OpenSMOGxml\n";
	}
	return "";
}

sub OpenSMOGwriteXMLconstants{
	my ($handle0,$type,$space)=@_;
	my $ones="$space";
	my $twos="$space$space";
	my $threes="$space$space$space";
	
	my $localxmlout = "$space<$type>\n";

	my $handle1=$handle0->{$type}->{constant};
	foreach my $name(sort keys %{$handle1}){
		$localxmlout .= "$twos<constant name=\"$name\" value=\"$handle1->{$name}->{\"value\"}\"/>\n";
	}

	$localxmlout .= "$ones</$type>\n";
	my @num=split(/\s+/,$localxmlout);
	my $num = @num; 
	if($num > 3){	
		# there must be some content, so write it.
		return $localxmlout;
	}else{
		return "";
	}
}

sub OpenSMOGwriteXMLinteractions{
	my ($inttype,$handle0,$type,$space)=@_;
	my $ones="$space";
	my $twos="$space$space";
	my $threes="$space$space$space";
	my @interactingindices;
	if($inttype eq "contacts"){
		@interactingindices=("i","j");
	}elsif($inttype eq "dihedrals"){
		@interactingindices=("i","j","k","l");
	}else{
		smog_quit("Internal error: XML writing error 1");
	}
	my $localxmlout = "$space<$type>\n";
	my $handle1=$handle0->{$type};
	foreach my $subtype(sort keys %{$handle1}){
	   	my $handle2=$handle1->{$subtype};
	   	foreach my $name(sort keys %{$handle2}){
	   		my $handle3=$handle2->{"$name"};
			my $anydefined=0;
	   	     	foreach my $param(@{$handle3->{interaction}}){
				if(defined $param){
					$anydefined=1;
					last;
				}
			}
			if($anydefined==0){
				# means none of the interactions exist in the final system.  So, don't write this type
				next;
			}

	   	     	$localxmlout .= "$twos<$subtype name=\"$name\">\n";
			my $expr=$handle3->{expression}->{"expr"};
	   	     	$localxmlout .= "$threes<expression expr=\"$expr\"/>\n";
			my @paramlist=@{$handle3->{parameter}};
	   	     	foreach my $param(@paramlist){
	   	     		$localxmlout .= "$threes<parameter>$param</parameter>\n";
	   	     	}
	   	     	foreach my $param(@{$handle3->{interaction}}){
				if(!defined $param){
					#must have been deleted
					next;
				}
	   	     		$localxmlout .="$threes<interaction";
	   	     		my %tmphash=%{$param};
				foreach my $key(@interactingindices,@paramlist){
   	     				my $fmt;
   	     				# write integers as integers.  Everything else as scientific notation
   	     				if($tmphash{$key} =~ m/^[0-9]*$/){
   	     					$fmt="%i";
   	     				}else{
   	     					$fmt="%7.5e";
   	     				}
   	     				my $val=sprintf("$fmt",$tmphash{$key});
   	     				$localxmlout .=" $key=\"$val\"";
				}
	   	     		$localxmlout .="/>\n";
	   	     	}
	   	     	$localxmlout .= "$twos</$subtype>\n";
           	 }
	}
	$localxmlout .= "$ones</$type>\n";
	my @num=split(/\s+/,$localxmlout);
	my $num = @num; 
	if($num > 3){	
		# there must be some content, so write it.
		return $localxmlout;
	}else{
		return "";
	}
}

sub OpenSMOGwriteXMLnonbond{
	my ($handle0,$type,$space,$ionnm)=@_;
	my $ones="$space";
	my $twos="$space$space";
	my $threes="$space$space$space";
	
	my $localxmlout = "$space<$type>\n";

        my $handle1=$handle0->{$type};
        foreach my $subtype(sort keys %{$handle1}){
		if($subtype ne "nonbond_bytype"){
			smog_quit("Only nonbond_bytype is currently supported for non-bonded custom potentials in OpenSMOG. Found $subtype");
		}
                my $handle3=$handle1->{$subtype};
		$localxmlout .= "$twos<$subtype>\n";
		if (defined $handle3->{"cutoff"}){
	   		$localxmlout .= "$threes<cutoff distance=\"$handle3->{cutoff}\"/>\n";
		}
		my $expr=$handle3->{expression}->{"expr"};
		if(! defined $expr){
			smog_quit("Custom non-bonded potential was not found in the OpenSMOG XML file. This can happen if you are adding ions that use custom potentials to a system that does not have a custom potential defined. Check your input XML file.");
		}
		$localxmlout .= "$threes<expression expr=\"$expr\"/>\n";
		my @paramlist=();
		my $paramlistnull=0;
		if(defined $handle3->{parameter}){
			@paramlist=@{$handle3->{parameter}};
		}else{
			# this force was not defined to have any parameters. So, defining it as null will tell SMOG2 to write all possible nb pairs with no parameters, so the same potential will be used for all interactions in the model. This can only return true when generating a SMOG model with smog2.  At that time, the parameter list is taken from the templates. smog_ion, or other programs will read a generated XML file, which will have a parameter defined, even if it is just "null". If there is anything, then we will not set everything to null=0.
			$paramlistnull=1;
			push(@paramlist,"null");
		}
		foreach my $param(@paramlist){
			$localxmlout .= "$threes<parameter>$param</parameter>\n";
		}
		if($paramlistnull == 1){
			# fill with null parameter values for all combinations of atom types
			my %used;
			foreach my $I(keys %NBtypespresent){
				foreach my $J(keys %NBtypespresent){
					if(defined $used{"$J,$I"}){
						next;
					}else{
						$localxmlout .= "   <nonbond_param type1=\"$I\" type2=\"$J\" null=\"0\"/>\n";
					}

					$used{"$I,$J"}=1;
				}
			}
		}else{
			my %typesinsystem;
			my %nbpairsinsystem;
			# try to populate with parameters from extras file
			foreach my $param(@{$handle3->{nonbond_param}}){
				if(!defined $param){
					#must have been deleted
					next;
				}
				$localxmlout .="$threes<nonbond_param";
				my %tmphash=%{$param};
				my $tmptype="";
				foreach my $key("type1","type2"){
					my $fmt="%s";
					my $val=sprintf("$fmt",$tmphash{$key});
					$localxmlout .=" $key=\"$val\"";
					$typesinsystem{$val}=1;
					if($tmptype eq ""){
						$tmptype ="$val-";
					}else{
						$tmptype .="$val";
					}
				}
				# keep track of the pairs that already have nb params
				$nbpairsinsystem{$tmptype}=1;
				foreach my $key(@paramlist){
					my $fmt;
					# write integers as integers.  Everything else as scientific notation
					if(!defined $tmphash{$key}){
						smog_quit("When writing the OpenSMOG XML file, there appears to be an incorrectly assigned\nnonbond_param. Specifically, value for $key was not found when trying to write parameters for\n$tmptype interactions. This generally occurs if your extras file is not properly formatted.\nFormatting is checked when running smog2, but the file you are using may differ from\nthat used to generate the original model.");
					}
					if($tmphash{$key} =~ m/^[0-9]*$/){
						$fmt="%i";
					}else{
						$fmt="%7.5e";
					}
					my $val=sprintf("$fmt",$tmphash{$key});
					$localxmlout .=" $key=\"$val\"";
				}
				$localxmlout .="/>\n";
			}
			# now check if paramlist[0]=="null" and $ionnm is defined. If both conditions are satisfied see if we need "null" parameters between the ion and existing types in the system.
			if(defined $ionnm and $paramlist[0] eq "null"){
				$typesinsystem{"$ionnm"}=1;
				foreach my $an(keys %typesinsystem){
					unless(exists $nbpairsinsystem{"$ionnm-$an"} || exists $nbpairsinsystem{"$an-$ionnm"} ){
						$localxmlout .="$threes<nonbond_param type1=\"$an\" type2=\"$ionnm\" null=\"0\"/>\n";
					}
				}				
			}
		}
		$localxmlout .= "$twos</$subtype>\n";
	}
	$localxmlout .= "$ones</$type>\n";
	my @num=split(/\s+/,$localxmlout);
	my $num = @num; 
	if($num > 3){	
		# there must be some content, so write it.
		return $localxmlout;
	}else{
		return "";
	}
}

sub OpenSMOGextractXML{
	my ($OSref,$OpenSMOGxml,$keepatoms,$typesinsystem,$header)=@_;
        # OSref is a handle to the hash holding all information to be written.
        # $OpenSMOGxml is the output file name
	# Only load the module if we are writing an OpenSMOG file
 	my $checkPackage=`\$perl4smog -e "use XML::LibXML" 2>&1`;
        if(length($checkPackage) > 0) { smog_quit("Perl module XML::LibXML not installed. Since you are using OpenSMOG, we can not continue...")}
        # this was a workaround to a cryptic shared variable error in perl
	use if 0==0 , "XML::LibXML";
	OpenSMOGextractContacts($OSref,$keepatoms);
	OpenSMOGextractDihedrals($OSref,$keepatoms);
	OpenSMOGextractNonBonds($OSref,$keepatoms,$typesinsystem);
	OpenSMOGwriteXML($OSref,$OpenSMOGxml,$header,);
	return \%OpenSMOGatoms2restrain;
}

sub OpenSMOGextractContacts{
	my ($OSref,$keepatoms)=@_;
	my $type="contacts";
	if(defined $OSref->{$type}){
		print "Contacts found in OpenSMOG XML file.  Will extract.\n";
		my $handle1=$OSref->{$type};
		foreach my $subtype(sort keys %{$handle1}){
		   	my $handle2=$handle1->{$subtype};
		   	foreach my $name(sort keys %{$handle2}){
		   		my $handle3=$handle2->{"$name"}->{interaction};
		   	     	#foreach my $param(@{$handle3->{interaction}}){
		   	     	for (my $I=0;$I<=$#{$handle3};$I++){
					if(OpenSMOGkeepContact(${$handle3}[$I],$keepatoms)){
						delete ${$handle3}[$I];
						# if evals to 1, then delete
					}
					# this renumbers, or removes the interaction
		   	     	}
        	   	 }
		}
	}
}


sub OpenSMOGextractDihedrals{
	my ($OSref,$keepatoms)=@_;
	my $type="dihedrals";
	if(defined $OSref->{$type}){
		print "Dihedrals found in OpenSMOG XML file.  Will extract.\n";
		my $handle1=$OSref->{$type};
		foreach my $subtype(sort keys %{$handle1}){
		   	my $handle2=$handle1->{$subtype};
		   	foreach my $name(sort keys %{$handle2}){
		   		my $handle3=$handle2->{"$name"}->{interaction};
		   	     	#foreach my $param(@{$handle3->{interaction}}){
		   	     	for (my $I=0;$I<=$#{$handle3};$I++){
					if(OpenSMOGkeepDihedral(${$handle3}[$I],$keepatoms)){
						delete ${$handle3}[$I];
						# if evals to 1, then delete
					}
					# this renumbers, or removes the interaction
		   	     	}
        	   	 }
		}
	}
}


sub OpenSMOGextractNonBonds{
	my ($OSref,$keepatoms,$typesinsystem)=@_;
	my %typesinsystem=%{$typesinsystem};
	my $type="nonbond";
	if(defined $OSref->{$type}){
		print "Nonbonded terms found in OpenSMOG XML file.  Will extract.\n";
        	my $handle1=$OSref->{$type};
        	foreach my $subtype(sort keys %{$handle1}){
			if($subtype ne "nonbond_bytype"){
				smog_quit("Only nonbond_bytype is currently supported for non-bonded custom potentials in OpenSMOG. Found $subtype");
			}
        	        my $handle3=$handle1->{$subtype}->{nonbond_param};
			for (my $I=0;$I<=$#{$handle3};$I++){
				my $type1=${$handle3}[$I]->{"type1"};
				my $type2=${$handle3}[$I]->{"type2"};

				unless(exists $typesinsystem{$type1} && exists $typesinsystem{$type2}){
					delete ${$handle3}[$I];
        	   	 	}
			}
		}
	}
}

sub OpenSMOGkeepContact {
	my ($tmphash,$keepatoms)=@_;
        if(exists $keepatoms->{$tmphash->{"i"}} && exists $keepatoms->{$tmphash->{"j"}}){
		$tmphash->{"i"}=$keepatoms->{$tmphash->{"i"}};
		$tmphash->{"j"}=$keepatoms->{$tmphash->{"j"}};
		return 0;
	}elsif(exists $keepatoms->{$tmphash->{"i"}}){
		$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"i"}}}=1;
	}elsif(exists $keepatoms->{$tmphash->{"j"}}){
		$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"j"}}}=1;
	}
	return 1;
}

sub OpenSMOGkeepDihedral {
	my ($tmphash,$keepatoms)=@_;
        if(exists $keepatoms->{$tmphash->{"i"}} && exists $keepatoms->{$tmphash->{"j"}} && exists $keepatoms->{$tmphash->{"k"}} && exists $keepatoms->{$tmphash->{"l"}}){
		$tmphash->{"i"}=$keepatoms->{$tmphash->{"i"}};
		$tmphash->{"j"}=$keepatoms->{$tmphash->{"j"}};
		$tmphash->{"k"}=$keepatoms->{$tmphash->{"k"}};
		$tmphash->{"l"}=$keepatoms->{$tmphash->{"l"}};
		return 0;
	}else{
		if(exists $keepatoms->{$tmphash->{"i"}}){
			$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"i"}}}=1;
		}
		if(exists $keepatoms->{$tmphash->{"j"}}){
			$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"j"}}}=1;
		}

		if(exists $keepatoms->{$tmphash->{"k"}}){
			$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"k"}}}=1;
		}
		if(exists $keepatoms->{$tmphash->{"l"}}){
			$OpenSMOGatoms2restrain{$keepatoms->{$tmphash->{"l"}}}=1;
		}
	}
	return 1;
}

sub newOpenSMOGfunction{
	my ($OpenSMOGhandle,$fh,$fN)=@_;
	# $fh is the functions handle
	# is the new function name to add

	if($fN =~ m/[\+\-\*\/\^]/){
		smog_quit("OpenSMOG function names can not have +, -, *, ^, or /.  Problematic definition: $fN");
	}
        checkOpenSMOGparam("function",$fN);

	if ( exists $fTypes{"$fN"}){
		smog_quit("Can not create new OpenSMOG function. $fN is a protected name.");
	}
	# the fType array is the gromacs function type. But, this is not used in OpenSMOG.  We still use the array to keep track of defined functions, and use it as a tag for determining that this is an OpenSMOG potential.  But we don't use the value in the output
	$fTypes{"$fN"}="-2";

	# set the number of required parameters
	my $parmstring=$fh->{$fN}->{"OpenSMOGparameters"};
	my $parmstringorig=$parmstring;
	$parmstring =~ s/\s+//g;
	if($parmstring =~ m/^\,|\,\,|\,$/){
		smog_quit("Incorrectly formatted parameter list given for function $fN. Found \"$parmstring\"\nCheck .sif file.");
	}
	my @parmarr=split(/\,/,$parmstring);
	$fTypesArgNum{"$fN"}=$#parmarr+1;

	# even though OpenSMOG doesn't use gromacs directives, various checks in the code use the directive name for checking validity of function definitions.
	if($fh->{$fN}->{"OpenSMOGtype"} eq "contact"){
		$fh->{$fN}->{"directive"}="pairs";
		# creating this element to keep track of the fact that it was a custom term.
		$fh->{$fN}->{"IsCustom"}=1;
	}elsif($fh->{$fN}->{"OpenSMOGtype"} eq "dihedral"){
		$fh->{$fN}->{"directive"}="dihedrals";
		# creating this element to keep track of the fact that it was a custom term.
		$fh->{$fN}->{"IsCustom"}=1;
	}else{
		my $name=$fh->{$fN}->{"OpenSMOGtype"};
		smog_quit("Function $fN has OpenSMOGtype of $name, which is not supported.");
	}

	my $pot=$fh->{$fN}->{"OpenSMOGpotential"};
	checkPotentialFunction($pot,$pot,\@parmarr,$OpenSMOGhandle,0);
	my $pind=0;
	my %seenparm;
	foreach my $param(@parmarr){
		checkOpenSMOGparam($fN,$param);
		if(exists $seenparm{$param}){
			smog_quit("Function $fN defines $param as a parameter more than once. See .sif file. Found $fh->{$fN}->{\"OpenSMOGparameters\"}");
		}else{
			$seenparm{$param}=1;
		}
		if($param =~ m/^weight$/){
			#keep track of which parameter is the weight, if any.
			$OpenSMOGhandle->{$fN}->{weight}=$pind;
		}
		$pind++;
		if($pot !~ m/([^0-9a-zA-Z]|^)$param([^0-9a-zA-Z]|$)/){
			smog_quit("Function $fN defines $param as a parameter, but it does not appear in the definition of the potential. Found $fh->{$fN}->{\"OpenSMOGpotential\"}");
		}
	}

	$OpenSMOGhandle->{$fN}->{expression}=$fh->{$fN}->{"OpenSMOGpotential"}  ;
	$OpenSMOGhandle->{$fN}->{cutoff}=$fh->{$fN}->{"cutoff"};
	foreach my $par(@parmarr){
		push(@{$OpenSMOGhandle->{$fN}->{parameters}},"$par");
	}

}

sub checkOpenSMOGparam{
	my ($term,$param)=@_;
	if($param !~ m/^[a-zA-Z]([a-zA-Z0-9_]+)?$/){
		my $message;
		if($term eq "function"){
			$message="Function names must begin with a letter, and be followed by any number of letters, numbers, or underscores. Offending parameter: $param";
		}else{
			$message="Parameters with $term custom potential must begin with a letter, and be followed by any number of letters, numbers, or underscores. Offending parameter: $param";
		}
		smog_quit("OpenSMOG issue. $message");
	}
	if($term ne "function" && exists $OSrestrict{$param}){
		smog_quit("Parameters \"$param\" used with $term custom potential is a reserved name (i.e. may denote a function, charge, index, distance, etc).");
	}
}

sub OpenSMOGAddNBsettoXML{
	my ($OpenSMOG,$OpenSMOGout,$AddCustomParmsToXML,$NBbuffer,$ionnm,$header)=@_;
	# read the input XML
	my $data=readOpenSMOGxml($OpenSMOG);
	# add content
	my @lines=split(/\n/,$NBbuffer);
	foreach my $line(@lines){
        	my @tarray=split(/\s+/,$line);
		splice(@tarray,2,1);
		AddNonbondOShash($data,\@tarray);
	}
        $data->{nonbond}->{nonbond_bytype}->{cutoff}=$data->{nonbond}->{nonbond_bytype}->{cutoff}->{'distance'};

	# write new XML
	OpenSMOGwriteXML($data,$OpenSMOGout,$header,$ionnm);
}

1;
