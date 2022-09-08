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
# SMOGglobals package provides things that are used all over the place
########################

package SMOGglobals;
use strict;
use warnings FATAL => 'all';
use smog_common;
use Exporter;
use XML::Simple qw(:strict);
use XML::LibXML;
our @ISA = 'Exporter';
our @EXPORT = qw(initOSrestrict %OSrestrict %NBtypespresent);
# make a list of names that can't be used as parameters in OpenSMOG.
our %OSrestrict;
our %NBtypespresent;

sub initOSrestrict {
	foreach my $i ("q1", "q2", "r", "r_c", "i", "j", "k", "l", "m", "n", "type1", "type2", "sqrt", "exp", "log", "sin", "cos", "sec", "csc", "tan", "cot", "asin", "acos", "atan", "sinh", "cosh", "tanh", "erf", "erfc", "min", "max", "abs", "floor", "ceil", "step", "delta", "select", "null"){ $OSrestrict{$i}=0;}
}
1;
