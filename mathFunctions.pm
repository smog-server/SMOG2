package mathFunctions;
use strict;
use warnings;

my $pi = 3.14159265358979;
sub pi {return $pi;}
sub deg_to_rad { ($_[0]/180) * $pi }
sub rad_to_deg { ($_[0]/$pi) * 180 }
sub asin { atan2($_[0], sqrt(1 - $_[0] * $_[0])); }
sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ); }
sub tan  { sin($_[0]) / cos($_[0]);  }
sub atan { atan2($_[0],1); }