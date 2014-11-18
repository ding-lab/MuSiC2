use strict;
use warnings;

use Test::Most;
use TGI::MuSiC2::Complicated;

TGI::MuSiC2::Complicated::print_stuff();
ok(1, '1 is definitely okay');

done_testing;

