# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl CoNVaDING.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More tests => 2;

BEGIN { use_ok('CoNVaDING') };
#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

ok(system('set -ex && perl -wc ./CoNVaDING.pl 1>/dev/null') == 0, "CoNVaDing.pl test");

