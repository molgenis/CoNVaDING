# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl CoNVaDING.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More tests => 3;

BEGIN { use_ok('CoNVaDING') };
#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

ok(system('bash -c "echo > /dev/stderr && set -ex && perl -wc ./CoNVaDING.pl  "') == 0, "CoNVaDing.pl syntax test");

ok(system('bash -c "echo > /dev/stderr &&
        set -ex &&
        perl -w ./CoNVaDING.pl \
                -mode StartWithMatchScore \
                -bed ./Test_dataset/bedfile/Test_dataset_bedfile.bed \
                -inputDir ./Test_dataset/sample \
                -controlsDir ./Test_dataset/controls \
                -outputDir ./results/"') == 0 , "CoNVaDing.pl StartWithMatchScore test");
#perl ./CoNVaDING.pl \
#  -mode StartWithMatchScore \
#  -inputDir $DATADIR/results/StartWithAvgCount \
#  -controlsDir $DATADIR/controls \
#  -outputDir $DATADIR/results/StartWithMatchScore
#  CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
#  DATADIR="/PATH/TO/Test_dataset/"

