# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl miRTRAP.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 1;
BEGIN { use_ok('miRTRAP') };

#########################

# Currently, there are no tests for miRTRAP. This will be added in the next update.
#
