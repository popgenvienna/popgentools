# Copyright (c) 2007-2009 Martin Becker.  All rights reserved.
# This package is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.
#
# $Id: 91_pod_cover.t 36 2009-06-08 11:51:03Z demetri $

# Check for POD coverage.
# This is a test for the distribution maintainer.

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl t/91_pod_cover.t'

use strict;
use warnings;
use lib 't/lib';
use Test::MyUtils;

BEGIN {
    maintainer_only();
    use_or_bail('Test::Pod::Coverage', '0.08');
    require Test::Pod::Coverage;        # redundant, but clue for CPANTS
}

all_pod_coverage_ok();

__END__
