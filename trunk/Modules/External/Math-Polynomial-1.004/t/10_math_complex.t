# Copyright (c) 2007-2009 Martin Becker.  All rights reserved.
# This package is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.
#
# $Id: 10_math_complex.t 36 2009-06-08 11:51:03Z demetri $

# Checking coefficient space compatibility with Math::Complex.

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl t/10_math_complex.t'

#########################

use strict;
use warnings;
use Test;
BEGIN { plan tests => 2 };
use Math::Complex;
use Math::Polynomial 1.000;
ok(1);  # 1

#########################

my $c0 = Math::Complex->new(0, 3);
my $c1 = Math::Complex->new(2, 1);
my $c2 = Math::Complex->new(1, -2);
my $x  = Math::Complex->new(-1, 1);
my $y  = Math::Complex->new(-7, 2);
my $p = Math::Polynomial->new($c0, $c1, $c2);
ok($p->evaluate($x) == $y);     # 2

__END__
