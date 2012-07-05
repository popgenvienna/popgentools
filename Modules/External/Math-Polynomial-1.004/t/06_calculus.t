# Copyright (c) 2007-2009 Martin Becker.  All rights reserved.
# This package is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.
#
# $Id: 06_calculus.t 36 2009-06-08 11:51:03Z demetri $

# Checking calculus operators.

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl t/06_calculus.t'

#########################

use strict;
use warnings;
use Test;
BEGIN { plan tests => 5 };
use Math::Polynomial 1.000;
ok(1);

#########################

sub has_coeff {
    my $p = shift;
    if (!ref($p) || !$p->isa('Math::Polynomial')) {
        print
            '# expected Math::Polynomial object, got ',
            ref($p)? ref($p): defined($p)? qq{"$p"}: 'undef', "\n";
        return 0;
    }
    my @coeff = $p->coeff;
    if (@coeff != @_ || grep {$coeff[$_] != $_[$_]} 0..$#coeff) {
        print
            '# expected coefficients (',
            join(', ', @_), '), got (', join(', ', @coeff), ")\n";
        return 0;
    }
    return 1;
}

my $p = Math::Polynomial->new(0, -5, -5, 5, 5, 1);
my $pd = $p->differentiate;
ok(has_coeff($pd, -5, -10, 15, 20, 5));

my $pi = $pd->integrate;
ok(has_coeff($pi, 0, -5, -5, 5, 5, 1));

$pi = $pd->integrate(-1);
ok(has_coeff($pi, -1, -5, -5, 5, 5, 1));

my $a = $pd->definite_integral(0, 1);
ok($a == 1);

__END__
