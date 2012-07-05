# Copyright (c) 2008-2010 Martin Becker.  All rights reserved.
# This package is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.
#
# $Id: Makefile.PL 82 2010-08-31 00:13:03Z demetri $

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME           => 'Math::Polynomial',
    VERSION_FROM   => 'lib/Math/Polynomial.pm',
    PREREQ_PM      => {
        overload   => 0,
        Carp       => 0,
    },
    ABSTRACT_FROM  => 'lib/Math/Polynomial.pm',
    AUTHOR         => 'Martin Becker <becker-cpan-mp@cozap.com>',
    $ExtUtils::MakeMaker::VERSION >= 6.18? (
        SIGN       => 1,
    ): (),
    $ExtUtils::MakeMaker::VERSION >= 6.31? (
        LICENSE    => 'perl',
    ): (),
    $ExtUtils::MakeMaker::VERSION >= 6.46? (
        META_MERGE => {
            recommends => {
                'Math::Complex'       => 0,
                'Math::BigRat'        => 0,
                'Math::ModInt'        => 0,
            },
        },
    ): (),
    $ExtUtils::MakeMaker::VERSION >= 6.47? (
        MIN_PERL_VERSION => '5.006',
    ): (),
    $ExtUtils::MakeMaker::VERSION >= 6.56? (
        BUILD_REQUIRES => {
            'Test'                => 0,
            'Math::Complex'       => 0,
            'File::Spec'          => 0,
            'File::Basename'      => 0,
            'FindBin'             => 0,
            'ExtUtils::MakeMaker' => 0,
        },
    ): (),
);
