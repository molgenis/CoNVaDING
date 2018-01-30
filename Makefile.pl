use strict;
use warnings;

use 5.22;

use ExtUtils::MakeMaker;
WriteMakefile
(
	NAME         => 'CoNVaDING_reload',
	VERSION_FROM => 'lib/CoNVaDING.pm',
	PREREQ_PM    => {
		'File::Glob' => '1.24',
		'File::Basename' => '2.85',
		'Getopt::Long' => '2.45',
                'List::Util' => '1.41',
                'Math::Complex' => '1.59',
                'POSIX' => '1.53_01',
                'Statistics::Normality' => '0.01',
                'File::Temp' => '0.2304',
	},
	BUILD_REQUIRES => {
		'Test::More' => '0.47'
	},
	EXE_FILES' => [
		'CoNVaDING.pl'
	],
);
