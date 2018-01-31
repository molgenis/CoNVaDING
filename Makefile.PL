use strict;
use warnings;

require 5.22.0;

use ExtUtils::MakeMaker;
WriteMakefile
(
	NAME         => 'CoNVaDING',
	VERSION_FROM => 'CoNVaDING.pl',
	PREREQ_PM    => {
		'File::Glob' => '1.24',
		'File::Basename' => '2.85',
		'Getopt::Long' => '2.45',
                'List::Util' => '1.41',
                'Math::Complex' => '1.59',
                'POSIX' => '1.53',
<<<<<<< HEAD
                'Statistics::Distributions' => '1.02',
=======
                'Statistics::Distributions' => '1.02'
>>>>>>> f91e9cdef69dc5c63aefaedc2372e16b85274b54
                'Statistics::Normality' => '0.01',
                'File::Temp' => '0.2304',
	},
	BUILD_REQUIRES => {
		'Test::More' => '0.47'
	},
	EXE_FILES => [
		'CoNVaDING.pl'
	],
);
