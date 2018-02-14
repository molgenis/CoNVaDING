use strict;
use warnings;

require 5.22.0;

use ExtUtils::MakeMaker;
WriteMakefile
(
	"NAME"         => 'CoNVaDING',
	"VERSION_FROM" => 'CoNVaDING.pl',
	"ABSTRACT" => 'Copy Number Variation Detection In Next-generation sequencing Gene panels',
	"LICENSE" => 'lgpl_3_0',
	"AUTHOR" => "Johansson LF, van Dijk F, de Boer EN, van Dijk-Bos KK, Jongbloed JD, 
van der Hout AH, Westers H, Sinke RJ, Swertz MA, Sijmons RH, Sikkema-Raddatz B.
CoNVaDING: Single Exon Variation Detection in Targeted NGS Data. Hum Mutat. 2016 May;37(5):457-64. doi: 10.1002/humu.22969. Epub 2016 Feb 24.",
	"PREREQ_PM"    => {
		'File::Glob' => '1.24',
		'File::Basename' => '2.85',
		'Getopt::Long' => '2.45',
                'List::Util' => '1.41',
                'Math::Complex' => '1.59',
                'POSIX' => '1.53',
                'Statistics::Distributions' => '1.02',
                'Statistics::Normality' => '0.01',
                'File::Temp' => '0.2304',
	},
	"BUILD_REQUIRES" => {
		'Test::More' => '0.47',
		'ExtUtils::MakeMaker' => '0'
	},
	"EXE_FILES" => [
		'CoNVaDING.pl'
	],
);
