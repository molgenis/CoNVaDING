use strict;
use warnings;

require 5.22.0;

use ExtUtils::MakeMaker;
#How to add a new file
# -create the script and put in into src
# -update the tests to do al least a syntax test
# -Update pipeline/util with newer version numer
#git add & commit it


my $extdeps;
$extdeps = [
	{	"dep" => "htslib",
		"test" => "(exit 0)",
		"hint" => "wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && tar -jxf htslib* && cd htslib* && make" 
	},{	"dep" => "samtools",
		"test" => "which samtools",
		"hint" => "wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar -jxf samtools* && cd samtools* && make"
	},{	"dep" => "bedtools",
		"test" => "bedtools --version",
		"hint" => "wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && tar -xzf https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && cd bed* && make"
	}
];

for my $dep (@{$extdeps}){
	warn "testing for ".$dep -> {'dep'}.".\n";
	CmdRunner($dep -> {'test'});
}


sub CmdRunner {
        my $ret;
        my $cmd = join(" ",@_);

        warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";

        @{$ret} = `($cmd )2>&1`;
	my $exitcode = $?;
        if ($exitcode == -1) {
                die localtime( time() ). " [ERROR] failed to execute: $!\n";
        }
	elsif ($exitcode & 127) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }
	elsif ($exitcode != 0) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }else {
               	warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $exitcode >> 8;
        }
		return @{$ret};
}

WriteMakefile(
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
