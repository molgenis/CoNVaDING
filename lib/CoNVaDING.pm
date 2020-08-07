package CoNVaDING;

require 5.22.0;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use CoNVaDING ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

#Look in the perl script
#our $VERSION = '0.01';


# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=encoding utf8


=head1 NAME

CoNVaDING - CoNVaDING (Copy Number Variation Detection In Next-generation sequencing Gene panels)


=head1 SYNOPSIS

Nothing in here is meant for public consumption. Use CoNVaDING.pl from the command line.

    ./CoNVaDING.pl --help


=head1 DESCRIPTION

CoNVaDING (Copy Number Variation Detection In Next-generation sequencing Gene panels) was designed 
for small (single-exon) copy number variation (CNV) detection in high coverage next-generation sequencing 
(NGS) data, such as obtained by analysis of smaller targeted gene panels.CoNVaDING makes use of a group 
of (at least 30) possible control samples from which the samples with the most similar overall pattern are
selected as control samples. These samples are then used for read-depth normalization on all (autosomal)
targets and on all targets per gene. CNV prediction is based on a combination of ratio scores and Z-scores
of the sample of interest compared to the selected controlsamples. Quality (QC) metrics are calculated per
sample and per analyzed target. Output is generated on three levels:

=over 

=item *

longlist	This list contains all calls, disregarding the target quality.

=item *

shortlist	This list contains a subset of the longlist, filtered on within sample target QC metrics.

=item *

final list	This list contains a subset of the shortlist, filtered on target QC metrics obtained from other samples.

=back

CoNVaDING has been written for use of CNV detection in high coverage NGS data (at 
least ~200x). This will fail with coverages >8000 due to samtools, then it is advised to change the setting
 -samtoolsdepthmaxcov 8000 to a higher value. With lower coverages it might still work, but more targets
 will fail QC metrics.

The program is written in perl and has dependencies on specific perl libraries
as well as on samtools version 1.3 or higher and bedtools 2.25.0. This is currently 
also tested on the versions served by apt on travis.


=head2 EXPORT

None by default.


=head1 SEE ALSO

https://github.com/molgenis/CoNVaDING


=head1 CITATION

If you use this tool in publications cite:

Johansson LF, van Dijk F, de Boer EN, van Dijk-Bos KK, Jongbloed JD, 
van der Hout AH, Westers H, Sinke RJ, Swertz MA, Sijmons RH, 
Sikkema-Raddatz B. CoNVaDING: Single Exon Variation Detection in Targeted NGS Data. 
Hum Mutat. 2016 May;37(5):457-64. doi: 10.1002/humu.22969. Epub 2016 Feb 24.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

original version at https://github.com/molgenis/CoNVaDING

current version at https://github.com/duartemolha/CoNVaDING_reload

=cut
