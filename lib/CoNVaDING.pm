package CoNVaDING;

use 5.022001;
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
longlist: This list contains all calls, disregarding the target quality.
=item *
shortlist: This list contains a subset of the longlist, filtered on within sample target QC metrics.
=item *
final list: This list contains a subset of the shortlist, filtered on target QC metrics obtained from other samples.
=back

CoNVaDING has been written for use of CNV detection in high coverage NGS data (at least ~200x). With lower coverages it might still work, but more targets will fail QC metrics.

The program is written in perl and has dependencies on specific perl libraries as well as on samtools version 1.3 or higher.

=head1 INSTALL

    perl Makefile.pl
    make
    make install


=head2 EXPORT

None by default.


=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

umcg-mterpstra, E<lt>umcg-mterpstra@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 by umcg-mterpstra

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.22.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
