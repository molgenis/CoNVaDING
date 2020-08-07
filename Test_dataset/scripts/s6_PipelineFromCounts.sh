CONVADINGDIR="../.."
DATADIR="../../Test_dataset"


perl $CONVADINGDIR/CoNVaDING.pl \
-mode PipelineFromCounts \
-inputDir $DATADIR/sample \
-bed $DATADIR/bedfile/Test_dataset_bedfile.bed \
-outputDir $DATADIR/results \
-controlsDir $DATADIR/controls
