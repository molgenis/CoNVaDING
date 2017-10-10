CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl -d $CONVADINGDIR/CoNVaDING.pl \
-mode PipelineFromCounts \
-inputDir $DATADIR/sample \
-bed $DATADIR/bedfile/Test_dataset_bedfile.bed \
-outputDir $DATADIR/results/pipelineFrom_counts \
-controlsDir $DATADIR/controls
