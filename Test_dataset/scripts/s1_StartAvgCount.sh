CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl $CONVADINGDIR/CoNVaDING.pl \
-mode StartWithAvgCount \
-inputDir $DATADIR/sample \
-bed $DATADIR/bedfile/Test_dataset_bedfile.bed \
-outputDir $DATADIR/results/StartWithAvgCount \
-controlsDir $DATADIR/controls

