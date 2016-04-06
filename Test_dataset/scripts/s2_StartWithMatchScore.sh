CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl $CONVADINGDIR/CoNVaDING.pl \
-mode StartWithMatchScore \
-inputDir $DATADIR/results/StartWithAvgCount \
-controlsDir $DATADIR/controls \
-outputDir $DATADIR/results/StartWithMatchScore
