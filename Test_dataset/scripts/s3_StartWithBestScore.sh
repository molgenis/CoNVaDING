CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl $CONVADINGDIR/CoNVaDING.pl \
-mode StartWithBestScore \
-outputDir $DATADIR/results/StartWithBestScore \
-controlsDir $DATADIR/controls \
-inputDir $DATADIR/results/StartWithMatchScore
