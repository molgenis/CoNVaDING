CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl $CONVADINGDIR/CoNVaDING.pl \
-mode GenerateTargetQcList \
-outputDir $DATADIR/results/GenerateTargetQcList \
-controlsDir $DATADIR/controls \
-inputDir $DATADIR/controls
