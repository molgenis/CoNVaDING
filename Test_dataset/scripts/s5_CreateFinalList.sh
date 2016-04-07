CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
DATADIR="/PATH/TO/Test_dataset/"

perl $CONVADINGDIR/CoNVaDING.pl \
-mode CreateFinalList \
-inputDir $DATADIR/results/StartWithBestScore \
-outputDir $DATADIR/results/CreateFinalList \
-targetQcList $DATADIR/results/GenerateTargetQcList/targetQcList.txt
