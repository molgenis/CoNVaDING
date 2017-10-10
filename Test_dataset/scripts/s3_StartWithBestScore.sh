CONVADINGDIR="../.."
DATADIR="../../Test_dataset"

perl $CONVADINGDIR/CoNVaDING_old.pl \
-mode StartWithBestScore \
-outputDir $DATADIR/results \
-controlsDir $DATADIR/controls \
-inputDir $DATADIR/results
