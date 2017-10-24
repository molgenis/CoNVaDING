CONVADINGDIR="../.."
DATADIR="../../Test_dataset"

perl $CONVADINGDIR/CoNVaDING_old.pl \
-mode StartWithMatchScore \
-inputDir $DATADIR/results \
-controlsDir $DATADIR/controls \
-outputDir $DATADIR/results
