CONVADINGDIR="../.."
DATADIR="../../Test_dataset"

perl $CONVADINGDIR/CoNVaDING_old.pl \
-mode GenerateTargetQcList \
-outputDir $DATADIR/results \
-controlsDir $DATADIR/controls \
-inputDir $DATADIR/controls
