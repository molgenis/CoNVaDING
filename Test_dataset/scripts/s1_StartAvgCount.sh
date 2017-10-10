CONVADINGDIR="../.."
DATADIR="../../Test_dataset"

perl $CONVADINGDIR/CoNVaDING_old.pl \
-mode StartWithAvgCount \
-inputDir $DATADIR/sample \
-bed $DATADIR/bedfile/Test_dataset_bedfile.bed \
-outputDir $DATADIR/results \
-controlsDir $DATADIR/controls

