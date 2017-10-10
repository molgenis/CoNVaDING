CONVADINGDIR="../.."
DATADIR="../../Test_dataset"

perl $CONVADINGDIR/CoNVaDING_old.pl \
-mode CreateFinalList \
-inputDir $DATADIR/results \
-outputDir $DATADIR/results \
-targetQcList $DATADIR/results/targetQcList.txt
