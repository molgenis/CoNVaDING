# CoNVaDING User Guide

Written for release CoNVaDING v1.2.1


# Introduction

CoNVaDING (Copy Number Variation Detection In Next-generation sequencing Gene panels) was designed for small (single-exon) copy number variation (CNV) detection in high coverage next-generation sequencing (NGS) data, such as obtained by analysis of smaller targeted gene panels. 

CoNVaDING makes use of a group of (at least 30) possible control samples from which the samples with the most similar overall pattern are selected as control samples. These samples are then used for read-depth normalization on all (autosomal) targets and on all targets per gene. CNV prediction is based on a combination of ratio scores and Z-scores of the sample of interest compared to the selected controlsamples.

Quality (QC) metrics are calculated per sample and per analyzed target. Output is generated on three levels:

- longlist: This list contains all calls, disregarding the target quality.

- shortlist: This list contains a subset of the longlist, filtered on within sample target QC metrics.
- final list: This list contains a subset of the shortlist, filtered on target QC metrics obtained from other samples.

CoNVaDING has been written for use of CNV detection in high coverage NGS data (at least ~200x). With lower coverages it might still work, but more targets will fail QC metrics.

The program is written in perl and has dependencies on specific perl libraries as well as on samtools version 1.3 or higher.



# Installation

The latest version of CoNVaDING can be downloaded from the [github](https://github.com/molgenis/CoNVaDING/releases) page.

CoNVaDING has several dependencies:

- [Perl](https://www.perl.org/)

- The [Statistics::Normality package](http://search.cpan.org/~mwendl/Statistics-Normality-0.01/lib/Statistics/Normality.pm) for perl  

- [Samtools](http://www.htslib.org/) (v1.3 or higher)

After installation the samtools executable has to be added to your local environment using the $PATH variable.

This version is known to be compatible with SAMtools version 1.3



# General comments

Targets, or the region of interest, usually consist of an exon with some flanking bases, but can be anything that you specify, for instance part of an exon.

Because of the fact that the average coverage of the specified targets is used for the calculations, the resolution of the calls is also dependent on those targets. This means that the exact break point can not be determined. It might be within the last target of the call, or in the intron flanking this exon. Sub-exonic CNV's can also be missed if they if they do not have enough effect on the average coverage of the target.



# Analysis

The analysis consists of three steps that have to be run separately and is based on a list of targets in a bed format.

CoNVaDING can be started using the following command:
```bash
  perl ./CoNVaDING.pl
```

If no options are used the help menu will appear.
```bash
  Usage: ./CoNVaDING.pl <mode> <parameters>
  -h	  This manual.
  -mode		Mode to run in, one of the following required:
			PipelineFromBams
			  	Run entire pipeline starting with a folder of BAM files as input, to enable duplicate
  				removal use the rmdup variable.
				REQUIRED:
  				[-inputDir, -outputDir, -bed, -controlsDir]
  				OPTIONAL:
  				[-rmDup, -useSampleAsControl --controlSamples]
				[-regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]
				[-percentageLessReliableTargets]
			PipelineFromCounts
  				Run entire pipeline starting with a folder of of normalised count files as input
				REQUIRED:
  				[-inputDir, -outputDir, -bed, -controlsDir]
  				OPTIONAL:
  				[-rmDup, -useSampleAsControl --controlSamples]
				[-regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]
				[-percentageLessReliableTargets]
			addToControls
  				REQUIRED:
				[-inputDir, -outputDir, -bed]
				OverWritten parameteres:
				[-useSampleAsControl} this is necessarily true on this mode
				[-controlsDir] This is set to the same as -outputDir
  				OPTIONAL:
				[-rmDup, ] 
					make sure that all the control files already present in the control folder where all calculated with
					or without -rmDup. A mix of the 2 types in the control folders might bring unexpected results. 
					In the future, I might make the control files contain an indication of how their where originally 
					calculated and get the software to automatically pick the control that where processed the same way 
					as the input samples
  			StartWithBam :
  				Start with BAM files as input, to enable duplicate
  				removal use the rmdup variable.
  				REQUIRED:
  				[-inputDir, -outputDir, -bed, -controlsDir]
  				OPTIONAL:
  				[-rmDup, -useSampleAsControl]
  			StartWithAvgCount :
  				Start with Average Count files as input. This is a five column text file
  				with predefined column names. Please read the manual for instructions.
  				REQUIRED:
  				[-inputDir, -outputDir, -bed, -controlsDir]
  				OPTIONAL:
  				[-useSampleAsControl]
  			StartWithMatchScore :
  				Start with Normalized Coverage files as input.
  				REQUIRED:
  				[-inputDir, -outputDir, -controlsDir]
  				OPTIONAL:
  				[-controlSamples, -sexChr]
  			StartWithBestScore :
  				Best score analysis using Match score files as input.
  				REQUIRED:
  				[-inputDir, -outputDir, -controlsDir]
  				OPTIONAL:
  				[-regionThreshold, -sexChr, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh]
  			GenerateTargetQcList :
  				Generate a target QC list to use as input for finallist creation.
  				REQUIRED:
  				[-outputDir, -controlsDir]
  				OPTIONAL:
  				[-controlSamples, -regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]
  			CreateFinalList :
  				Creates the final list using the target QC list for filtering.
  				REQUIRED:
  				[-inputDir, -targetQcList, -outputDir]
  				OPTIONAL:
  				[-percentageLessReliableTargets]
  PARAMETERS:
  -inputDir             Input directory, depending on the analysis mode this contains
                        BAM, AvgCount, normalized coverage or match score files.
  -bed                  Input file specifying regions to analyze in BED format.
  -outputDir            Output directory to write results to.
  -controlsDir          Directory containing control samples.
  -targetQcList         Path to file containing target QC values.
  -controlSamples       Number of samples to use in Match score analysis. DEFAULT: 30
  -regionThreshold      Percentage of all control samples differing more than 3
                        standard deviations from mean coverage of a region in the specified
                        BED file to exlude from sample CV calculation. DEFAULT: 20
  -rmDup                Switch to enable duplicate removal when using BAM files as input.
  -sexChr               Switch to include sex chromosomes in analysis.
  -useSampleAsControl   Switch to use samples as control. Example: when using BAM
                        files to create count files and subsequentially use the
                        generated count files as controls.
  -ratioCutOffLow       Lower ratio cutoff value. Region ratio values below this
                        threshold are marked as deletion. DEFAULT: 0.65
  -ratioCutOffHigh      Higher ratio cutoff value. Region ratio values above this
                        threshold are marked as duplication. DEFAULT: 1.4
  -zScoreCutOffLow      Lower Z-score cutoff value. Regions with a Z-score below
                        this threshold are marked as deletion. DEFAULT: -3
  -zScoreCutOffHigh     Higher Z-score cutoff value. Regions with a Z-score above
                        this threshold are marked as duplication. DEFAULT: 3
  -sampleRatioScore     Sample CV z-score cutoff value. Sample with a ratio
                        score below this value are excluded from analysis. DEFAULT: 0.09. NOTE: this variable now is named "SAMPLE CV" in all outputs!
  -percentageLessReliableTargets	Target labelled as less reliable in percentage
                        of control samples. DEFAULT: 20
```

## Create normalized count files

The first step in the analysis is to create normalized count files. This can be done in two ways, from a bam file or from a text file including mean coverage per target.


### StartWithBam

If a bam file is used CoNVaDING will use samtools to calculate the mean coverage for each target. For this type of analysis the StartWithBam mode has to be selected.

The basic analysis starts as follows: 
```bash  
  perl ./CoNVaDING.pl \
  -mode StartWithBam \
  -inputDir /PATH/TO/INPUTDIR \
  -controlsDir /PATH/TO/CONTROLSDIR \
  -outputDir /PATH/TO/OUTPUTDIR \
  -bed /PATH/TO/DIR/target_bedfile.bed
```

All bamfiles should be stored in the same folder, which can be specified with the 'inputDir' option
The 'outputDir' option should specify the path to the folder in which the normalized coverage files *.aligned.only.normalized.coverage.txt* should be stored.

The bed file should contain the regions of interest seperated in four columns specifying the chromosome, start position, stop position and the gene. No headers should be included. Make sure the style of the bed file matches the style of the bam file. UCSC style files show 'chr1, chr2, etc.' in the header, while in other styles this is '1, 2, etc.'. The naming of the chromosomes in the bed file should be the same as in the bam file header.


It is important that the gene column has the exact same gene name for every target of the same gene, 
because these names are used to cluster targets for a normalization based on the targets belonging to the same gene. 
The bedfile should be sorted on chromosome and start position.

You Can use the 5th column to further describe the targets and these will be displayed in the output report files

Bed file example: 
```bash
  2       96919506        96919893        TMEM127	Exon1
  2       96920531        96920775        TMEM127	Exon2
  2       96930836        96931159        TMEM127	Exon3
  2       215593360       215593772       BARD1	Exon1
  2       215595095       215595272       BARD1	Exon2
  2       215609751       215609923       BARD1	Exon3
  2       215610406       215610618       BARD1	Exon4
  2       215617131       215617319       BARD1	Exon5
  2       215632166       215632418       BARD1	Exon6
  2       215633916       215634076       BARD1	Exon7
  2       215645244       215646273       BARD1	Exon8
  2       215656981       215657209       BARD1	Exon9
  2       215661745       215661881       BARD1	Exon10
  2       215674096       215674333       BARD1	Exon11
  3       10183492        10183911        VHL	Exon1
  3       10188158        10188360        VHL	Exon2
  3       10191431        10191689        VHL	Exon3
```

The analysis options can be further extended:

If a control set is not yet present, or if the samples that are analyzed have to be added to the control set the following options should be added:
```bash
  -useSampleAsControl
  -controlsDir /PATH/TO/CONTROLSDIR
```
The 'useSampleAsControl' option specifies that the samples have to be used as a control sample later on. The 'controlsDir' is the location where the normalized coverage files of the control samples will be stored.


If duplicates have to be removed before coverage calculations use the following option:
```bash
  -rmdup
```
This is advisable for capturing data, but should not be done for amplicon data.




### StartWithAvgCount

If no bam files are present the analysis can also start with a text file specifying average counts per target. This file shoud contain the headers as shown in the example below.
This enables the use of alternative analysis software.


An example of a text file that can be used in this mode:
```bash
  CHR     START           STOP            GENE    REGION_COV
  2       96919506        96919893        TMEM127 209.606 
  2       96920531        96920775        TMEM127 230.959 
  2       96930836        96931159        TMEM127 127.735 
  2       215593360       215593772       BARD1   273.726 
  2       215595095       215595272       BARD1   297.522 
  2       215609751       215609923       BARD1   230.191
  2       215610406       215610618       BARD1   224.822
  2       215617131       215617319       BARD1   204.979 
  2       215632166       215632418       BARD1   211.352 
  2       215633916       215634076       BARD1   240.627 
  2       215645244       215646273       BARD1   281.97  
  2       215656981       215657209       BARD1   137.293 
  2       215661745       215661881       BARD1   264.81  
  2       215674096       215674333       BARD1   127.689 
  3       10183492        10183911        VHL     174.233 
  3       10188158        10188360        VHL     230.704 
  3       10191431        10191689        VHL     226.012 
```

When this mode is used the analysis is started as follows:
```bash
  perl ./CoNVaDING.pl \
  -mode StartWithAvgCount \
  -inputDir /PATH/TO/INPUTDIR \
  -outputDir /PATH/TO/OUTPUTDIR \
  -bed /PATH/TO/DIR/target_bedfile.bed
```

Also here the following options can be used if the samples should be used as a control set in later steps:
```bash
  -useSampleAsControl
  -controlsDir /PATH/TO/CONTROLSDIR
```
The -rmdup option is not available in this mode. If necessary the duplicates should have been removed before calculating the mean coverage per target.


## Selecting the most informative control samples


### StartWithMatchScore

The next step in the analysis is selecting the control samples estimated to be the most informative.


If a bam file is used CoNVaDING will use samtools to calculate the mean coverage for each target. For this type of analysis the StartWithBam mode has to be selected.

The basic analysis starts as follows: 
```bash  
  perl ./CoNVaDING.pl \
  -mode StartWithMatchScore \
  -inputDir /PATH/TO/INPUTDIR \
  -outputDir /PATH/TO/OUTPUTDIR \
  -controlsDir /PATH/TO/CONTROLSDIR
```

The 'inputDir' option should specify the path to the folder in which the normalized coverage files *.aligned.only.normalized.coverage.txt* are be stored (the outputfolder of the previous step).

The 'outputDir' option should specify the path to output folder. The script will produce two types of output files:

*.best.match.score.txt* shows the matchscore and the paths to the selected control samples.

*.normalized.autosomal.coverage.all.controls.txt* show the normalized coverage for all possible control samples


The 'controlsDir' option should show the directory in which the control samples that have to be used are stored.


The analysis options can be further extended:

On default only targets located on autosomal chromosomes will be analyzed. If some targets are located on the sex chromosomes the following option should be added:
```bash
  -sexChr
```
Note that for this option only samples of the same sex as the sample of interest can be used as possible control samples.

On default 30 samples are selected to create the control group. If you wish to use a different number of control samples this can be indicated with the option:
```bash
  -controlSamples 40
```
to select for instance the 40 best matching samples.


## CNV Detection


### StartWithBestScore

The last step in the analysis is the CNV detection itself.

The basic analysis starts as follows: 
```bash  
  perl ./CoNVaDING.pl \
  -mode StartWithBestScore \
  -inputDir /PATH/TO/INPUTDIR \
  -outputDir /PATH/TO/OUTPUTDIR \
  -controlsDir /PATH/TO/CONTROLSDIR
```

The 'inputDir' option should specify the path to the folder in which the *.best.match.score.txt* files are stored (the outputfolder of the previous step).

The 'outputDir' option should specify the path to output folder. The script will produce four types of output files:

*.best.score.log* show the used control samples, the sample ratio score and the omitted regions for the sample ratio score calculation.

*.best.score.longlist.txt* contains all calls, regardless of the target quality

*.best.score.shortlist.txt* contains the high quality calls based on within sample target QC

*.best.score.totallist.txt* contains information about all targets (ratio scores, Z-scores, QC) 
 

The analysis options can be further extended:


The sample ratio calculation is based on calculation the coefficient of variation of the normalized targets of the sample of interest. In this calculation highly variable targets are excluded. On default a target is considered highly variable if after transforming the normalized target ratio's of all samples in the possible control group 20 percent or more of the samples have a Z-score outside the -3 to 3 range. This percentage can be altered using the 'regionTreshold' option.
For a threshold of 30 percent of the samples for instance the following option can be used:
```bash
  -regionThreshold 30
```
To alter the ratio thresholds when making a call for a deletion of duplication for a region during the analysis, the ratioCutOffLow and ratioCutOffHigh parameters can be used.
To apply a threshold of ratio score below 0.65 for a deletion and above 1.4 for duplication use:
```bash
  -ratioCutOffLow 0.65
  -ratioCutOffHigh 1.4
```
The same thresholds for calling a deletion or duplication can also be applied using the Z-score value cutoff.
To call a deletion when the Z-score is below -3 or duplication when the Z-score is above 3 use:
```bash
  -zScoreCutOffLow=-3
  -zScoreCutOffHigh 3
```

Note that to use the minus variable in the -zScoreCutOffLow an '=' sign is needed instead of a space.

To finetune the variant list one can generate a list of targets which in general are of lower quality in all possible controlsamples and apply this as a filter to generate a final list of high quality calls. This can be done by executing two steps:


### GenerateTargetQcList

To generate the list of targets and corresponding quality thresholds run:
```bash
  perl ./CoNVaDING.pl \
  -mode GenerateTargetQcList \
  -inputDir /PATH/TO/CONTROLSDIR \
  -outputDir /PATH/TO/OUTPUTDIR \
  -controlsDir /PATH/TO/CONTROLSDIR
```
For this analysis, the same region threshold, ratio cutoffs and Z-score cutoffs as explained above can be altered using their corresponding parameters.


### CreateFinalList

To apply the generated list of sample target QCs to the *.best.score.shortlist.txt* files execute:
```bash
  perl ./CoNVaDING.pl \
  -mode CreateFinalList \
  -inputDir /PATH/TO/BESTSCOREOUTPUT \
  -targetQcList /PATH/TO/TARGETQCLISTFILE \
  -outputDir /PATH/TO/OUTPUTDIR
```

To change the percentage of samples in which a target can be labelled as less reliabe, for example in 20 percent of the samples, use the option:
```bash
  -percentageLessReliableTargets 20
```
This produces the following output file:

*.finallist.txt* contains all final calls, basically a filtered shortlist file.


# Test dataset

A test dataset can be downloaded from [here](https://github.com/molgenis/CoNVaDING/tree/master/Test_dataset).
The dataset contains one sample and fifty control samples. The scripts described above are also included in the test dataset folder.

### StartWithAvgCount

Analysis starts from the StartWithAvgCount step:

```bash
  CONVADINGDIR="/PATH/TO/CoNVaDINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
  
  perl $CONVADINGDIR/CoNVaDING.pl \
  -mode StartWithAvgCount \
  -inputDir $DATADIR/sample \
  -bed $DATADIR/bedfile/Test_dataset_bedfile.bed \
  -outputDir $DATADIR/results/StartWithAvgCount \
  -controlsDir $DATADIR/controls
```

Running this step should create a resultsfolder within the Test_dataset folder containing a normalized coverage file. 

### StartWithMatchScore

Subsequently the best matching samples can be determined using the StartWithMatchScore option:

```bash
  CONVADINGDIR="/PATH/TO/CoNVaDINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
  
  perl $CONVADINGDIR/CoNVaDING.pl \
  -mode StartWithMatchScore \
  -inputDir $DATADIR/results/StartWithAvgCount \
  -controlsDir $DATADIR/controls \
  -outputDir $DATADIR/results/StartWithMatchScore
  CONVADINGDIR="/PATH/TO/CoNVaINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
```

The best thirty samples should be selected, as shown below:

```bash
  #######################################
  Selecting best 30 control samples for analysis..
  #######################################
  Control: Control39.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0528891755854009
  Control: Control44.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0534887805765252
  Control: Control43.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0542475244113742
  Control: Control36.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0554800052287148
  Control: Control35.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0562685136912045
  Control: Control42.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0574222886538507
  Control: Control40.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0598553222944651
  Control: Control48.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0600383118556296
  Control: Control50.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0608676166830126
  Control: Control45.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0625870879199154
  Control: Control10.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0666569275813975
  Control: Control34.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0675543014488717
  Control: Control28.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0679211728595409
  Control: Control38.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0691059245241418
  Control: Control07.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0745447949437231
  Control: Control23.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0781103186115519
  Control: Control02.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0792689738060643
  Control: Control47.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0811643176924381
  Control: Control15.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.081433371903405
  Control: Control12.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0821882398310768
  Control: Control22.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.084360960603957
  Control: Control24.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0866925360424663
  Control: Control09.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0894002915947557
  Control: Control41.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0896640102541329
  Control: Control37.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0899251768650165
  Control: Control16.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0915173887217015
  Control: Control46.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0939866051324493
  Control: Control27.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.0966658138493035
  Control: Control32.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.100018174789515
  Control: Control17.aligned.only.normalized.coverage.txt                 Avg abs diff: 0.101351812101807
  #######################################
```

### StartWithBestScore

Now the CNV detection can be performed, using the StartWithBestScore option:

```bash
  CONVADINGDIR="/PATH/TO/CoNVaDINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
  
  perl $CONVADINGDIR/CoNVaDING.pl \
  -mode StartWithBestScore \
  -outputDir $DATADIR/results/StartWithBestScore \
  -controlsDir $DATADIR/controls \
  -inputDir $DATADIR/results/StartWithMatchScore
```

A StartWithBestScore folder is created in the results folder, containing four files.

- Sample1.average.counts.best.score.log
- Sample1.average.counts.best.score.longlist.txt
- Sample1.average.counts.best.score.shortlist.txt
- Sample1.average.counts.best.score.totallist.txt

The log file will show the Sample QC: SAMPLE_CV: 0.0749944561392519

and the Match QC: MEAN_AVERAGE_BEST_MATCHSCORE: 0.0748225246685803

The longlist should contain six calls:


*Sample1.average.counts.best.score.longlist.txt*
```bash
  CHR START      STOP       GENE    GENETARGETS	NUMBER_OF_TARGETS   NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST  ABBERATION
  1   156104958  156105124  LMNA    exon1-exon2	       1                   1                                         DEL
  2   179511192  179511307  TTN     exon1-exon10	   1                   1                                         DEL
  6   7576507    7578140    DSP     exon1	           3                   3                                         DUP
  12  21997397   21997507   ABCC9   exon5              1                   1                                         DUP
  18  28647961   28681955   DSC2    regionA           17                  17                                         DEL
  18  29078195   29102233   DSG2    intron1	           6                   5                                         DEL
```

The totallist contains the information of all targets and shows the ratio's and Z-scores and the coefficient of variation of each target of the control set ratio's. If the coefficient of variation of the target (AUTO_VC) is too high (above 0.10) the target QC fails and the target is labelled low quality. 


*Sample1.average.counts.best.score.totallist.txt*
```bash
  CHR	START	   STOP	       GENE	  .. .. AUTO_VC     .. .. .. ABBERATION   QUALITY
  1     156104958   156105124  LMNA	  .. .. 0.06058213  .. .. .. DEL          .
  ..						
  2     179511192   179511307  TTN	  .. .. 0.095184824 .. .. .. DEL          .
  ..
  6     7576507	   7576710	  DSP	  .. .. 0.051947043 .. .. .. DUP          .
  6	   7577172	  7577296   DSP	  .. .. 0.094469481 .. .. .. DUP          . 
  6	   7577992	  7578140   DSP	  .. .. 0.061275822 .. .. .. DUP          .
  ..						
  12    21997397   21997507   ABCC9   .. .. 0.120684521 .. .. .. DUP         LOW_QUALITY
  ..						
  18    28647961   28648199   DSC2	  .. .. 0.063042636 .. .. .. DEL          .
  18    28648255   28648331   DSC2	  .. .. 0.160488643 .. .. .. DEL          LOW_QUALITY
  18    28648840   28649138   DSC2	  .. .. 0.048139583 .. .. .. DEL          .
  18    28650672   28650837   DSC2	  .. .. 0.080681945 .. .. .. DEL          .
  18    28651551   28651828   DSC2	  .. .. 0.061021188 .. .. .. DEL          .
  18    28654629   28654894   DSC2	  .. .. 0.052091653 .. .. .. DEL          .
  18    28659793   28659976   DSC2	  .. .. 0.118304686 .. .. .. DEL          LOW_QUALITY
  18    28660042   28660339   DSC2	  .. .. 0.062031772 .. .. .. DEL          .
  18    28662184   28662410   DSC2	  .. .. 0.16236608  .. .. .. DEL          LOW_QUALITY
  18    28662872   28663047   DSC2	  .. .. 0.097764448 .. .. .. DEL          .
  18    28666519   28666727   DSC2	  .. .. 0.07031019  .. .. .. DEL          .
  18    28667612   28667797   DSC2	  .. .. 0.081027693 .. .. .. DEL          .
  18    28669382   28669578   DSC2	  .. .. 0.074671788 .. .. .. DEL          .
  18    28670971   28671133   DSC2	  .. .. 0.064327342 .. .. .. DEL          .
  18    28672044   28672284   DSC2	  .. .. 0.067765317 .. .. .. DEL          .
  18    28673502   28673627   DSC2	  .. .. 0.136768702 .. .. .. DEL          LOW_QUALITY
  18    28681846   28681955   DSC2	  .. .. 0.188210597 .. .. .. DEL          LOW_QUALITY
  18    29078195   29078280   DSG2	  .. .. 0.10492824  .. .. .. DEL          LOW_QUALITY
  18    29098182   29098258   DSG2	  .. .. 0.069570847 .. .. .. DEL          .
  18    29099746   29099921   DSG2	  .. .. 0.058676083 .. .. .. DEL          .
  18    29100746   29100948   DSG2	  .. .. 0.148595096 .. .. .. DEL          LOW_QUALITY
  18	29101042   29101227   DSG2	  .. .. 0.102009337 .. .. .. DEL          LOW_QUALITY
  18	29102026   29102233   DSG2	  .. .. 0.091922155 .. .. .. DEL          .
```

Since the ABCC9 exon in which a duplication was detected has a low quality, this exon is filtered from the shortlist


*Sample1.average.counts.best.score.shortlist.txt*
```bash
  CHR START      STOP       GENE   GENETARGETS	 NUMBER_OF_TARGETS   NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST  ABBERATION
  1   156104958  156105124  LMNA   exon1-exon2	        1                   1                                         DEL
  2   179511192  179511307  TTN    exon1-exon10         1                   1                                         DEL
  6   7576507    7578140    DSP    exon1	            3                   3                                         DUP
  18  28647961   28681955   DSC2   exon5               17                  17                                         DEL
  18  29078195   29102233   DSG2   regionA              6                   5                                         DEL
```                                	    

### GenerateTargetQcList

For final filtering the TargetQCList is made using all control samples:

```bash
  CONVADINGDIR="/PATH/TO/CoNVaDINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
  
  perl $CONVADINGDIR/CoNVaDING.pl \
  -mode GenerateTargetQcList \
  -outputDir $DATADIR/results/GenerateTargetQcList \
  -controlsDir $DATADIR/controls \
  -inputDir $DATADIR/controls
```

### CreateFinalList

Finally the shortlist is filtered usint the targetQClist:

```bash
  CONVADINGDIR="/PATH/TO/CoNVaDINGDIR/"
  DATADIR="/PATH/TO/Test_dataset/"
  
  perl $CONVADINGDIR/CoNVaDING.pl \
  -mode CreateFinalList \
  -inputDir $DATADIR/results/StartWithBestScore \
  -outputDir $DATADIR/results/CreateFinalList \
  -targetQcList $DATADIR/results/GenerateTargetQcList/targetQcList.txt
```

The call of the titin exon had sufficient quality within the analysed sample. However, the exon performed poorly in a large portion of the control samples. Therefore, the call is filtered from the final list, leaving four calls. Notice that the whole gene deletion of DSC2 and the six exon deletion of DSG2 consist of consecutive targets. It is possible that there is one big deletion, containing both genes. However, CoNVaDING will always treat CNVs in different genes as seperate calls.  


*Sample1.average.counts.best.score.shortlist.finallist.txt*
```bash
  CHR START      STOP       GENE  GENETARGETS	   NUMBER_OF_TARGETS   NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST  ABBERATION
  1   156104958  156105124  LMNA  exon1-exon2	          1                   1                                         DEL
  6   7576507    7578140    DSP   exon1-exon10            3                   3                                         DUP
  18  28647961   28681955   DSC2  exon1	                 17                  17                                         DEL
  18  29078195   29102233   DSG2  exon5                   6                   5                                         DEL
```                                      


 
# QC Thresholds

If the following thresholds are exceeded using default settings the CNV calling is less reliable.
The target ratio is also used to filter calls for the shortlist. Both QC metrics are used for filtering the final list.

Sample CV: 0.09

Target CV: 0.10


#FAQ

Which are the criteria of adding a call to the longlist, shortlist or finallist?
	
- All calls passing threshold criteria in the StartWithBestScore function are added to the longlist, regardless of there quality. These will be all targets showing DEL, DUP, HOM_DEL or HOM_DUP in the totallist. This list provides a high sensitivity, but has a relatively low specificity. 
- The shortlist takes into account the variability of the target in the selected control samples. If the target CV in the used control samples is higher than 0.10 the target is considered LOW_QUALITY, as shown in the totallist. Calls consisting only of low quality targets are filtered from the longlist and are not present on the shortlist.
- The finallist is a further filtering of the shortlist. If one or more targets of a call are not considered low quality in the sample, the call will be on the shortlist. However, if all targets of the call are considered low quality in a higher percentage of samples having a SAMPLE_CV lower than 0.09, than specified in the -regionThreshold option (default 20%), the call will not be present on the finallist. This list has the highest specificity, but CNV's in low quality regions can be missed.

If a target has a coverage of 0 CoNVaDING crashes. For this target the normalized coverage file shows -NaN as region coverage.

- This is the result of a bug in CoNVaDING v1.1.6. This should be fixed in version 1.2.0.

A call I expect is present on the shortlist, but not on the finallist.

- If the sample containing the CNV is present in the ControlsDir a call can be filtered because of itself. The reason for this is that the sample will be sometimes selected as a control sample. The deviating ratio of reads in the (here control) sample can cause a target to be considered as low quality. This might increase the number of control samples considering the target as low quality to above the threshold, causing the call to be filtered.
- This can be prevented by performing a preliminary analysis using all the samples in the ControlsDir as input as well as controls for the StartWithMatchScore and StartWithBestScore option. Samples showing calls in the shortlist should be removed from the ControlsDir. After this the real analysis can be started.

What is the difference between LOW_QUALITY and LOW_QUALITY, FAILED?

- Both are low quality targets. However, for the LOW_QUALITY, FAILED targets all targets belonging to the same gene are of low quality. In other words, analysis has failed for the entire gene. In LOW_QUALITY targets at least one of the other targets belonging to the same gene has a good quality. 

Can CoNVaDING be used on exome data?

- The performance of CoNVaDING is mainly dependent on the number of reads covering a target. If the exome coverage is high enough, CoNVaDING should work well. The optimum coverage is 200x or higher. Lower coverages will result in a higher number of low quality targets. Increasing the target size to include more reads in the target, for instance by adding 100bp flanking regions, may increase the quality.
- CoNVaDING can take some time to run when used on large bam files.
- The analysis speed of the StartWithBam function can be increased by running several instances of CoNVaDING in parallel. All bam files in the inputDir will be processed in sequence. Creating two inputDirs each containing half of the bam files can half the analysis time in case enough recources are available.

All the targets on the X chromosome are of low quality.

- Make sure all control samples are of the same gender. 

When setting the -zScoreCutOffLow parameter on X I get 'Unknown option :X'

- If the parameteres in the print to the screen are set correctly the output is correct. However, to prevent the error message instead of a space use an = sign to connect the variable. For instance:\

```
-zScoreCutOffLow=-3
```


# Changelog

##v1.3 CoNVaDING_reload
- Added Pipeline modes
- added support for target description via 5th column on bed file input
- added bed file outputs
- lots of close cleanup

##v1.2.1
- Enabled merging HOM_DEL/HOM_DUP and DEL/DUP calls. Now everything works correctly.  

##v1.2.0.1
- Disabled merging HOM_DEL/HOM_DUP and DEL/DUP calls as implemented in v1.2.0, because of a bug resulting in showing the start position of the last exon in the call in the longlist instead of the first exon in the call.

##v1.2.0

###Major
- In version 1.1.6 after running the StartWithBam function analyzing a region having a coverage of 0, the normalized coverage would have a value of -NaN. This issue is fixed in version 1.2.0, by adding the -a flag to the samtools depth function. However, this causes a dependency on samtools version 1.3 or higher. 
- In version 1.1.6 HOM_DEL and HOM_DUP calls were interpreted as separate events as DEL and DUP calls. This could cause calls to be cut into several pieces and create gaps in the calls because of downstream filtering. Now, HOM_DEL calls are interpreted as part of a DEL call as well as a separate call. For instance if exons 1 and 3 have a DEL call in the totallist and exon 2 a HOM_DEL call. A deletion of exons 1 to 3 will be called as well as a homozygous deletion of exon 2. Similarly HOM_DUP events are considered as part of DUP events as well as a separate call. Note that in case all targets of a call are considered HOM_DEL or HOM_DUP the call will appear twice on the longlist, shortlist and/or finallist.

###Minor
- In the StartWithMatchScore function the 'Avg abs diff incl sex' was showing the same value as the 'Avg abs diff autosomal'. Now the values are correct.
- A bedfile containing UCSC style notation 'chr1, chr2, etc.' is now be compatible with a header notation of '1, 2, etc.'
- The finallist is now called *finallist.txt instead of *shortlist.finallist.txt.
-  In the log file produced during the StartWithBestScore function SAMPLE_RATIO has been renamed SAMPLE_CV.
- The inputDir is not a requirement anymore for the GenerateTargetQCList function.
- In the GenerateTargetQC function the searchpattern is changed to prevent the step trying to process files other than normalized coverage files.
- The used commandline options and values are now printed to the stdout by default.


# Literature
[CoNVaDING: Single Exon Variation Detection in Targeted NGS Data](http://www.ncbi.nlm.nih.gov/pubmed/26864275)

# Tutorial
A video tutorial can be found [here](https://youtu.be/-geFWkvKZzE).

# Contact

Mailto: 

Duarte Molha <duartemolha@gmail.com> (CoNVaDING_reload fork)

Lennart Johansson <l.johansson@umcg.nl>

Freerk van Dijk <f.van.dijk02@umcg.nl>
