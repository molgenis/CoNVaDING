#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use List::Util qw(sum);
use Math::Complex;
use POSIX qw(ceil);
use POSIX qw(floor);
use Statistics::Normality 'shapiro_wilk_test';
use File::Temp qw/ tempfile tempdir /;

######CHANGE VERSION PARAMETER IF VERSION IS UPDATED#####
#my $version_reload = "1.3";
#my $version = "1.2.1" ;
our $VERSION = '1.3.8';
my $version_reload = $VERSION;
my $version = $VERSION;

##############################################################################################
##############################################################################################
##   CoNVaDING reload V1, copy number variation detecting in next-generation                ##
##   sequencing gene panels                                                                 ##
##                                                                                          ##
##   this is a fork from the original program written by                                    ##
##   Copyright (C) 2015  Freerk van Dijk & Lennart Johansson                                ##
##                                                                                          ##
##                                                                                          ##
##   This file is part of CoNVaDING.                                                        ##
##                                                                                          ##
##   CoNVaDING is free software: you can redistribute it and/or modify                      ##
##   it under the terms of the GNU Lesser General Public License as published by            ##
##   the Free Software Foundation, either version 3 of the License, or                      ##
##   (at your option) any later version.                                                    ##
##                                                                                          ##
##   CoNVaDING is distributed in the hope that it will be useful,                           ##
##   but WITHOUT ANY WARRANTY; without even the implied warranty of                         ##
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                          ##
##   GNU Lesser General Public License for more details.                                    ##
##                                                                                          ##
##   You should have received a copy of the GNU Lesser General Public License               ##
##   along with CoNVaDING.  If not, see <http://www.gnu.org/licenses/>.                     ##
##############################################################################################
##############################################################################################

#Commandline variables
my $params = {};

#set defaults
$params->{regionThreshold}               = 20;
$params->{ratioCutOffLow}                = 0.65;
$params->{ratioCutOffHigh}               = 1.4;
$params->{zScoreCutOffLow}               = -3;
$params->{zScoreCutOffHigh}              = 3;
$params->{sampleRatioScore}              = 0.09;
$params->{percentageLessReliableTargets} = 20;
$params->{numBestMatchSamplesCmdL}       = 30;
$params->{mode}                          = "StartWithBam";
$params->{samtoolsdepthmaxcov}           = 8000;

#### get options
GetOptions(
    "mode:s"                          => \$params->{mode}, #For options, see list below (default="StartWithBam")
    "inputDir:s"                      => \$params->{inputdir},# required for all modes except GenerateTargetQcList
    "outputDir=s"                     => \$params->{outputdir}, # required always
    "controlsDir:s"                   => \$params->{controlsdir}, #optional
    "bed:s"                           => \$params->{bedfile}, #optional
    "controlSamples:s"                => \$params->{numBestMatchSamplesCmdL}, #optional
    "regionThreshold:s"               => \$params->{regionThreshold}, #optional
    "rmDup:s"                         => \$params->{rmdup},  #optional
    "sexChr:s"                        => \$params->{sexchr}, #optional
    "useSampleAsControl:s"            => \$params->{sampleAsControl}, #optional
    "ratioCutOffLow:s"                => \$params->{ratioCutOffLow}, #optional
    "ratioCutOffHigh:s"               => \$params->{ratioCutOffHigh}, #optional
    "zScoreCutOffLow:s"               => \$params->{zScoreCutOffLow}, #optional
    "zScoreCutOffHigh:s"              => \$params->{zScoreCutOffHigh}, #optional
    "sampleRatioScore:s"              => \$params->{sampleRatioScore}, #optional
    "targetQcList:s"                  => \$params->{targetQcList}, #optional
    "percentageLessReliableTargets:s" => \$params->{percentageLessReliableTargets}, #optional
    "samtoolsdepthmaxcov:i"           => \$params->{samtoolsdepthmaxcov},
    "h|help"                          => sub { usage() and exit(1)},
    "version"                         => sub { print "CoNVaDING relaod v".$version_reload." modified fork from CoNVaDING v".$version."\n" and exit(1)}
);
#Obligatory args
usage() and exit(1) unless $params->{mode};
usage() and exit(1) unless $params->{outputdir};
#Add more parameters later

#Check input parameters
unless ($params->{mode} eq "PipelineFromBams"      ||
        $params->{mode} eq "PipelineFromCounts"   ||
        $params->{mode} eq "addToControls"        ||
        $params->{mode} eq "StartWithBam"         ||
        $params->{mode} eq "StartWithAvgCount"    ||
        $params->{mode} eq "StartWithMatchScore"  ||
        $params->{mode} eq "StartWithBestScore"   ||
        $params->{mode} eq "GenerateTargetQcList" ||
        $params->{mode} eq "CreateFinalList") {
    die "Mode ".$params->{mode}." is not supported. Please read the manual.\n";
}

#If controlsdir is specified while running best score mode, throw warning
#Set default threshold values
if ($params->{mode} eq "StartWithBestScore"){
    if (defined $params->{controlsdir}){
        print "\n##### WARNING ##### WARNING #####\n";
        print "User specified parameter \"controlsDir\" not used in best score analysis.";
        print "\n##### WARNING ##### WARNING #####\n";
    }
}

#Check if output directory exists, otherwise create it
if (!-d $params->{outputdir}) { #Output directory does not exist, create it
    mkdir $params->{outputdir};
}

my @bedfile;
#Check if BED file is specified for StartWithBam and StartWithAvgCount, else die and throw error message
if ($params->{mode} eq "PipelineFromBams"    ||
    $params->{mode} eq "PipelineFromCounts" ||
    $params->{mode} eq "addToControls"      ||
    $params->{mode} eq "StartWithBam"       ||
    $params->{mode} eq "StartWithAvgCount"){
    if (not defined $params->{bedfile}){
        die "Required BED file not specified, please specify to continue analysis.\n";
    }
}elsif ($params->{mode} eq "StartWithMatchScore" ||
        $params->{mode} eq "StartWithBestScore"){ #Throw warning saying BED file is not used in this analysis.
    if (defined $params->{bedfile}){
        print STDERR "\n##### WARNING ##### WARNING #####\n";
        print STDERR "User specified parameter \"bed\" not used in this analysis.";
        print STDERR "\n##### WARNING ##### WARNING #####\n";
    }
}

if ($params->{mode} eq "addToControls"){ 
    if (not defined $params->{sampleAsControl}){
        print STDERR  "#####\n#Parameter sampleAsControl forcebly defined\n#####\n";
        $params->{sampleAsControl} = 1;
    }
    if (not defined $params->{controlsDir}){
        $params->{controlsDir} = $params->{outputDir};
    }
}

#Check if controlsdir exists for four param options
#if ($params->{mode} eq "StartWithBam" || $params->{mode} eq "StartWithAvgCount" || $params->{mode} eq "StartWithMatchScore" || $params->{mode} eq "GenerateTargetQcList"){ #Check if controlsdir is specified
if ($params->{mode} eq "PipelineFromBams"     ||
    $params->{mode} eq "PipelineFromCounts"  ||
    $params->{mode} eq "StartWithBam"        ||
    $params->{mode} eq "StartWithAvgCount"   ||
    $params->{mode} eq "StartWithMatchScore" ||
    $params->{mode} eq "addToControls"       ||
    $params->{mode} eq "GenerateTargetQcList"
    ){
    
    #Check if controlsdir is specified
    if (not defined $params->{controlsdir}) { #Throw error when controlsdir is not specified
        die "Directory for controlsamples (-controlsDir) is not specified, please specify to continue analysis.\n";
    }elsif (!-d $params->{controlsdir}) {
        #Check if controls directory exists, otherwise create it except in pipeline modes
        if ( $params->{mode} eq "StartWithAvgCount"   ||
             $params->{mode} eq "StartWithMatchScore" ||
             $params->{mode} eq "addToControls"){
            mkdir $params->{controlsdir};
        }elsif($params->{mode} eq "PipelineFromBams"     ||
               $params->{mode} eq "PipelineFromCounts"  ||
               $params->{mode} eq "GenerateTargetQcList"){
            die "Directory for controlsamples (-controlsDir) does not exist please specify to continue analysis.\n";
        }
    }
    if ($params->{mode} eq "PipelineFromBams" || $params->{mode} eq "PipelineFromCounts"){
        if (not defined $params->{targetQcList}){
           $params->{targetQcList} = $params->{outputdir}."/"."targetQcList.txt";
        }
    }
    #Check if controldir contains avg count files when useSampleAsControl paramater is not specified
    if (not defined $params->{sampleAsControl}){
        opendir(DIRHANDLE, $params->{controlsdir}) || die "Couldn't open directory ".$params->{controlsdir}.": $!\n";
        my @txt = grep {  /\.normalized.coverage.txt$/ && -f $params->{controlsdir}."/".$_ } readdir DIRHANDLE;
        closedir DIRHANDLE;
        if (@txt) { #not empty
            #continue
            my $numFiles = scalar(@txt); #If match score or best score mode is used, check if number of control sample files is greater than or equal to what user requested via parameter.
            if ($params->{mode} eq "StartWithMatchScore") {
                if ( $numFiles < $params->{numBestMatchSamplesCmdL}){
                    die "The number of controlsamples in ".$params->{controlsdir}." is less than specified number of ".$params->{numBestMatchSamplesCmdL}." controlsamples by \"controlSamples\" parameter.\n";
                }
            }
        }else{
            die "Directory specified by \"controlsDir\" parameter does not contain normalized coverage files.\n";
        }
    }    
}

#Check if input directory exists for all modes that require it
unless ($params->{mode} eq "GenerateTargetQcList"){
    if (not defined $params->{inputdir}) { #Throw error when inputdir is not specified
        die "Directory for input samples (-inputDir) is not specified, this is a requirement on all modes\n."
            ."with the exception GenerateTargetQcList. Please specify input folder to continue analysis.\n";
    }
}else{
    if (not defined $params->{inputdir}){
        $params->{inputdir} = $params->{controlsdir};
    }
}


#Checks when running in CreateFinalList mode
if ($params->{mode} eq "CreateFinalList"){ #Check if targetQcList is specified
    if (not defined $params->{targetQcList}) { #Throw error when targetQcList is not specified
        die "Target QC list (-targetQcList) is not specified, please specify to continue analysis.\n";
    }
}

#Retrieve and print starttime
my $starttime = localtime();

print STDERR "\nStarting analysis $starttime\n";

print STDERR "\n#######################################\n";
print STDERR "Parameteres in effect in this Run:\n";
foreach my $key (keys %{$params}){
    if(defined($params->{$key})){
        warn "\t$key:\t".$params->{$key}."\n";
    }else{
        warn "\t$key:\tundef=(FLAG in off position)\n";
    }
}
print STDERR "#######################################\n";


#sometimes the script want to modify  the outputdir and input dir to temp folders.
# these are therefore stored as additional parameteres
$params->{outputdirOriginal} = $params->{outputdir};
$params->{inputdirOriginal}  = $params->{inputdir};


#pipeline mode
if ($params->{mode} eq "PipelineFromBams" || $params->{mode} eq "PipelineFromCounts"){

    #mode to run complete pipeline from fresh bam inputs
    if ($params->{mode} eq "PipelineFromBams"){
        startWithBamMode();
     }   
    #mode to run complete pipeline from previous counts
    if ($params->{mode} eq "PipelineFromCounts"){
        startWithAvgCountMode()
    }
    $params->{inputdir} = $params->{outputdir} ;
    startWithMatchScoreMode();
    startWithBestScoreMode();
    
    unless (-e $params->{targetQcList}){
        $params->{inputdir} = $params->{controlsdir};
        generateTargetQcListMode();
        $params->{inputdir} = $params->{outputdir};
    }
    createFinalListMode();
#Start analysis from BAM file
}elsif ($params->{mode} eq "StartWithBam" || $params->{mode} eq "addToControls"){ 
    startWithBamMode();
#Start analysis from average count files
}elsif ($params->{mode} eq "StartWithAvgCount"){
    startWithAvgCountMode();
#Start analysis from match score
}elsif ($params->{mode} eq "StartWithMatchScore"){
    startWithMatchScoreMode();
#Start analysis from best score
}elsif ($params->{mode} eq "StartWithBestScore"){
    startWithBestScoreMode();
#Generate target QC list from all controlsamples
}elsif ($params->{mode} eq "GenerateTargetQcList"){
    generateTargetQcListMode();
#Create final list based on target QC file
}elsif ($params->{mode} eq "CreateFinalList"){ #Apply the target filtering using the file created in previous step
    createFinalListMode();
}

#Retrieve and print end time
my $endtime = localtime();
print STDERR "\nFinished analysis $endtime\n";



########################################################################################################
########## SUBS #################### SUBS #################### SUBS #################### SUBS ##########
########################################################################################################


##########################################################
## Main subs for each mode                              ##
##########################################################

sub startWithBamMode{
    if (defined $params->{rmdup}) { #Remove duplicate switch added in cmdline
        print "\n############\nrmdup switch detected, duplicate removal included in analysis\n############\n\n";
        print "Starting removing duplicates and creating new BAM files..\n";
    }
    #Read BAM files
    print "Reading BAM files to process..\n";
    my @inputfiles = readFile($params->{inputdir}, ".bam");  
    #Start analysis from BAM file
    startWithBam(\@inputfiles);
}

sub startWithAvgCountMode{
    #Read count TXT files
    print "Reading count files..\n";
    my @inputfiles = readFile($params->{inputdir}, ".txt");
    print "Starting counts analysis..\n";
    startWithAvgCount(\@inputfiles);
}

sub startWithMatchScoreMode{
    if (defined $params->{sexchr}) { #Use sex chromosomes switch added in cmdline
        print "\n\n############\nsexchr switch detected, sex chromosomes included in analysis\n############\n\n";
    }
    #Read count TXT files
    print "Starting search for best match scores..\n";
    print "Reading count files..\n";
    my @inputfiles   = readFile($params->{inputdir}, ".normalized.coverage.txt"); #Read files in input directory
    my @controlfiles = readFile($params->{controlsdir}, ".normalized.coverage.txt"); #Read files in controls directory
    print "Starting match score analysis..\n";
    #Start analysis from match score file
    startWithMatchScore( ".normalized.coverage.txt", \@inputfiles, \@controlfiles);
}

sub startWithBestScoreMode{
    #continue
    if (defined $params->{sexchr}) { #Use sex chromosomes switch added in cmdline
        print "\n\n############\nsexchr switch detected, sex chromosomes included in analysis\n############\n\n";
    }
    #Read count TXT files
    print "Starting search for best scores..\n";
    #Read all normalized autosomal coverage control files into array
    my @normAutoControls = readFile($params->{inputdir}, ".normalized.autosomal.coverage.all.controls.txt"); #Read files in input directory
    #Read best match score files
    print "Reading best match score files..\n";
    my @inputfiles = readFile($params->{inputdir}, ".best.match.score.txt"); #Read files in input directory
    print "Starting best score analysis..\n";
    #Start CNV detection analysis
    startWithBestScore(".best.match.score.txt", \@inputfiles, \@normAutoControls);
}

sub generateTargetQcListMode{
    #Read count TXT files
    print "Reading controls directory..\n";
    #Read all *.txt files in control directory into array
    my @inputfiles   = readFile($params->{inputdir}, ".normalized.coverage.txt");    #Read files in input directory
    my @controlfiles = readFile($params->{controlsdir}, ".normalized.coverage.txt"); #Read files in controls directory
    
    #Set temporary dirs to do intermediate work
    my $tmpStartWithMatchScore = File::Temp->newdir();
    my $tmpStartWithBestScore  = File::Temp->newdir();

    $params->{outputdir} = $tmpStartWithMatchScore->dirname;
    
    #Start analysis from match score file
    startWithMatchScore(".normalized.coverage.txt", \@inputfiles, \@controlfiles);
    
    #Read count TXT files
    print "Starting search for best scores..\n";
    #Read all normalized autosomal coverage control files into array
    my @normAutoControls = readFile($tmpStartWithMatchScore->dirname, ".normalized.autosomal.coverage.all.controls.txt"); #Read files in input directory
    #Read best match score files
    print "Reading best match score files..\n";
    @inputfiles = readFile($tmpStartWithMatchScore->dirname, ".best.match.score.txt"); #Read files in input directory
    print STDERR "Starting best score analysis..\n";
    
    $params->{inputdir}  = $tmpStartWithMatchScore->dirname;
    $params->{outputdir} = $tmpStartWithBestScore->dirname;
    
    #Start CNV detection analysis
    startWithBestScore(".normalized.coverage.txt", \@inputfiles, \@normAutoControls);
    
    
    #Read all *.log files in control directory into array
    $params->{inputdir} = $tmpStartWithBestScore->dirname;
    my @logfiles = readFile($tmpStartWithBestScore->dirname, ".log"); #Read files in controls directory
    
    #Extract sampleratio from *.log file

    print STDERR "\n#######################################\n";
    print STDERR "\n\nSamples failing sample CV threshold of ".$params->{sampleRatioScore}.":\n\n";
    
    my @passSampleRatioSamples;
    foreach my $logfile (@logfiles){
        my $dirname= $tmpStartWithBestScore->dirname;
        my $grep = `grep SAMPLE_CV: $dirname/$logfile`;
        chomp $grep;
        if ($grep =~ m/SAMPLE_CV: ([0-9].+)/gs){ #Check if sampleratio is below threshold(default 0.09), otherwise don't use it in targetlist analysis
            my $sampleRatio = $1;
            if ($sampleRatio <= $params->{sampleRatioScore}) {
                $logfile =~ s/.log/.totallist.txt/gs; #Change *.log extension to *.totallist.txt
                push(@passSampleRatioSamples, $dirname."/".$logfile);
            }else{
                #Print samples failing sample CV score to stdout
                print "$logfile\n";
            }
        }
    }
    
    print "\n\n#######################################\n\n";
    print "\nGenerating target QC list..\n";
    
        
    #Calculate target autoVC
    generateTargetQcList(@passSampleRatioSamples);
    print "\nDone generating target QC list\n";
}


sub createFinalListMode{
    #Read shortlist TXT files
    print "Reading input directory..\n";
    #Read all *.txt files in input directory into array
    my @inputfiles = readFile($params->{inputdir}, ".shortlist.txt"); #Read files in controls directory
    
    #Generate the final list
    createFinalList(\@inputfiles);
}           

###############################################
## Main code to generate final list using    ##
## the target QC list                        ##
sub createFinalList{
    my ($inputfiles) = @_;                                              
    #Read target QC list into array;
    open(FILE, $params->{targetQcList}) or die("Unable to open file: $!"); #Read count file
    my @file= <FILE>;
    close(FILE);
    #Retrieve chr, start, stop and genename from file by searching column indices
    my $header = uc($file[0]); #Header in uppercase
    chomp $header;
    #Extract columns from header and use as index
    my @colNames = qw(CHR START STOP GENE);
    my @indices = getColumnIdx($header, \@colNames);
    my $chrIdx = $indices[0];
    my $startIdx = $indices[1];
    my $stopIdx = $indices[2];
    my $geneIdx = $indices[3];
    my $lastFileIdx=$#file;
    
    #Push every index into hash, to exclude from analysis later
    my %toExclude;
    foreach my $index (@indices){
        $toExclude{ $index } = $index;
    }
    
    my @indicesCalculation;
    #Iterate over header to extract indices to calculate 20% quality over
    my @Header = split("\t", $header);
    my $lastHeaderIdx = $#Header;
    for (my $j=0; $j <= $lastHeaderIdx; $j++){
        if (exists $toExclude{ $j }){
            #Exclude from calculation
        }else{
            push(@indicesCalculation, $j);
        }
    }
    
    my %targetQC;
    #Retrieve chr, start, stop, gene for each line in target QC list
    for (my $i=1; $i<=$lastFileIdx; $i++){ #Iterate over lines in file
        my $line=$file[$i];
        chomp $line;
        my @lines=split("\t",$line);
        my $chr = $lines[$chrIdx];
        my $start = $lines[$startIdx];
        my $stop = $lines[$stopIdx];
        my $gene = $lines[$geneIdx];
        #Calculate number of samples passing cutoff, setting target to PASS or FAIL
        my $key = "$chr\t$start\t$stop\t$gene";
        my $total=0;
        my $highQual=0;
        my $lowQual=0;
        foreach my $element (@indicesCalculation){
            my $autoVC = $lines[$element];
            $total++;
            if ($autoVC eq "NA"){
                $lowQual++;
            }elsif ($autoVC <= 0.10) { #Good quality
                $highQual++;
            }else{ #Low quality target
                $lowQual++;
            }
        }
        my $perc = (($lowQual/$total)*100); #percentage low quality targets
        if ($perc > $params->{percentageLessReliableTargets}){ #More than X percentage (default 20%) of targets are low quality, so FAIL
            $targetQC{ $key } = "FAIL";
        }else{
            $targetQC{ $key } = "PASS";
        }
    }

    #Iterate over input files, asses if failing targets are within the calls in the shortlist file
    foreach my $shortlist (@$inputfiles){
        print STDERR"\n\n#######################################\n";
        print STDERR"Analyzing sample: $shortlist..\n";
        #Read file into array;
        my ($file,$dir,$ext) = fileparse($shortlist, qr/\.[^.]*/);
        $file =~ s/.shortlist//gs;
        open(FILE, "$params->{inputdir}/$shortlist") or die("Unable to open file: $!"); #Read count file
        my @file= <FILE>;
        close(FILE);
        #Retrieve chr, start, stop, genename and region coverage from file by searching column indices
        my $header = uc($file[0]); #Header in uppercase
        chomp $header;
        my $outputToWrite = "$header\n";
        my $outputfileBedToWrite = "#gfftags\n";
        my @colNames = qw(CHR START STOP GENE);
        my $indices  = getColumnIdx($header, \@colNames);
        my $chrIdx = $indices[0];
        my $startIdx = $indices[1];
        my $stopIdx = $indices[2];
        my $geneIdx = $indices[3];
        my $lastFileIdx=$#file;
        #Retrieve chr, start, stop, gene and regcov for each line in avg count file
        for (my $k=1; $k<=$lastFileIdx; $k++){
            my $line=$file[$k];
            chomp $line;
            my @lines=split("\t",$line);
            my $chr=$lines[$chrIdx];
            $chr =~ s/X/23/g;
            $chr =~ s/Y/24/g;
            $chr =~ s/Z/25/g;
            my $start=$lines[$startIdx];
            my $stop=$lines[$stopIdx];
            my $gene=$lines[$geneIdx];
            
            my $totalTargets=0;
            my $targetsFail=0;
            #Check if all targets within a single event fail, if true filter event out of finallist.txt
            foreach my $target (sort keys %targetQC){
                my @targets = split("\t", $target);
                my $chrQC = $targets[0];
                $chrQC =~ s/X/23/g;
                $chrQC =~ s/Y/24/g;
                $chrQC =~ s/Z/25/g;
                my $startQC = $targets[1];
                my $stopQC = $targets[2];
                my $geneQC = $targets[3];
                if ($chrQC == $chr && $startQC >= $start && $stopQC <= $stop) { #chromosomes matching and start and stop of targetQC falls within the region extracted from *.shortlist.txt file
                    $totalTargets++;
                    #Check if target is pass or fail
                    my $check = $targetQC{ $target };
                    if ($check eq "FAIL") { #Target failed QC, count it
                        $targetsFail++;
                    }
                }
            }
            #Check if number of failing targets equals number of total targets within this event, if true event fails quality threshold and is removed from *.finallist.txt
            if ($totalTargets == $targetsFail) {
                #event fails QC, don't write it to output
                print STDERR"\nEvent failing target QC: $line\n";
            }else{
                $outputToWrite .= "$line\n";
                my @fields = split /\t/, $line;
                my $chr = shift @fields;
                my $start = shift @fields;
                my $end = shift @fields;
                my $gene = shift @fields;
                my $gene_targets = shift @fields;
                my $n_targets = shift @fields;
                my $n_targets_SHAPIRO = shift @fields;
                my $abberation = shift @fields;

                my $namestring =  "Name=".$abberation.":".$gene_targets.";".
                                  "Note=NUMBER_OF_TARGETS:".$n_targets.",".
                                  "NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST:".$n_targets_SHAPIRO.";";
                $outputfileBedToWrite .= join "\t", $chr,
                                                    $start,
                                                    $end,
                                                    $gene,
                                                    $namestring."\n";
            }
        }
        #Write output to *.finallist.txt file
        my $outputfile = $params->{outputdir}."/".$file.".finallist.txt"; #Output filename
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        my $nr_of_lines = ($outputfileBedToWrite =~ tr/\n//);
        if ($nr_of_lines > 1) {
            my $outputfileBed = $params->{outputdir}."/".$file.".finallist.bed"; #Output filename
           $outputfileBedToWrite =~ s/ /%20/g;
            writeOutput($outputfileBed, $outputfileBedToWrite); #Write output to above specified file
        }
        print STDERR "#######################################\n\n";
        undeff($outputToWrite);
    }
}
###############################################
## Main code to generate target QC list from ##
## files in controls directory               ##
sub generateTargetQcList{
    my @passSampleRatioSamples = @_;

    #Set header for output file
    my $outputToWrite = "CHR\tSTART\tSTOP\tGENE"; #Instantiate header to write
    my %autoVcHash;
    my @keyFiles;
    my $lastInputfilesIdx = $#passSampleRatioSamples; #Idx of last inputfile
    my @autoVCsToOutput;
    for (my $j=0; $j <= $lastInputfilesIdx; $j++){ #Iterate over inputfiles
        my $totallist = $passSampleRatioSamples[$j];
        #Read file into array;
        my ($file,$dir,$ext) = fileparse($totallist, qr/\.[^.]*/);
        #Append filename to outputToWrite to create header
        $outputToWrite .= "\t$file";
        open(FILE, "$totallist") or die("Unable to open file: $!"); #Read count file
        my @file= <FILE>;
        close(FILE);
        #Retrieve chr, start, stop, genename and region coverage from file by searching column indices
        my $header = uc($file[0]); #Header in uppercase
        chomp $header;
        #Extract columns from header and use as index
        my @colNames = qw(CHR START STOP GENE AUTO_VC);
        my @indices = getColumnIdx($header, \@colNames);
        my $chrIdx = $indices[0];
        my $startIdx = $indices[1];
        my $stopIdx = $indices[2];
        my $geneIdx = $indices[3];
        my $autoVcIdx = $indices[4];
        my $lastFileIdx=$#file;
        #Retrieve chr, start, stop, gene and regcov for each line in avg count file
        my @autoVCs;
        for (my $i=1; $i<=$lastFileIdx; $i++){ #Iterate over lines in file
            my $line=$file[$i];
            chomp $line;
            my @lines=split("\t",$line);
            my $chr = $lines[$chrIdx];
            my $start = $lines[$startIdx];
            my $stop = $lines[$stopIdx];
            my $gene = $lines[$geneIdx];
            my $autoVc = $lines[$autoVcIdx];
            if ($j == 0){ #If first file to pass, add chr,start,stop,gene to array, to use in output later
                push(@autoVCsToOutput, "$chr\t$start\t$stop\t$gene");
            }
            push(@autoVCs, $autoVc); #Push autoVC in array
        }
        my $autoVcString = join("_", @autoVCs); #Join all autoVCs and push into string
        $autoVcHash{ $file } = $autoVcString; #Push filename as key and string as value into hash
        push(@keyFiles, $file); #Save filenames in particular order for output later
    }
    $outputToWrite .= "\n";
    
    my $lastKeyFilesIdx = $#keyFiles; #Retrieve last idx from files
    my @autoVcOutput;
    #Loop through files to create matrix of autoVCs
    for (my $i=0; $i <= $lastKeyFilesIdx; $i++){
        my $key = $keyFiles[$i]; #Loop through filenames, take value from hash, split into array and loop through array per line
        my @autoVCsToPrint = split("_", $autoVcHash{ $key }); #Push splitted value into array to iterate over later
        
        my $lastAutoVCsToPrintIdx = $#autoVCsToPrint; #Last index from autoVC array
        
        #Loop through all values and add them to original value in autoVCsToOutput array
        for( my $k=0; $k<=$lastAutoVCsToPrintIdx; $k++){ #Iterate over all autoVC values for sample
            my $autoVcToAppend = $autoVCsToPrint[$k]; #Retrieve autoVC for target k from sample
            my $newElement;
            if ($i == $lastKeyFilesIdx) { #last file to process, so add additional linebreak to value
                $newElement = $autoVCsToOutput[$k] . "\t$autoVcToAppend\n"; #Take original value from array and append new autoVC for target k to it
            }else{
                $newElement = $autoVCsToOutput[$k] . "\t$autoVcToAppend";
            }
            #Splice output element k from array and update with new (above created) element
            splice(@autoVCsToOutput, $k, 1, $newElement);
        }
    }

    #Concat array with autoVCs to output string
    my $outputString = join("", @autoVCsToOutput);
    $outputToWrite .= $outputString;
    
    
     #reset outputs and inputs to original folders
    $params->{outputdir} = $params->{outputdirOriginal};
    $params->{inputdir} = $params->{inputdirOriginal};   
    
    #Write output to file
    my $outputfile = $params->{targetQcList}; #Output filename
    writeOutput($outputfile, $outputToWrite); #Write output to above specified file
    

}                                                                                                                                     



###############################################
## Main code to extract region coverage from ##
## BAM file(s)                               ##
sub startWithBam{
    my ($inputfiles) = @_;                     
    ###############################################
    ## Main code to extract region coverage from ##
    ## BAM file(s)                               ##
    foreach my $bam (@$inputfiles){
    
        #Check if *.bam.bai or *.bai file exist, otherwise skip this bam file
        my ($file,$dir,$ext) = fileparse($bam, qr/\.[^.]*/);
        print STDERR "Processing $file..\n";
        my $bai_file = $file.".bai";
        my $file_to_count;
        unless (-e $params->{inputdir}."/".$bai_file) {
            #Don't process this file, because it doesn't have an indexfile
            print  STDERR  "##### WARNING #####WARNING #####\n".
                            "Cannot find an index file for file: ".$params->{inputdir}."/".$bam."\n".
                            "\t,skipping this file from analysis\n".
                            "##### WARNING ##### WARNING #####\n\n";
        }else{
            #set temp files... if necessary 
            my $tmp_dir     = File::Temp->newdir();
            my $rmdup_bam   = File::Temp->new( TEMPLATE => 'tempXXXXX',DIR => $tmp_dir, SUFFIX => '.rmdup.bam'  );

            #Check if duplicates need to be removed
            if (defined $params->{rmdup}){
                #Process BAM files generating duplicate removed BAM files
                rmDupBam($bam, $rmdup_bam, $tmp_dir);
                print STDERR "Starting counts analysis..\n"; #Start to count regions
                $file_to_count = $rmdup_bam;
            }else{ #BAM files are already rmdupped, add them to list of files to process (retrieve them from inputdir cmdline)
                $file_to_count = $params->{inputdir}."/".$bam;
            }
            
            print STDERR "Starting counts analysis..\n";
            countFromBam($file_to_count, $file);
        }
    }    
}    
sub startWithAvgCount{
    my ($inputfiles) = @_;
     
    foreach my $txt (@$inputfiles){
        #Set header for output file
        my $outputToWrite = "CHR\tSTART\tSTOP\tGENE\tTARGET\tREGION_COV\tAVG_AUTOSOMAL_COV\tAVG_TOTAL_COV\tAVG_GENE_COV\tNORMALIZED_AUTOSOMAL\tNORMALIZED_TOTAL\tNORMALIZED_GENE\n";
        #Read file into array;
        my ($file,$dir,$ext) = fileparse($txt, qr/\.[^.]*/);
        open(FILE, $params->{inputdir}."/".$txt) or die("Unable to open file: $!"); #Read count file
        my @file= <FILE>;
        close(FILE);
        #Retrieve chr, start, stop, genename and region coverage from file by searching column indices
        my $header = uc($file[0]); #Header in uppercase
        chomp $header;
        my @colNames = qw(CHR START STOP GENE TARGET REGION_COV);
        my @indices = getColumnIdx($header, \@colNames);
        my $args;
        $args->{chrIdx} = $indices[0];
        $args->{startIdx} = $indices[1];
        $args->{stopIdx} = $indices[2];
        $args->{geneIdx} = $indices[3];
        $args->{targetIdx} = (defined $indices[4] ?  $indices[4] : undef);
        $args->{regcovIdx} = $indices[5];
        $args->{lastFileIdx} =$#file;
        #Calculate coverage on gene
        my ($coverage, $genehash, $genes, $covchrauto, $covchrsex) = calcGeneCov($args, @file);

        #Calculate coverage including sex chromomsomes
        my ($covchrautoval, $covchrall, $counts) = calcCovAutoSex($genes, $covchrauto, $covchrsex);
        
        #Foreach line in input file write away all calculated stats/values
        $outputToWrite .= writeCountFile($args, $covchrautoval, $covchrall, $genehash, $counts, $coverage, @file);
        my $outputfile = $params->{outputdir}."/".$file.".normalized.coverage.txt"; #Output filename
        print STDERR "Writing normalized coverage counts to: $outputfile\n\n\n";
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        print STDERR "Finished processing file: $txt\n\n\n";
    }
}  
  
  
##########################################################
## Main code to select most informative control samples ##
##########################################################
sub startWithMatchScore{
    my $extension = shift;
    my ($inputfiles) = $_[0];
    my ($controlfiles) = $_[1];
    foreach my $inputfile (@$inputfiles){ #Open sample file
        my $outputExtnsn = "normalized.autosomal.coverage.all.controls.txt"; #Specify extension for output normalized coverage file
        my $colsToExtract = "CHR START STOP GENE TARGET REGION_COV NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL"; #Specify Columns to extract
        
        my ($autodiff, $sexdiff, $perfectMatch) = createNormalizedCoverageFiles($inputfile, $extension, $outputExtnsn, $colsToExtract, \@$controlfiles); #Give input file to analyze, extension, columns to extract and list of controlfiles to use to function to retrieve normalized coverage
        
        my %autodiffh = %$autodiff;
        my %sexdiffh = %$sexdiff;

        #Open output bestmatch file
        my ($outputPostfixRemoved,$dir,$ext) = fileparse($inputfile, $extension);

        #$outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
        my $outputfile = $params->{outputdir}."/".$outputPostfixRemoved.".best.match.score.txt"; #Output filename
        #Do additional checks to detect perfect match samples and check when this occurs if there still are enough samples available in controlsdirectory to continues analysis
        my $numBestMatchSamplesToProcess = keys %autodiffh; #Count number of samples to process (all - perfectMatch)
        my $numBestMatchSamples=$params->{numBestMatchSamplesCmdL};
        if ( $perfectMatch > 0){
            if ($numBestMatchSamplesToProcess == $numBestMatchSamples){ #Detected perfect matches, but number of samples to process is the same as number of samples requested on cmdline
                print STDERR "\n##### WARNING ##### WARNING #####\n";
                print STDERR "Detected $perfectMatch perfect match(es) between sample and controlsamples, excluding these samples from the analysis.";
                print STDERR "\n##### WARNING ##### WARNING #####\n";
            }elsif ($numBestMatchSamplesToProcess >= ($numBestMatchSamples+$perfectMatch)){ #If there are more or equal #samples in controlsdir then specified on cmdline, use best #cmdline samples
                $numBestMatchSamples=$numBestMatchSamples;
            }else{ #Less samples than cmdline specified available, continue analysis but throw warning
                print STDERR "\n##### WARNING ##### WARNING #####\n";
                print STDERR "Detected $perfectMatch perfect match(es) between sample and controlsamples, which means only $numBestMatchSamplesToProcess instead of $numBestMatchSamples samples from controls directory are used for Match score analysis.";
                print STDERR "\n##### WARNING ##### WARNING #####\n";
                $numBestMatchSamples = $numBestMatchSamplesToProcess;
            }
        }
        $perfectMatch=0; #Reset perfectMatch variable
        print STDERR "\n#################################################################\n";
        print STDERR "Selecting best $numBestMatchSamples control samples for analysis..\n";
        print STDERR "###################################################################\n";
        my @keys;
        my @vals;
        if (defined $params->{sexchr}) { #Use sex chromosomes switch, print all abs diffs including sex chromosomes
            @keys = sort { $sexdiffh{$a} <=> $sexdiffh{$b} } keys(%sexdiffh);
            @vals = @sexdiffh{@keys};
        }else { #Only use autosomal chrs
            @keys = sort { $autodiffh{$a} <=> $autodiffh{$b} } keys(%autodiffh);
            @vals = @autodiffh{@keys};
        }
        my $outputToWrite= "SAMPLE\tSAMPLE_PATH\tCONTROL_SAMPLE\tCONTROL_SAMPLE_PATH\tAVERAGE_BEST_MATCH_SCORE\n"; #Assign header to output best match file
        
        for (my $k=0; $k < $numBestMatchSamples; $k++){
            print STDERR "Control: " . $keys[$k] . "\t\t\tAvg abs diff: " . $vals[$k] . "\n";
            my $lin = join "\t",    $inputfile,
                                    $params->{inputdir}."/".$inputfile,$keys[$k],
                                    $params->{controlsdir}."/".$keys[$k],$vals[$k]."\n";
            $outputToWrite .= $lin; #concatenate full generated line to files
        }
        print STDERR "#######################################\n\n";
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
    }
}

###############################################
## Main code to detect CNVs                  ##
###############################################
sub startWithBestScore{
    my $extension = shift;
    my ($inputfiles_ref) = shift;
    my ($normAutoControls_ref) =  shift;
    
    my @inputfiles       = @$inputfiles_ref;
    my @normAutoControls = @$normAutoControls_ref;
    
    my $n_files = scalar @inputfiles;

    for (my $m=0; $m < $n_files; $m++){
        my $inputfilename = $inputfiles[$m];

        my @sampleFile;
        my @controlFile;
        my @controlChr;
        my @controlStart;
        my @controlStop;
        my @controlGene;
        my @controlTarget;
        my @sampleRatio;
        print STDERR "\nAnalyzing sample: $inputfilename..\n";
        
        my $outputToWrite= "CHR\tSTART\tSTOP\tGENE\tTARGET\t$inputfilename\t"; #Assign header to output best match file
        open(INPUTFILE, $params->{inputdir}."/".$inputfilename) or die("Unable to open file: $!"); #Read best match file
        my @inputfiledata= <INPUTFILE>;
        close(INPUTFILE);
        my $header = uc($inputfiledata[0]); #Header in uppercase
        chomp $header;
        
        my @colNames = qw(SAMPLE SAMPLE_PATH CONTROL_SAMPLE CONTROL_SAMPLE_PATH);
        my @indices = getColumnIdx($header, \@colNames);
        my $sampleIdx      = $indices[0];
        my $samplePathIdx  = $indices[1];
        my $controlIdx     = $indices[2];
        my $controlPathIdx = $indices[3];
        my @array = split("\t", $header);
        my $samplename  = $array[$sampleIdx];
        my $samplefile  = $array[$samplePathIdx];
        my $controlname = $array[$controlIdx];
        my $controlfile = $array[$controlPathIdx];
        
        my $lastLine=$#inputfiledata;
        
        #Initialize arrays for arrayreferences
        my @NORMAUTOSOMAL;
        my @NORMTOTAL;
        my @NORMGENE;
        
        # Read samplefile into array (only do this once)
        my $currentline = $inputfiledata[1]; #Open first samplefile from inputfile
        my @crntlnArray = split("\t", $currentline);
        # Read input samplefile into array
        open(SAMPLEFILE, $crntlnArray[$samplePathIdx]) or die("Unable to open file: $!"); #Read SAMPLE best match file
        @sampleFile= <SAMPLEFILE>;
        close(SAMPLEFILE);
        my $headerSample = uc($sampleFile[0]); #Header in uppercase
        chomp $headerSample;
        
        # Read normalized_autosomal values file
        #Open norm auto file, further on calculate if sample is within 3SD from mean, otherwise exclude in sampleratio calculation
        my @normAutoCon = @normAutoControls;
        my $normautofile = $normAutoCon[$m];
        open(NORMAUTOFILE, $params->{inputdir}."/".$normautofile) or die("Unable to open file: $!"); #Read count file
        my @normautofile= <NORMAUTOFILE>;
        close(NORMAUTOFILE);
        
        #Extract sample file header indices
        @colNames = qw(CHR START STOP GENE TARGET NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL NORMALIZED_GENE);
        my @indices_control = getColumnIdx($headerSample, \@colNames);
        my $chrIdxControl     = $indices_control[0];
        my $startIdxControl   = $indices_control[1];
        my $stopIdxControl    = $indices_control[2];
        my $geneIdxControl    = $indices_control[3];
        my $targetIdxControl  = $indices_control[4];
        my $normAutoIdxSample = $indices_control[5];
        my $normSexIdxSample  = $indices_control[6];
        my $normGeneIdxSample = $indices_control[7];
        #Read through samplefile and push all three values into array of arrays;
        my $lastLineSampleFile = $#sampleFile;
        my @sampleNormAuto;
        my @sampleNormTotal;
        my @sampleNormGene;

        for (my $i=1; $i<=$lastLineSampleFile; $i++){
            my $currentline = $sampleFile[$i];
            my @crntlnArray = split("\t", $currentline);
            push(@controlChr     , $crntlnArray[$chrIdxControl]    );
            push(@controlStart   , $crntlnArray[$startIdxControl]  );
            push(@controlStop    , $crntlnArray[$stopIdxControl]   );
            push(@controlGene    , $crntlnArray[$geneIdxControl]   );
            push(@controlTarget  , (defined $targetIdxControl ? $crntlnArray[$targetIdxControl] : "-" ));
            push(@sampleNormAuto , $crntlnArray[$normAutoIdxSample]);
            push(@sampleNormTotal, $crntlnArray[$normSexIdxSample] );
            push(@sampleNormGene , $crntlnArray[$normGeneIdxSample]);
        }
        my $refNormAuto  = \@sampleNormAuto; #Create array references
        my $refNormTotal = \@sampleNormTotal;
        my $refNormGene  = \@sampleNormGene;
        push(@NORMAUTOSOMAL, $refNormAuto);
        push(@NORMTOTAL, $refNormTotal);
        push(@NORMGENE, $refNormGene);
        #Finalize header line for outputfile
        my $lin = "AUTO_RATIO\tAUTO_ZSCORE\tAUTO_VC\tGENE_RATIO\tGENE_ZSCORE\tGENE_VC\tSHAPIRO-WILK\n";
        $outputToWrite .= $lin; #concatenate full generated line to files
        
        #Iterate over controlfiles and put into array of arrays
        for (my $i=1; $i<=$lastLine; $i++){
            # Read line from inputfile
            my $currentline = $inputfiledata[$i];
            my @crntlnArray = split("\t", $currentline);
            #Read control file into array
            open(CONTROLFILE, $crntlnArray[$controlPathIdx]) or die("Unable to open file: $!"); #Read CONTROL best match file
            @controlFile= <CONTROLFILE>;
            close(CONTROLFILE);
            #Control file is now in an array
            my $headerControl = uc($controlFile[0]); #Header in uppercase
            chomp $headerControl;
            #Extract control header columnindices
            my @colNames = qw(NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL NORMALIZED_GENE);
            my @indices_control = getColumnIdx($headerControl, \@colNames);
            my $normAutoIdxControl  = $indices_control[0];
            my $normSexIdxControl   = $indices_control[1];
            my $normGeneIdxControl  = $indices_control[2];
            
            my $lastLineControlFile = $#controlFile;
            #Walk through contents of control file
            my @controlNormAuto;
            my @controlNormTotal;
            my @controlNormGene;
            for (my $i=1; $i<=$lastLineControlFile; $i++){
                my $currentline = $controlFile[$i];
                my @crntlnArray = split("\t", $currentline);
                push(@controlNormAuto , $crntlnArray[$normAutoIdxControl]);
                push(@controlNormTotal, $crntlnArray[$normSexIdxControl]);
                push(@controlNormGene , $crntlnArray[$normGeneIdxControl]);
            }
            my $refNormAuto = \@controlNormAuto; #Create array references
            my $refNormTotal = \@controlNormTotal;
            my $refNormGene = \@controlNormGene;
            push(@NORMAUTOSOMAL, $refNormAuto);
            push(@NORMTOTAL, $refNormTotal);
            push(@NORMGENE, $refNormGene);
        }

        #Iterate through controls in array of arrays
        for (my $i=0; $i<$lastLineSampleFile; $i++){ #$numBestMatchSamples
            my @controlAutoArray;
            my @controlGeneArray;
            my @controlInputArray = @NORMAUTOSOMAL; #Use default NORMAUTOSOMAL array
            my $sampleValue;
            #Iterate through all normalized values
            if (defined $params->{sexchr}) { #Use sex chromosomes switch, so NORMTOTAL needs to be used
                @controlInputArray = @NORMTOTAL;
            }
            #samplevalue is always first array (index 0) to interate over
            $sampleValue = $controlInputArray[0]->[$i]; #Use sampleValue normalized on autosomal (or total) controls
            for (my $k=1; $k<=$lastLine; $k++){ #$lastLineSampleFile
                #Push all control values in array
                push (@controlAutoArray, $controlInputArray[$k]->[$i]);
                push (@controlGeneArray, $NORMGENE[$k]->[$i]);
            }
            
            #####SAMPLE NORMALIZATION#####
            my $lin = $controlChr[$i] . "\t" . $controlStart[$i] . "\t" . $controlStop[$i] . "\t" . $controlGene[$i] . "\t" . $controlTarget[$i];
            
            my ($autoMean, $autoSD, $autoRatio, $autoZscore, $autoVc, $sampleR) = calcAutoRatioZscoreVc($sampleValue, $lin, \@controlAutoArray);
            push @sampleRatio, $sampleR;
            $lin .= "\t$sampleValue\t$autoRatio\t$autoZscore\t$autoVc";
            
            my $TNsi = $NORMGENE[0]->[$i];
            my $TNciMean = $autoMean;
            my $TNciSD=  $autoSD;
            my $ATNci = $TNciMean;
            my $TNNsi;
            if ($ATNci == 0) { #Check if normalized total depth for target in sample equals 0 (In very rare occasions this is true)
                $TNNsi = 0;
            }else{
                $TNNsi = ($TNsi/$ATNci);
            }
            
            #####GENE NORMALIZATION#####
            $sampleValue = $NORMGENE[0]->[$i]; #Use sampleValue normalized gene
            #Calculate mean and SD for Gene array
            my ($geneMean, $geneSD) = calcMeanSD(\@controlGeneArray);
            #calculate ratio, z-score and variation coefficient
            #ratio, observed devided by mean
            my $geneRatio;
            my $geneZscore;
            my $geneVc;
            if ($geneMean == 0 || $geneSD == 0) {
                $geneRatio = "NA";
                $geneZscore = "NA";
                $geneVc = "NA";
                print STDERR "\n##### WARNING ##### WARNING #####\n";
                print STDERR $controlChr[$i] . ":" . $controlStart[$i] . "-" . $controlStop[$i] . "\t" . $controlGene[$i] . "\t" . $controlTarget[$i]."\n";
                print STDERR "Gene normalization route not available.\nMean or Standard Deviation is 0.\nCan not calculate ratio, zscore and variation coefficient on this region\n";
            }else{
                $geneRatio = ($sampleValue/$geneMean);
                #z-score, observed minus mean devided by sd
                $geneZscore = (($sampleValue-$geneMean)/$geneSD);
                #variation coefficient, sd devided by mean
                $geneVc = ($geneSD/$geneMean);
            }
            
            #####SHAPIRO-WILK TEST#####
            #Calculate z-score foreach control sample and perform shapiro-wilk test
            my @zscores;
            my $pVal;
            if ($autoMean == 0 || $autoSD == 0) {
                $pVal = "NA";
            }else{
                foreach my $controlRatio (@controlAutoArray){
                    my $Zscore = (($controlRatio-$autoMean)/$autoSD);
                    push(@zscores, $Zscore);
                }
                #Do shapiro-wilk test
                $pVal = shapiro_wilk_test([@zscores]);
            }
            
            #Concatenate outputs
            $lin .= "\t$geneRatio\t$geneZscore\t$geneVc\t$pVal\n";
            $outputToWrite .= $lin; #concatenate full generated line to files
            
            undef($autoMean); undef($autoSD); undef($autoRatio); undef($autoZscore); undef($autoVc);
        }
        #Open output bestmatch file
        my ($outputPostfixRemoved,$dir,$ext) = fileparse($inputfilename, $extension);
        #my $outputPostfixRemoved = $inputfile;
        #$outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
        my $outputfile = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.txt"; #Output filename
        print STDERR "#######################################\n\n";
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        
        #Write sample CV score to log file in output directory
        $outputfile = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.log"; #Output filename
        my $lastLineIdx = $#inputfiledata;
        $header = uc($inputfiledata[0]); #Header in uppercase
        chomp $header;
        #Extract control header columnindices
        @colNames = qw(AVERAGE_BEST_MATCH_SCORE);
        my @ind = getColumnIdx($header, \@colNames);
        my $avgBestMatchScoreValIdx = $ind[0];
        my @avgBestMatchScores;
        for (my $i=1; $i<=$lastLineIdx; $i++){
            my $line = $inputfiledata[$i];
            chomp $line;
            my @array = split("\t", $line);
            my $avgBestMatchScoreVal = $array[$avgBestMatchScoreValIdx];
            push(@avgBestMatchScores, $avgBestMatchScoreVal);
            $outputToWrite .= "$line\n"; #concatenate full generated line to files
        }
        
        ######
        my $failedRegionsToWrite = "\n\n###REGIONS FAILING USER SPECIFIED QUALITY THRESHOLD OF $params->{regionThreshold} PERCENT###\n###THESE REGIONS ARE OMMITTED FROM SAMPLE_CV CALCULATION###\n";
        my $lastLin=$#normautofile; #lastlineIdx of file
        my @idxToKeep;
        for (my $i=1; $i<=$lastLin; $i++){ #Loop through norm auto file
            my $currentline = $normautofile[$i];
            chomp $currentline;
            my @crntlnArray = split("\t", $currentline); #Complete line containing region, sample and control values
            my $region = $crntlnArray[0];
            shift(@crntlnArray); #Remove first element (in this case region) from array
            
            my ($mean, $sd) = calcMeanSD(\@crntlnArray); #returns mean and sd
            # Iterate over elements in line
            my $covPass=0;
            my $covFail=0;
            my $lastCrntlnArray = $#crntlnArray; #retrieve last element Idx from array
            for (my $k=0; $k<=$lastCrntlnArray; $k++){ # $k is norm_auto per sample
                my $currentVal = $crntlnArray[$k]; # Where k=0 is sample, others are controls
                chomp $currentVal;
                if ($mean != 0 && $sd != 0) { #Check if returned mean and SD are not zero, if they are count it as failed coverage too
                    my $toCheck = ($currentVal/$mean);
                    my $threeSD = (3*$sd);
                    my $lowLimit = ($mean-$threeSD);
                    my $upLimit = ($mean+$threeSD);
                    if ($toCheck >= $lowLimit && $toCheck <= $upLimit) { #Value within 3SD from mean
                        $covPass++;
                    }else{
                        $covFail++;
                    }
                }else{
                    $covFail++;
                }
            }
            my $percentageFail;
            if ($covPass == 0) { #If all targets fail, set percentage fail to 100%, else calculate it
                $percentageFail = 100;
            }else{
                $percentageFail = (($covFail/$covPass) * 100);
            }
            if ($percentageFail >= $params->{regionThreshold}) {
                #Fail, don't count in sample CV
                #Remove this region from sampleRatio array
                my $idxToRm = ($i-1); #Subtract 1, since header is missing now
                $failedRegionsToWrite .= "$region\n";
            }else{ #Pass threshold/percentage
                push(@idxToKeep, ($i-1));
            }
        }
        my @calcSampleRatio;
        foreach my $idxKeep (@idxToKeep){ #Iterate over indices to keep
            if (defined $sampleRatio[$idxKeep]){
                push(@calcSampleRatio, $sampleRatio[$idxKeep]); #Push sample CVs to calculate mean and sd on into new array
            }else{
                print STDERR "Sample index undef\n";
            }
        }
        #Select 95% samples (exclude low and high) for sampleRatio calculation
        my @sortedNormVal = sort { $a <=> $b } @calcSampleRatio;
        my $numValues = scalar(@sortedNormVal);
        my $low25perc = ceil((0.025*$numValues));
        my $high25perc = floor((0.975*$numValues));
        my @sliceNormVal = @sortedNormVal[($low25perc) .. ($high25perc-1)]; #Slice values out of array (why not $low25perc-1??)
        
        #Calculate mean average best match score over all control samples
        my ($mean, $sd) = calcMeanSD(\@avgBestMatchScores);
        my $meanAvgBestMatchScore = $mean;
        #Calculate sample CV
        ($mean, $sd) = calcMeanSD(\@sliceNormVal); #Calculate mean and sd
        my $sampleRatio = ($sd/$mean);
        $lin = "\n\nSAMPLE_CV: $sampleRatio\nMEAN_AVERAGE_BEST_MATCHSCORE: $meanAvgBestMatchScore\n"; #Add mean average best match score and sample CV to output logfile
        $outputToWrite .= $lin;
        $outputToWrite .= $failedRegionsToWrite; #Add failed regions to output logfile
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        print STDERR "Sample CV: $sampleRatio\n";
        print STDERR "Mean average best match score of all control samples: $meanAvgBestMatchScore\n";
        print STDERR "#######################################\n\n";
        
        createOutputLists( $extension, $inputfilename);
    }
    
}

sub calcAutoRatioZscoreVc{
    my $sampleValue = $_[0];
    my $lin = $_[1];
    my ($controlAutoArray) = $_[2];
    my $sampleRatio = 0;
    
    #Calculate mean and SD for Auto array
    my ($autoMean, $autoSD) = calcMeanSD(\@$controlAutoArray);
    #mean and sd are returned by function, calculate ratio, z-score and variation coefficient
    #ratio, observed devided by mean
    my ($autoRatio, $autoZscore, $autoVc);
        
    if ($autoMean == 0 || $autoSD == 0) {
        $autoRatio = "NA";
        $autoZscore = "NA";
        $autoVc = "NA";
        print STDERR "\n##### WARNING ##### WARNING #####\n";
        print STDERR $lin."\n";
        print STDERR "Mean or Standard Deviation is 0.\nCan not calculate ratio, zscore and variation coefficient on this region\n";
    }else{
        $autoRatio = ($sampleValue/$autoMean);
        #z-score, observed minus mean devided by sd
        $autoZscore = (($sampleValue-$autoMean)/$autoSD);
        #variation coefficient, sd devided by mean
        $autoVc = ($autoSD/$autoMean);
        #push autoRatio in array, to calculate the VC ratio for the complete sample
        $sampleRatio = $autoRatio;
    }
    return($autoMean, $autoSD, $autoRatio, $autoZscore, $autoVc, $sampleRatio);
}

sub createOutputLists{
    my $extension = shift;
    my $inputfile = shift;
    my ($outputPostfixRemoved,$dir,$ext) = fileparse($inputfile, $extension);
    #my $outputPostfixRemoved = $inputfile;
    #$outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
    #Open best.score.txt file to extract all targets and mark them
    open(BESTSCORE, $params->{outputdir}."/".$outputPostfixRemoved.".best.score.txt") or die("Unable to open file: $!"); #Read best match file
    my @bestScoreFile= <BESTSCORE>;
    close(BESTSCORE);
    unlink $params->{outputdir}."/".$outputPostfixRemoved.".best.score.txt"; #Remove *.best.score.txt file, since it is almost the same as other files produced
    my @targets;
    my @genes;
    my @genetargets;
    my @values;
    my @autoVCs;
    my @shaps;
    my $header = uc($bestScoreFile[0]); #Header in uppercase
    chomp $header;
    my @colNames = qw(CHR START STOP GENE TARGET AUTO_RATIO AUTO_ZSCORE AUTO_VC GENE_RATIO GENE_ZSCORE GENE_VC SHAPIRO-WILK);
    my @indices = getColumnIdx($header, \@colNames);
    my $chrIdx = $indices[0];
    my $startIdx = $indices[1];
    my $stopIdx = $indices[2];
    my $geneIdx = $indices[3];
    my $targetIdx = $indices[4];
    my $autoRatioIdx = $indices[5];
    my $autoZscoreIdx = $indices[6];
    my $autoVcIdx = $indices[7];
    my $geneRatioIdx = $indices[8];
    my $geneZscoreIdx = $indices[9];
    my $geneVcIdx = $indices[10];
    my $shapIdx = $indices[11];
    my $lastIdx = $#bestScoreFile;
    for (my $i=1; $i <= $lastIdx; $i++){
        my @line = split("\t", $bestScoreFile[$i]);
        my $chr=$line[$chrIdx];
        $chr =~ s/^X/23/;
        $chr =~ s/^Y/24/;
        $chr =~ s/^MT/25/;
        my $start=$line[$startIdx];
        my $stop=$line[$stopIdx];
        my $gene=$line[$geneIdx];
        my $geneTarget=$line[$targetIdx];
        my $autoRatio=$line[$autoRatioIdx];
        my $autoZscore=$line[$autoZscoreIdx];
        my $autoVc=$line[$autoVcIdx];
        my $geneRatio=$line[$geneRatioIdx];
        my $geneZscore=$line[$geneZscoreIdx];
        my $geneVc=$line[$geneVcIdx];
        my $shap=$line[$shapIdx];
        my $target = "$chr\t$start-$stop";
        my $value = "$autoRatio\t$autoZscore\t$autoVc\t$geneRatio\t$geneZscore\t$geneVc";
        chomp $value;
        chomp $shap;
        push(@targets, $target);
        push(@genes, $gene);
        push(@genetargets, $geneTarget);
        push(@values, $value);
        push(@autoVCs, $autoVc);
        push(@shaps, $shap);
    }
    my @abberation;
    my @llValues;
    my $ref1 = \@targets;
    my $ref2 = \@genes;
    my $ref3 = \@values;
    my $ref4 = \@abberation;
    my $ref5 = \@llValues;
    my $ref6 = \@autoVCs;
    my $ref7 = \@shaps;
    my $ref8 = \@genetargets;

    my @arrayRefs = ($ref1, $ref2, $ref3, $ref4, $ref5, $ref6, $ref7, $ref8);
    my $lastTargetIdx = $#targets;
    
    my $currentGene;
    my $currentTarget;
    my $currentValue;
    my @idxs;
    my %geneCounts;
    my %totalGeneCounts;
    my %geneCountsForGeneCorrection;
    my %totalGeneCountsForGeneCorrection;
    #Iterate over targets to generate list with LOW and HIGH targets
    for (my $m=0; $m <= $lastTargetIdx; $m++){
        $currentTarget=$arrayRefs[0][$m];
        $currentGene=$arrayRefs[1][$m];
        my @geneTargetArray;
        my @geneGeneArray;
        my @geneValueArray;
        my $abberationValue = ".";
        my $toCompare = $arrayRefs[1][$m+1]; #Compare with gene belonging to next target
        if ($m==$lastTargetIdx){
           $toCompare = "LAST"; #If last gene
        }
        if ($currentGene eq $toCompare){ #If current gene matches next one, push id of current gene into array
            push(@idxs, $m);
        }else { #Else start with gene analysis
            my $geneAbberationCount=0;
            push(@idxs, $m);
            my $abberationCount = 0;
            my @geneResultLongList;
            my $longListRef = \@geneResultLongList;
            #Do first round of iteration to detect abberations
            foreach my $l (@idxs){ #Foreach elements ID for this gene
                $abberationValue = ".";
                my $numTargets = scalar(@idxs);
                $currentValue=$arrayRefs[2][$l];
                my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
                my $autoratio = $vals[0];
                my $autozscore = $vals[1];
                my $genezscore = $vals[4];
                my $line = $arrayRefs[0][$l] . "\t" . $arrayRefs[1][$l] . "\t" . $arrayRefs[2][$l] . "\t";
                my $firstElement = $idxs[0];
                my $lastidxsIdx = $#idxs;
                my $lastElement = $idxs[$lastidxsIdx];
                
                #Label: if auto ratio, auto zscore and gene zscore don't pass threshold for either negative or positive increase the abberation count.
                if ($autoratio eq "NA" || $autozscore eq "NA" || $genezscore eq "NA") { #One of the values is NA, cannot calculate anything
                    $abberationValue = ".";
                }elsif ($autoratio < $params->{ratioCutOffLow} && $autozscore < $params->{zScoreCutOffLow} && $genezscore < $params->{zScoreCutOffLow}) { #Target labelled as low
                    $abberationValue = "DEL";
                }elsif ($autoratio > $params->{ratioCutOffHigh} && $autozscore > $params->{zScoreCutOffHigh} && $genezscore > $params->{zScoreCutOffHigh}) { #Target labelled as high
                    $abberationValue = "DUP";
                }else{
                    #no aberration
                    $abberationCount=0;
                }
                #Push abberation value belonging to this target in array ref
                $arrayRefs[3]->[$l] = $abberationValue;
            }
            #Do second round of detection to see whether 2 sequential targets are del/dup
            my $lastIdx=$#idxs;
            for (my $i=0; $i <= $lastIdx; $i++){ #For each element in the array do checks
            #foreach my $l (@idxs){ #Foreach elements ID for this gene
                my $l = $idxs[$i]; #Element from array, in this case the id of line
                $abberationValue = ".";
                my $fwdCounter=0;
                my $numTargets = scalar(@idxs);
                $currentValue=$arrayRefs[2][$l];
                my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
                my $autoratio = $vals[0];
                my $autozscore = $vals[1];
                my $genezscore = $vals[4];
                my $currentAbberation = $arrayRefs[3][$l];
                my $nextAbberation;
                if ($l == $idxs[$lastIdx]) { #If next element doesn't exist, create fake one
                    $nextAbberation = "LAST";
                }else{
                    $nextAbberation = $arrayRefs[3][$l+1];
                }
                if ($currentAbberation eq "DEL" || $currentAbberation eq "DUP") { #If DEL or DUP, increase counter
                    #Check next element to see if it also is DEL or DUP
                    if ($nextAbberation eq "DEL" || $nextAbberation eq "DUP") { #Sequential target is also DEL or DUP
                        if ($currentAbberation eq "DEL" && $nextAbberation eq "DEL" || $currentAbberation eq "DUP" && $nextAbberation eq "DUP") { #Check if both targets are DEL or DUP, else we don't count further
                            #Check if autoRatio, autoZscore or geneZscore from previous events meet thresholds
                            for (my $c=$l; $c >= $idxs[0]; $c--){
                                my $currentValue=$arrayRefs[2][$c];
                                my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
                                my $autoratio = $vals[0];
                                my $autozscore = $vals[1];
                                my $genezscore = $vals[4];
                                if ($currentAbberation eq "DEL") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio < $params->{ratioCutOffLow} || $autozscore < $params->{zScoreCutOffLow} || $genezscore < $params->{zScoreCutOffLow}) {
                                            $abberationValue = "DEL";
                                        }else{
                                            $abberationValue = ".";
                                            last;
                                        }
                                    }
                                }
                                if ($currentAbberation eq "DUP") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio > $params->{ratioCutOffHigh} || $autozscore > $params->{zScoreCutOffHigh} || $genezscore > $params->{zScoreCutOffHigh}) {
                                            $abberationValue = "DUP";
                                        }else{
                                            $abberationValue = ".";
                                            last;
                                        }
                                    }
                                }
                                $arrayRefs[3]->[$c] = $abberationValue;
                            }
                            #Check forward autoRatio, autoZscore and geneZscore
                            for (my $d=$l+1; $d <= $idxs[$lastIdx]; $d++){
                                $fwdCounter++;
                                my $currentValue=$arrayRefs[2][$d];
                                my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
                                my $autoratio = $vals[0];
                                my $autozscore = $vals[1];
                                my $genezscore = $vals[4];
                                if ($currentAbberation eq "DEL") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio < $params->{ratioCutOffLow} || $autozscore < $params->{zScoreCutOffLow} || $genezscore < $params->{zScoreCutOffLow}) {
                                            $abberationValue = "DEL";
                                        }else{
                                            $abberationValue = ".";
                                            $d=($idxs[$lastIdx]+1);
                                        }
                                    }
                                }
                                if ($currentAbberation eq "DUP") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio > $params->{ratioCutOffHigh} || $autozscore > $params->{zScoreCutOffHigh} || $genezscore > $params->{zScoreCutOffHigh}) {
                                            $abberationValue = "DUP";
                                        }else{
                                            $abberationValue = ".";
                                            $d=($idxs[$lastIdx]+1);
                                        }
                                    }
                                }
                                $arrayRefs[3]->[$d] = $abberationValue;
                            }
                        }
                    }
                }
            #Set $i to correct value to continue analysis
            $i = ($i+$fwdCounter);
            }
            #Calculate number of abberations detected in gene
            foreach my $element (@idxs){
                my $abberationValue = $arrayRefs[3][$element];
                if ($abberationValue eq "DEL" || $abberationValue eq "DUP" ){
                    $geneAbberationCount++;
                }
            }
            undef(@idxs);
        }
    }
    
    #Print values for all targets
    for (my $m=0; $m <= $lastTargetIdx; $m++){
        my $currentValue=$arrayRefs[2][$m];
        my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
        my $autoratio = $vals[0];
        my $autozscore = $vals[1];
        if ($autoratio ne "NA" && $autozscore ne "NA"){
            $totalGeneCountsForGeneCorrection{$arrayRefs[1][$m]}++; #Count all occurences of genes to use in later calculation, while auto ratio and zscore are not 0. When they are 0 target is not counted
            if ($autoratio < $params->{ratioCutOffLow} || $autoratio > $params->{ratioCutOffHigh}){ #Count if auto ratio is either LOW or HIGH
                $geneCountsForGeneCorrection{$arrayRefs[1][$m]}++; #Abberation detected, count it
            }
        }
    }
    
    #Check auto VC for quality
    for (my $m=0; $m <= $lastTargetIdx; $m++){ #iterate over all targets
        my $shortListValue = ".";
        my $currentTarget=$arrayRefs[0][$m];
        my $currentLlValue = $arrayRefs[3][$m];
        my $currentAutoValue = $arrayRefs[5][$m];
        chomp $currentLlValue;
        if ($currentAutoValue eq "NA" || $currentAutoValue > 0.10){ #if auto vc is "NA" or higher than 0.10, assign low quality to it (this target)
            $shortListValue = "LOW_QUALITY";
            $geneCounts{$arrayRefs[1][$m]}++; #Low quality detected, count it
        }
        $totalGeneCounts{$arrayRefs[1][$m]}++; #Count all occurences of genes to use in later calculation
        
        #Nothing has to be updated
        $arrayRefs[4]->[$m] = $shortListValue; #push values into shortlist array
    }
    
    # Count if half or more of targets belonging to one gene fail thresholds, if so, do gene correction
    my $highQualityTargets=0;
    for (my $g=0; $g <= $lastTargetIdx; $g++){
        my $gene = $arrayRefs[1][$g];
        my $nextGene;
        if ($g == $lastTargetIdx) { #If next element doesn't exist, create fake one
            $nextGene = "LAST";
        }else{
            $nextGene = $arrayRefs[1][$g+1];
        }
        if ($gene ne $nextGene) { #if next gene is different than current gene set high quality counter back to 0
            $highQualityTargets=0;
        }
        
        if (exists $geneCountsForGeneCorrection{ $gene }){ #Foreach target retrieve gene, if exists continue
            my $toCheck = ($geneCountsForGeneCorrection{ $gene } / $totalGeneCountsForGeneCorrection{ $gene });
            if ($toCheck >= 0.5){ #if more than 50 % of targets have abberation, do gene check
                my $currentValue=$arrayRefs[2][$g];
                my @vals=split("\t", $currentValue);
                my $autoratio = $vals[0];
                my $autozscore = $vals[1];
                my $currentQuality = $arrayRefs[4][$g];
                my $abberationValue;
                if ($currentQuality eq ".") { #current targets quality is high
                    #If LOW or HIGH thresholds passed, update the abberation value
                    if ($autoratio < $params->{ratioCutOffLow} && $autozscore < $params->{zScoreCutOffLow}){
                        $abberationValue = "DEL";
                        $arrayRefs[3]->[$g] = $abberationValue;
                        $highQualityTargets++;
                    }elsif ($autoratio > $params->{ratioCutOffHigh} && $autozscore > $params->{zScoreCutOffHigh}){
                        $abberationValue = "DUP";
                        $arrayRefs[3]->[$g] = $abberationValue;
                        $highQualityTargets++;
                    }else {
                        #No updates to be done;
                    }
                }else{ #Low quality target, check if more high quality targets are marked as DEL or DUP
                    if ($autoratio < $params->{ratioCutOffLow} || $autoratio > $params->{ratioCutOffHigh} && $highQualityTargets > 0) { #if autoratio is either LOW or HIGH and other high quality targets passed ratio and zscore
                        if ($autoratio < $params->{ratioCutOffLow}){
                            $abberationValue = "DEL";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }elsif ($autoratio > $params->{ratioCutOffHigh}){
                            $abberationValue = "DUP";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }
                    }else{#Just use normal thresholds, autoratio and autozscore
                        if ($autoratio < $params->{ratioCutOffLow} && $autozscore < $params->{zScoreCutOffLow}){
                            $abberationValue = "DEL";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }elsif ($autoratio > $params->{ratioCutOffHigh} && $autozscore > $params->{zScoreCutOffHigh}){
                            $abberationValue = "DUP";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }
                    }
                    
                }
            }
        }
    }
    
    #iterate over all targets, detect genes of which all targets have low quality
    for (my $m=0; $m <= $lastTargetIdx; $m++){
        my $currentTarget=$arrayRefs[0][$m];
        my $currentGene=$arrayRefs[1][$m];
        if (exists $geneCounts{ $currentGene }){ #in current gene at least one target with low quality detected
            if ( ($geneCounts{ $currentGene } / $totalGeneCounts{ $currentGene }) == 1) { #All targets in gene have low quality, alter the shortlist column in output
                $arrayRefs[4]->[$m] = "LOW_QUALITY,FAILED";
            }
        }
    }
    
    #iterate over all targets, detect homozygous deletions or duplications
    for (my $h=0; $h <= $lastTargetIdx; $h++){
        my $currentValue=$arrayRefs[2][$h];
        my @vals=split("\t", $currentValue); #Extract auto ratio, auto zscore and gene zscore
        my $autoratio = $vals[0];
        my $currentAbberation = $arrayRefs[3][$h];
        if ($currentAbberation eq "DEL" && $autoratio < 0.10){
            #Update abberation value to homozygous
            $arrayRefs[3]->[$h] = "HOM_DEL";
        }elsif ($currentAbberation eq "DUP" && $autoratio > 1.75){
            #Update abberation value to homozygous duplication
            $arrayRefs[3]->[$h] = "HOM_DUP";
        }else{
            # Don't update anything
        }
    }
    
    my $outputfileTotal = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.totallist.txt"; #Output total filename
    my $outputfileLong  = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.longlist.txt"; #Output long filename
    my $outputfileShort = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.shortlist.txt"; #Output short filename
    my $outputfileTotalBed = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.totallist.bed"; #Output total filename
    my $outputfileLongBed  = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.longlist.bed"; #Output long filename
    my $outputfileShortBed = $params->{outputdir}."/".$outputPostfixRemoved.".best.score.shortlist.bed"; #Output short filename
    
    
    my $lin = "CHR\tSTART\tSTOP\tGENE\tTARGET\tAUTO_RATIO\tAUTO_ZSCORE\tAUTO_VC\tGENE_RATIO\tGENE_ZSCORE\tGENE_VC\tABBERATION\tQUALITY\tSHAPIRO-WILK\n"; #Set header for output files
    my $longShortLin = "CHR\tSTART\tSTOP\tGENE\tTARGET\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION\n";
    my $outputTotalToWrite .= $lin; #concatenate full generated line to files
    my $outputLongToWrite .= $longShortLin;
    my $outputShortToWrite .= $longShortLin;
    
    my $bed_line = "#gfftags\n";
    my $outputfileLongBedToWrite   = $bed_line;
    my $outputfileShortBedToWrite  = $bed_line;
    
    my $abberationCount = 0;
    my $shapPassCount = 0;
    my $highQualAbberationCount = 0;
    my @chrStarts;
    my @chrStartsHOM;
    my $abberationCountHOM = 0;
    my $shapPassCountHOM = 0;
    my $highQualAbberationCountHOM = 0;
    for (my $m=0; $m <= $lastTargetIdx; $m++){ #Iterate over targets
        my $target = $arrayRefs[0][$m]; #Extract current target
        $target =~ s/:/\t/g;
        $target =~ s/-/\t/g;
        $target =~ s/^23/X/gs;
        $target =~ s/^24/Y/gs;
        $target =~ s/^25/MT/gs;
        my $gene = $arrayRefs[1][$m];
        my $genetargets = $arrayRefs[7][$m];
        my $vals = $arrayRefs[2][$m];
        $vals =~ s/&&/\t/g;
        my @vals_array = split /\t/, $vals;
        my $abberation = $arrayRefs[3][$m];
        my $quality = $arrayRefs[4][$m];
        my $shap = $arrayRefs[6][$m];
        my @array = split("\t", $target);
        my $chr = $array[0];
        my $start = $array[1];
        my $stop = $array[2];
        
        #Add failed regions to output total file
        $outputTotalToWrite .= join "\t", $chr,
                                          $start,
                                          $stop,
                                          $gene,
                                          $genetargets,
                                          $vals,
                                          $abberation,
                                          $quality,
                                          $shap."\n"; 
        
        if ($abberation eq "DEL" || $abberation eq "DUP" || $abberation eq "HOM_DEL" || $abberation eq "HOM_DUP") { #If abberation detected, check next abberation, afterwards write to longlist file
            my $nextAbberation = $arrayRefs[3][$m+1];
            my $nextGene = $arrayRefs[1][$m+1];
            if ($m == $lastTargetIdx){
                $nextAbberation = "LAST";
                $nextGene = "LAST";
            }
            
            if ($abberationCount == 0) {
                #$outputLongToWrite .= "$chr\t$start\t"; #Write event chr and start
            }
            push(@chrStarts, $chr);
            push(@chrStarts, $start);
            if ($quality eq ".") { #If abberation is not filtered out due to low quality autoVC or failing gene low quality, push into shortlist
                $highQualAbberationCount++; #Increase number of high quality abberation counts
            }
            #Count number of targets passing shapiro-wilk test
            if ($shap >= 0.05){
                $shapPassCount++;
            }
            
            $nextAbberation =~ s/HOM_//gs;
            $abberation =~ s/HOM_//gs;
            
            if ($abberation eq $nextAbberation && $gene eq $nextGene){ #If next abberation is DEL or DUP within same gene, same event extended
                $abberationCount++;
            }else { #Not the same event anymore, write away stop position, gene and abberation
                my $abberationCountToPrint = ($abberationCount+1);
                
                #Retrieve first start to print from previous event
                my $targetPrev = $arrayRefs[0][$m-$abberationCount];
                $targetPrev =~ s/:/\t/g;
                $targetPrev =~ s/-/\t/g;
                $targetPrev =~ s/^23/X/gs;
                $targetPrev =~ s/^24/Y/gs;
                $targetPrev =~ s/^25/MT/gs;
                #Extract chr and start position from target
                my @arrayPrev = split("\t", $targetPrev);
                my $chrPrev = $arrayPrev[0];
                my $startPrev = $arrayPrev[1];
                my $gene_target_Prev = $arrayRefs[7][$m-$abberationCount];
                my $target_gene_abberation_interval = $gene_target_Prev." - ".$genetargets;
                if ($target_gene_abberation_interval eq "- - -") {
                    $target_gene_abberation_interval = "-";
                }

                #Write event
                $outputLongToWrite .= join "\t", $chrPrev,
                                                 $startPrev,
                                                 $stop,
                                                 $gene,
                                                 $target_gene_abberation_interval,
                                                 $abberationCountToPrint,
                                                 $shapPassCount,
                                                 $abberation."\n";

                my $namestring =  "Name=".$abberation.":".$target_gene_abberation_interval.";".
                                  "Note=NUMBER_OF_TARGETS:".$abberationCountToPrint.",".
                                  "NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST:".$shapPassCount.";";
                                                             
                $outputfileLongBedToWrite .= join "\t", $chrPrev,
                                                        $startPrev,
                                                        $stop,
                                                        $namestring."\n";
                
                
                if ($highQualAbberationCount > 0) { #If total abberation counts is equal to high quality calls all target of an abberation are PASS, so the event can be written to the shortlist
                    my $chr_print = $chr;
                    my $start_print = $start;
                    if ($abberationCount != 0 ) {
                        $chr_print   = $chrStarts[0];
                        $start_print = $chrStarts[1];
                    }
                    
                    my $target_gene_abberation_interval = $gene_target_Prev." - ".$genetargets;
                    if ($target_gene_abberation_interval eq "- - -") {
                        $target_gene_abberation_interval = "-";
                    }
                    
                    #Write high quality events
                    $outputShortToWrite .= join "\t", $chr_print,
                                                      $start_print,
                                                      $stop,
                                                      $gene,
                                                      $target_gene_abberation_interval,
                                                      $abberationCountToPrint,
                                                      $shapPassCount,
                                                      $abberation."\n";
                                                      
                    my $namestring =  "Name=".$abberation.":".$target_gene_abberation_interval.";".
                                      "Note=NUMBER_OF_TARGETS:".$abberationCountToPrint.",".
                                      "NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST:".$shapPassCount.";";

                    $outputfileShortBedToWrite .= join "\t", $chr_print,
                                                             $start_print,
                                                             $stop,
                                                             $namestring."\n";
                
                }
                $abberationCount = 0; #Reset $abberationCount to 0
                $shapPassCount = 0;
                $highQualAbberationCount = 0; #Reset $highQualAbberationCount to 0
                undef(@chrStarts); undef($chr); undef($start); undef($stop);
            }
        }else{
            #Undef chr, start array and variables
            undef(@chrStarts); undef($chr); undef($start); undef($stop);
            $abberationCount = 0; #Reset $abberationCount to 0
            $shapPassCount = 0;
        }
        
        
        #Perform second iteration to also write Homozygous events to output
        $target = $arrayRefs[0][$m];
        $target =~ s/:/\t/g;
        $target =~ s/-/\t/g;
        $target =~ s/^23/X/gs;
        $target =~ s/^24/Y/gs;
        $target =~ s/^25/MT/gs;
        $gene = $arrayRefs[1][$m];
        $genetargets = $arrayRefs[7][$m];
        $vals = $arrayRefs[2][$m];
        $vals =~ s/&&/\t/g;
        $abberation = $arrayRefs[3][$m];
        $quality = $arrayRefs[4][$m];
        $shap = $arrayRefs[6][$m];
        @array = split("\t", $target);
        if ($abberation eq "HOM_DEL" || $abberation eq "HOM_DUP") { #If abberation detected, check next abberation, afterwards write to longlist file
            my $nextAbberation = $arrayRefs[3][$m+1];
            my $nextGene = $arrayRefs[1][$m+1];
            if ($m == $lastTargetIdx){
                $nextAbberation = "LAST";
                $nextGene = "LAST";
            }
            
            #Extract chr and start position from target
            my @array = split("\t", $target);
            $chr = $array[0];
            $start = $array[1];
            $stop = $array[2];
            if ($abberationCountHOM == 0) {
                #$outputLongToWrite .= "$chr\t$start\t"; #Write event chr and start
            }
            push(@chrStartsHOM, $chr);
            push(@chrStartsHOM, $start);
            if ($quality eq ".") { #If abberation is not filtered out due to low quality autoVC or failing gene low quality, push into shortlist
                $highQualAbberationCountHOM++; #Increase number of high quality abberation counts
            }
            #Count number of targets passing shapiro-wilk test
            if ($shap >= 0.05){
                $shapPassCountHOM++;
            }
            
            if ($abberation eq $nextAbberation && $gene eq $nextGene){ #If next abberation is DEL or DUP within same gene, same event extended
                $abberationCountHOM++;
            }else { #Not the same event anymore, write away stop position, gene and abberation
                my $abberationCountToPrint = ($abberationCountHOM+1);
                
                #Retrieve first start to print for event
                $target = $arrayRefs[0][$m-$abberationCountHOM];
                $target =~ s/:/\t/g;
                $target =~ s/-/\t/g;
                $target =~ s/^23/X/gs;
                $target =~ s/^24/Y/gs;
                $target =~ s/^25/MT/gs;
                #Extract chr and start position from target
                my @arrayP = split("\t", $target);
                my $chrP = $arrayP[0];
                my $startP = $arrayP[1];
                my $gene_target_Prev = $arrayRefs[7][$m-$abberationCount];
                my $target_gene_abberation_interval = $gene_target_Prev." - ".$genetargets;
                if ($target_gene_abberation_interval eq "- - -") {
                    $target_gene_abberation_interval = "-";
                }
                #$outputLongToWrite .= "$chr\t$start\t$stop\t$gene\t$abberationCountToPrint\t$shapPassCountHOM\t$abberation\n"; #Write event end and details away
                if ($highQualAbberationCountHOM > 0) { #If total abberation counts is equal to high quality calls all target of an abberation are PASS, so the event can be written to the shortlist
                    if ($abberationCountHOM == 0 ) {
                        #$outputShortToWrite .= "$chr\t$start"; #Write event chr and start
                    }else{
                        #push(@chrStarts, $chr);
                        #push(@chrStarts, $start);
                        #$outputShortToWrite .= $chrStarts[0] . "\t" . $chrStarts[1];
                    }
                    #Write high quality events
                    $outputShortToWrite .= join "\t", $chrP,
                                                      $startP,
                                                      $stop,
                                                      $gene,
                                                      $target_gene_abberation_interval,
                                                      $abberationCountToPrint,
                                                      $shapPassCountHOM,
                                                      $abberation."\n";

                    my $namestring =  "Name=".$abberation.":".$target_gene_abberation_interval.";".
                                      "Note=NUMBER_OF_TARGETS:".$abberationCountToPrint.",".
                                      "NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST:".$shapPassCountHOM.";";
                                                             
                    $outputfileShortBedToWrite .= join "\t", $chrP,
                                                             $startP,
                                                             $stop,
                                                             $namestring."\n";
                }
                $abberationCountHOM = 0; #Reset $abberationCount to 0
                $shapPassCountHOM = 0;
                $highQualAbberationCountHOM = 0; #Reset $highQualAbberationCount to 0
                undef(@chrStartsHOM); undef($chr); undef($start); undef($stop);
            }
        }else{
            #Undef chr, start array and variables
            undef(@chrStartsHOM); undef($chr); undef($start); undef($stop);
            $abberationCountHOM = 0; #Reset $abberationCount to 0
            $shapPassCountHOM = 0;
        }
    }
    writeOutput($outputfileTotal, $outputTotalToWrite); #Write output to above specified file
    writeOutput($outputfileLong, $outputLongToWrite); #Write output to above specified file
    writeOutput($outputfileShort, $outputShortToWrite); #Write output to above specified file
    
    my $nr_of_lines = ($outputfileLongBedToWrite =~ tr/\n//);
    if ($nr_of_lines > 1) {
        $outputfileLongBedToWrite =~ s/ /%20/g;
        writeOutput($outputfileLongBed, $outputfileLongBedToWrite); #Write output to above specified file
    }

    $nr_of_lines = ($outputfileShortBedToWrite =~ tr/\n//);
    if ($nr_of_lines > 1) {
        $outputfileShortBedToWrite =~ s/ /%20/g;
        writeOutput($outputfileShortBed, $outputfileShortBedToWrite); #Write output to above specified file
    }
    
    
    undef(@arrayRefs); undef($outputTotalToWrite); undef($outputLongToWrite); undef($outputShortToWrite);
    undef(%geneCounts);
    undef(%totalGeneCounts);
}

sub allTargetNormalization {
    my $file = shift;
    my $choose = shift;
    my ($samplesToSlct, $targetsToSlct) = @_;
    
    print "#######################################\n";
    print "Starting target normalization..\n";
    
    #Retrieve number of control samples
    my $N = scalar(@$samplesToSlct);
    #Put array into hash, so search if element exists
    my %targetsToSlct = map { $_ => 1 } @$targetsToSlct; #Put array into
    #Read forward control file
    print "$file.normalized.$choose.coverage.fwd.controls.txt\n";
    print "$file.normalized.$choose.coverage.rvrs.controls.txt\n";
    open(FWDCONTROLS, "$file.normalized.$choose.coverage.fwd.controls.txt") or die("Unable to open file: $!"); #Read best match file
    my @normFwdControls= <FWDCONTROLS>;
    close(FWDCONTROLS);
    my $headerFwdControls = uc($normFwdControls[0]); #Header in uppercase
    chomp $headerFwdControls;
    #Split header
    my @headerFwdControlsArray = split("\t", $headerFwdControls);
    my @headerFwdControlsArrayIdx;
    #Retrieve indices for samples to use in further analysis
    foreach my $element (@$samplesToSlct){
        my $sample = uc($element);
        my( $index )= grep { $headerFwdControlsArray[$_] eq $sample } 0..$#headerFwdControlsArray;
        #print "$sample\t$index\n";
        push (@headerFwdControlsArrayIdx, $index);
    }
    
    #Read reverse control file
    open(RVRSCONTROLS, "$file.normalized.$choose.coverage.rvrs.controls.txt") or die("Unable to open file: $!"); #Read best match file
    my @normRvrsControls= <RVRSCONTROLS>;
    close(RVRSCONTROLS);
    my $headerRvrsControls = uc($normRvrsControls[0]); #Header in uppercase
    chomp $headerRvrsControls;
    
    my @headerRvrsControlsArray = split("\t", $headerRvrsControls);
    my @headerRvrsControlsArrayIdx;
    #Retrieve indices for samples to use in further analysis
    foreach my $element(@$samplesToSlct){
        my $sample = uc($element);
        my( $index )= grep { $headerRvrsControlsArray[$_] eq $sample } 0..$#headerRvrsControlsArray;
        push (@headerRvrsControlsArrayIdx, $index);
    }

    my $lastFwdControlsIdx = $#normFwdControls;
    my $lastRvrsControlsIdx = $#normRvrsControls;
    #Check if both files contain same number of lines, if not quit with error
    if ($lastFwdControlsIdx != $lastRvrsControlsIdx) {
        die("ERROR: files containing forward and reverse reads do not have the same amount of lines!\n");
    }
    
    #Iterate through lines in fwd and rvrs files
    for (my $i=1; $i <= $lastFwdControlsIdx; $i++){
        #Check if target passes Degrees of Freedom check
        my $fwdControlLine = $normFwdControls[$i];
        my @fwdControlLineArray = split("\t", $fwdControlLine);
        my $target = $fwdControlLineArray[0];
        if(exists($targetsToSlct{$target})) { #If target exists continue further analysis
            #Only select columns from samples we want to use
            my @fwdControlLineArrayValues;
            foreach my $idx (@headerFwdControlsArrayIdx){ #For every index value obtained before, extract value
                push(@fwdControlLineArrayValues, $fwdControlLineArray[$idx]);
            }
            my $rvrsControlLine = $normRvrsControls[$i];
            my @rvrsControlLineArray = split("\t", $rvrsControlLine);
            #Only select columns from samples we want to use
            my @rvrsControlLineArrayValues;
            foreach my $idx (@headerRvrsControlsArrayIdx){ #For every index value obtained before, extract value
                push(@rvrsControlLineArrayValues, $rvrsControlLineArray[$idx]);
            }
            
            my ($mean, $sd) = calcMeanSD(\@fwdControlLineArrayValues); #Calculate mean and SD for forward reads in target
            my $fwdmean = $mean;
            my $fwdSD = $sd;
            ($mean, $sd) = calcMeanSD(\@rvrsControlLineArrayValues); #Calculate mean and SD for reverse reads in target
            my $rvrsmean = $mean;
            my $rvrsSD = $sd;
            
            #Mean fwd + rev for target i
            my $AFRNci = (( $fwdmean + $rvrsmean ) / 2);
            
            my $fwdNormSample = $fwdControlLineArray[1];
            my $rvrsNormSample = $rvrsControlLineArray[1];

            #foreach controlsample do fwdNormSample - control value and square, push into array afterwards
            my @up;
            foreach my $val (@fwdControlLineArrayValues){ #forward values
                my $minus = ($fwdNormSample - $val);
                my $squared = ($minus * $minus);
                push(@up, $squared);
            }
            foreach my $val (@rvrsControlLineArrayValues){ #reverse values
                my $minus = ($rvrsNormSample - $val);
                my $squared = ($minus * $minus);
                push(@up, $squared);
            }
            #Sum all values
            my $sum = 0;
            map { $sum += $_ } @up;
            
            #Calculate SD for fwd and rev 
            #my $SDFRNci = sqrt((((($fwdNormSample - $fwdmean) * ($fwdNormSample - $fwdmean)) + (($rvrsNormSample - $rvrsmean) * ($rvrsNormSample - $rvrsmean))) / (($N*2)-1)));
            my $SDFRNci = sqrt((( $sum ) / (($N*2)-1)));
            #Calculate mean for fwd and rev in sample
            my $AFRNsi = (( $fwdNormSample + $rvrsNormSample ) / 2);
            #Calculate SD for fwd and rev in target of sample
            my $SDFRNsi = sqrt((($fwdNormSample - $AFRNsi) * ($fwdNormSample - $AFRNsi)) + (($rvrsNormSample - $AFRNsi) * ($rvrsNormSample - $AFRNsi)));
            #Calculate combined SD for t-test for target
            my $SDFRNccsi = sqrt((($SDFRNsi * $SDFRNsi) + (( (($N * 2) -1) * $SDFRNci * $SDFRNci) / (2*$N))));
            #Calculate
            my $tsi = (($AFRNsi-$AFRNci) / ($SDFRNccsi / sqrt(((2 * $N) +2))));
            
            #print "TSI: $tsi\n";
        }
    }
    undef(%targetsToSlct);
    print STDERR "Finished target normalization\n";
    print STDERR "#######################################\n\n";
}

################################################################################
# the next 2 subs are not called in the code at the moment
# so I am commenting them out for the time being
################################################################################
#sub targetAudit {
#    my $file = shift;
#    my $choose = shift;
#    my ($samplesToSlct) = @_;
#    
#    Read forward control file
#    print STDERR "$file.normalized.$choose.coverage.fwd.controls.txt\n";
#    print STDERR "$file.normalized.$choose.coverage.rvrs.controls.txt\n";
#    open(FWDCONTROLS, "$file.normalized.$choose.coverage.fwd.controls.txt") or die("Unable to open file: $!"); #Read best match file
#    my @normFwdControls= <FWDCONTROLS>;
#    close(FWDCONTROLS);
#    my $headerFwdControls = uc($normFwdControls[0]); #Header in uppercase
#    chomp $headerFwdControls;
#    Split header
#    my @headerFwdControlsArray = split("\t", $headerFwdControls);
#    my @headerFwdControlsArrayIdx;
#    Retrieve indices for samples to use in further analysis
#    foreach my $element (@$samplesToSlct){
#        my $sample = uc($element);
#        my( $index )= grep { $headerFwdControlsArray[$_] eq $sample } 0..$#headerFwdControlsArray;
#        push (@headerFwdControlsArrayIdx, $index);
#    }
#    
#    Read reverse control file
#    open(RVRSCONTROLS, "$file.normalized.$choose.coverage.rvrs.controls.txt") or die("Unable to open file: $!"); #Read best match file
#    my @normRvrsControls= <RVRSCONTROLS>;
#    close(RVRSCONTROLS);
#    my $headerRvrsControls = uc($normRvrsControls[0]); #Header in uppercase
#    chomp $headerRvrsControls;
#    
#    my @headerRvrsControlsArray = split("\t", $headerRvrsControls);
#    my @headerRvrsControlsArrayIdx;
#    Retrieve indices for samples to use in further analysis
#    foreach my $element(@$samplesToSlct){
#        my $sample = uc($element);
#        my( $index )= grep { $headerRvrsControlsArray[$_] eq $sample } 0..$#headerRvrsControlsArray;
#        push (@headerRvrsControlsArrayIdx, $index);
#    }
#
#    my $lastFwdControlsIdx = $#normFwdControls;
#    my $lastRvrsControlsIdx = $#normRvrsControls;
#    
#    Check if both files contain same number of lines, if not quit with error
#    if ($lastFwdControlsIdx != $lastRvrsControlsIdx) {
#        die("ERROR: files containing forward and reverse reads do not have the same amount of lines!\n");
#    }
#    my %resultAuditTtest;
#    
#    Iterate through lines in fwd and rvrs files
#    for (my $i=1; $i <= $lastFwdControlsIdx; $i++){
#        Check if target passes Degrees of Freedom check
#        my $fwdControlLine = $normFwdControls[$i];
#        my @fwdControlLineArray = split("\t", $fwdControlLine);
#        my $target = $fwdControlLineArray[0];
#        Only select columns from samples we want to use
#        my @fwdControlLineArrayValues;
#        foreach my $idx (@headerFwdControlsArrayIdx){ #For every index value obtained before, extract value
#            push(@fwdControlLineArrayValues, $fwdControlLineArray[$idx]);
#        }
#        my $rvrsControlLine = $normRvrsControls[$i];
#        my @rvrsControlLineArray = split("\t", $rvrsControlLine);
#        Only select columns from samples we want to use
#        my @rvrsControlLineArrayValues;
#        foreach my $idx (@headerRvrsControlsArrayIdx){ #For every index value obtained before, extract value
#            push(@rvrsControlLineArrayValues, $rvrsControlLineArray[$idx]);
#        }
#        
#        my ($mean, $sd) = calcMeanSD(\@fwdControlLineArrayValues); #Calculate mean and SD for forward reads in target
#        my $fwdmean = $mean;
#        my $fwdSD = $sd;
#        ($mean, $sd) = calcMeanSD(\@rvrsControlLineArrayValues); #Calculate mean and SD for reverse reads in target
#        my $rvrsmean = $mean;
#        my $rvrsSD = $sd;
#        my $numSamples = scalar(@fwdControlLineArrayValues);
#        Run audit for T-test to determine which targets can be used for further downstream analysis
#        my $resultTAi = auditTtest($fwdmean, $fwdSD, $rvrsmean, $rvrsSD, $numSamples);
#        $resultAuditTtest{ $target } = $resultTAi;
#    }
#    return(%resultAuditTtest);
#}
#
#sub auditTtest {
#    my $forwardMean = shift;
#    my $forwardSD = shift;
#    my $reverseMean = shift;
#    my $reverseSD = shift;
#    my $N = shift; #number of samples
#    first part
#    my $SDFRN = sqrt(((($forwardSD * $forwardSD)/$N) + (($reverseSD * $reverseSD)/$N)));
#    
#    second part, calc TAi
#    my $TAi = (abs(($forwardMean-$reverseMean))) / $SDFRN;
#    return($TAi);
#}

sub createNormalizedCoverageFiles {
    my $inputfile = shift;
    my $extension = shift;
    my $outputExtnsn = shift;
    my $colsToExtract = shift;
    my ($controfiles) = @_;
    
    my %autodiff;
    my %sexdiff;
    my $perfectMatch = 0;

    print STDERR "\nAnalyzing sample: $inputfile..\n";
    open(INPUTFILE, $params->{inputdir}."/".$inputfile) or die("Unable to open file: $!"); #Read count file
    my @inputfile= <INPUTFILE>;
    close(INPUTFILE);
    #Set counter for number of perfect matches between sample and controls
    #Retrieve chr, start, stop, genename, region coverage, normalized autosomal coverage and normalized coverage incl. sex chrs from file by searching column indices
    my $header = uc($inputfile[0]); #Header in uppercase
    chomp $header;
    my @colNames = split(/ /, $colsToExtract);
    my @indices = getColumnIdx($header, \@colNames);
    my $chrIdxSample      = $indices[0];
    my $startIdxSample    = $indices[1];
    my $stopIdxSample     = $indices[2];
    my $geneIdxSample     = $indices[3];
    my $targetIdxSample   = (defined $indices[4] ?  $indices[4] : undef);
    my $regcovIdxSample   = $indices[5];
    my $normAutoIdxSample = $indices[6];
    my $normSexIdxSample  = $indices[7];
    my $lastFileIdxSample =$#inputfile;
    
    #Check if output match score file containing normalized_autosomal value per region for all control samples already exists
    my ($outputPostfixRemoved,$dir,$ext) = fileparse($inputfile, $extension);
   
    # my $outputPostfixRemoved = $inputfile;
    # $outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
    my $normautofile    = $params->{outputdir}."/".$outputPostfixRemoved.".".$outputExtnsn; #Output filename
    #my $normautofiletmp = $params->{outputdir}."/".$outputPostfixRemoved.$outputExtnsn.".tmp"; #Output filename
    my $normautofiletmp   = File::Temp->new( TEMPLATE => 'tempXXXXX', SUFFIX => ".".$outputExtnsn.".tmp"  );
    my $normautofiletmp_filename = $normautofiletmp->filename;

    my $normAutoToWrite= "CHR\tSTART\tSTOP\tSAMPLE\n"; #Generate outputfile header
    for (my $i=1; $i<=$lastFileIdxSample; $i++){ #Extract information from sample file
        my $lineSample=$inputfile[$i];
        chomp $lineSample;
        my @linesSample=split("\t",$lineSample);
        my $chrSample=$linesSample[$chrIdxSample];
        $chrSample =~ s/^X/23/;
        $chrSample =~ s/^Y/24/;
        $chrSample =~ s/^MT/25/;
        my $startSample   = $linesSample[$startIdxSample];
        my $stopSample    = $linesSample[$stopIdxSample];
        my $normautoSample= $linesSample[$normAutoIdxSample];
        my $normsexSample = $linesSample[$normSexIdxSample];
        my $targetsSample = (defined $targetIdxSample ? $linesSample[$targetIdxSample] : "-");
        
        #Change parameter name later, since the selection for total or autosomal only is made earlier
        if (defined $params->{sexchr}) {
            $normAutoToWrite .= "$chrSample\t$startSample\t$stopSample\t$normsexSample\n";
        } else{
           $normAutoToWrite .= "$chrSample\t$startSample\t$stopSample\t$normautoSample\n"; 
        }
    }
    writeOutput($normautofile, $normAutoToWrite); #Write output to above specified file
    

    foreach my $ctrlfile (@$controfiles){ #Open control file

        my @absDiffsAuto; #Store absolute differences per control file
        my @absDiffsSex;
        if ($inputfile ne $ctrlfile){ #If inputfile and controlfile have the same name, do skip file

            open(CONTROLFILE, $params->{controlsdir}."/".$ctrlfile) or die("Unable to open file: $!"); #Read count file
            my @controlfile= <CONTROLFILE>;
            close(CONTROLFILE);
            #Retrieve chr, start, stop, genename, region coverage, normalized autosomal coverage and normalized coverage incl. sex chrs from file by searching column indices
            my $header = uc($controlfile[0]); #Header in uppercase
            chomp $header;
            my @colNames = split(" ", $colsToExtract);
            my @indices = getColumnIdx($header, \@colNames);
            my $chrIdxControl      = $indices[0];
            my $startIdxControl    = $indices[1];
            my $stopIdxControl     = $indices[2];
            my $geneIdxControl     = $indices[3];
            my $targetIdxControl   = (defined $indices[4] ?  $indices[4] : undef);
            my $regcovIdxControl   = $indices[5];
            my $normAutoIdxControl = $indices[6];
            my $normSexIdxControl  = $indices[7];
            my $lastFileIdxControl =$#controlfile;
            
            #Open norm auto file
            open(NORMAUTOFILE, "$normautofile") or die("Unable to open file: $!"); #Read count file
            my @normautofile= <NORMAUTOFILE>;
            close(NORMAUTOFILE);
            
            #generate first part output line
            my $firstLine = $normautofile[0];
            chomp $firstLine;
            $firstLine =~ s/^23/X/gs;
            $firstLine =~ s/^24/Y/gs;
            $firstLine =~ s/^25/MT/gs;
            my $normAutoToWrite .= "$firstLine\t$ctrlfile\n"; #concatenate full generated line to files
            #Check if sample and control file contain same amount of lines
            if ($lastFileIdxSample == $lastFileIdxControl) { #Continue
                for (my $i=1; $i<=$lastFileIdxSample; $i++){ #Extract information from sample file
                    my $lineSample=$inputfile[$i];
                    chomp $lineSample;
                    my @linesSample_fields = split("\t",$lineSample);
                    my $chrSample      = $linesSample_fields[$chrIdxSample];
                    $chrSample =~ s/^X/23/;
                    $chrSample =~ s/^Y/24/;
                    $chrSample =~ s/^MT/25/;
                    my $startSample    = $linesSample_fields[$startIdxSample];
                    my $stopSample     = $linesSample_fields[$stopIdxSample];
                    my $geneSample     = $linesSample_fields[$geneIdxSample];
                    my $targetSample   = (defined $targetIdxSample ? $linesSample_fields[$targetIdxSample] : "");
                    my $regcovSample   = $linesSample_fields[$regcovIdxSample];
                    my $normautoSample = $linesSample_fields[$normAutoIdxSample];
                    my $normsexSample  = $linesSample_fields[$normSexIdxSample];
                    my $keySample      = $lineSample;
                    #Extract information from control file
                    my $lineControl=$controlfile[$i];
                    chomp $lineControl;
                    my @linesControl_fields = split("\t",$lineControl);
                    my $chrControl      = $linesControl_fields[$chrIdxControl];
                    $chrControl =~ s/^X/23/;
                    $chrControl =~ s/^Y/24/;
                    $chrControl =~ s/^MT/25/;
                    my $startControl    = $linesControl_fields[$startIdxControl];
                    my $stopControl     = $linesControl_fields[$stopIdxControl];
                    my $geneControl     = $linesControl_fields[$geneIdxControl];
                    my $targetControl   = (defined $targetIdxControl ? $linesControl_fields[$targetIdxControl] : "");
                    my $regcovControl   = $linesControl_fields[$regcovIdxControl];
                    my $normautoControl = $linesControl_fields[$normAutoIdxControl];
                    my $normsexControl  = $linesControl_fields[$normSexIdxControl];
                    my $keyControl = $lineControl;
                    if ($chrSample == $chrControl && $startSample == $startControl && $stopSample == $stopControl){ #Check if chr, start and stop match, if not throw error and skip this file from analysis
                        my $absDiffAuto = abs($normautoSample-$normautoControl); #Calculate absolute difference autosomal coverage
                        my $absDiffSex = abs($normsexSample-$normsexControl); #Calculate absolute difference all coverage
                        push(@absDiffsAuto, $absDiffAuto);
                        push(@absDiffsSex, $absDiffSex);
                    }else{
                        #Throw error and continue with next file
                        print STDERR "File ".$params->{controlsdir}."/".$ctrlfile." and sample file do not have the same chromosome, start and stop position for regions, therefore skipping this file from the analysis\n";
                        print STDERR "Error at sample line: $chrSample\t$startSample\t$stopSample\n\n";
                        last;
                    }
                    
                    #Add original line and new norm_auto value to output line
                    my $existingLine = $normautofile[$i];
                    $existingLine =~ s/^23/X/gs;
                    $existingLine =~ s/^24/Y/gs;
                    $existingLine =~ s/^25/MT/gs;
                    chomp $existingLine;
                    #Change this parameter name later, since distinction between autosomal or total is made earlier
                    if (defined $params->{sexchr}) {
                        $normAutoToWrite .= "$existingLine\t$normsexControl\n"; #concatenate full generated line to files
                    } else{
                        $normAutoToWrite .= "$existingLine\t$normautoControl\n"; #concatenate full generated line to files
                    }
                }
                #Remove 5% highest values (biggest abs diff) from calculation
                #Determine how many elements to remove from array using 5% highest values
                my $numElem = scalar(@absDiffsAuto);
                my $numRemove = int(($numElem/100)*5);
                my @sortedAuto = sort {$b <=> $a} @absDiffsAuto; #Sort array descending
                my @sortedSex = sort {$b <=> $a} @absDiffsSex; #Sort array descending
                #Splice out first 5% elements from array
                splice(@sortedAuto, 0, $numRemove);
                splice(@sortedSex, 0, $numRemove);
                #Calculate avg abs diff per target/region
                my $sumAuto = sum(@sortedAuto);
                my $sumSex = sum(@sortedSex);
                my $numAuto = scalar(@sortedAuto);
                my $numSex = scalar(@sortedSex);
                my $avgDiffAuto = ($sumAuto/$numAuto);
                my $avgDiffSex = ($sumSex/$numSex);
                print STDERR "Control: $ctrlfile\t\t\tAvg abs diff autosomal: $avgDiffAuto\tAvg abs diff incl sex: $avgDiffSex\n";
                if ($avgDiffAuto != 0 || $avgDiffSex != 0) { #If 0 then perfect match with itself, so exclude from analysis
                    $autodiff{ $ctrlfile } = $avgDiffAuto;
                    $sexdiff{ $ctrlfile } = $avgDiffSex;
                }else{ #perfect match between sample and control, throw warning and use 1 sample less in match score output
                    $perfectMatch++;
                }
            }else{
                #Different number of regions in both files, comparison can't be made
                #Throw error and continue with next file
                print STDERR "File ".$params->{controlsdir}."/".$ctrlfile." does not contain the same number of regions as the sample file, therefore skipping this file from the analysis\n";
                #Continue with next element in controls array
                next;
            }
            
            #Write all lines to file, mv tmp file to permanent
            writeOutput($normautofiletmp_filename, $normAutoToWrite); #Write output to above specified file
            `cp $normautofiletmp_filename $normautofile`;#mv tmp file to orig file
        }
    }
    return(\%autodiff, \%sexdiff, $perfectMatch);
}

#Grep column indices and push into array
sub getColumnIdx {
    my $header = shift;
    my ($colnames) = @_;
    my @headerarray = split("\t", $header);
    my @indices;
    foreach my $columnName (@$colnames){
        my( $idx )= grep { $headerarray[$_] eq $columnName } 0..$#headerarray;
        push(@indices, $idx);
    }
    return(@indices);
}

#Calculate Mean and SD from values in array
sub calcMeanSD {
    my ($values) = @_;
    my $total=0;
    my $countEle = 0;
    my @sqdiff;
    foreach my $ele (@$values){ # Sum all values
        $total = $ele + $total;
        $countEle++;
    }
    #Calculate mean
    my $mean = ($total/$countEle);
    #Calculate variance
    #Foreach element take difference from mean and square it, push into array
    foreach my $ele (@$values){
        my $diff= ($ele-$mean);
        my $square= (abs($diff)*abs($diff));
        push(@sqdiff, $square);
    }
    #Calculate variance, sum all square and calc average
    my $totalSqDiff=0;
    foreach my $ele (@sqdiff){
        $totalSqDiff = $ele + $totalSqDiff;
    }
    my $variance = ($totalSqDiff/$countEle);
    #Calculate Standard deviation, square root of variance
    my $sd = sqrt($variance);
    return($mean, $sd);
}

#Read inputfiles by extension and push into array
sub readFile {
    my $inputdir =shift;
    my $extension=shift;
    opendir(DIR, $inputdir) or die "Cannot open directory ".$inputdir;
    my @inputfiles = grep(/$extension$/,readdir(DIR));
    return @inputfiles;
}

#Open output file and write contents from string away to file
sub writeOutput {
    #$outputfile = shift;
    #$outputToWrite = shift;
    #$params->{controlsdir} = shift;
    my ($outputfile, $outputToWrite) = @_;
    open (OUTPUT, ">$outputfile") or die "Cannot open outputfile: $outputfile\n";
    print OUTPUT $outputToWrite;
    close(OUTPUT);
    #If $params->{sampleAsControl} variable is specified write output files to controlsdir too
    if (defined $params->{sampleAsControl} && defined $params->{controlsdir} && $params->{controlsdir} ne $params->{outputDir}){
        my $c_dir = $params->{controlsdir};
        `cp $outputfile $c_dir/`;
    }
}                                                             

#Create avg count files from BAM
sub countFromBam {
    my $bam = shift;
    my $ori_file_name = shift;
    #Specify header to write in outputfile
    my $outputToWrite = "CHR\tSTART\tSTOP\tGENE\tTARGET\tREGION_COV\tAVG_AUTOSOMAL_COV\tAVG_TOTAL_COV\tAVG_GENE_COV\tNORMALIZED_AUTOSOMAL\tNORMALIZED_TOTAL\tNORMALIZED_GENE\n";
    my ($file,$dir,$ext) = fileparse($ori_file_name, qr/\.[^.]*/);
    open (BED, "< ".$params->{bedfile}) or die "Cannot open file: ".$params->{bedfile}."\n";
    my @bedfile = <BED>;
    close(BED);
    
    #Check if style of regions in BED file is normal or UCSC (so incl. "chr" in chromosomename)
    getChrMatch(\@bedfile,$bam);
     
    my $counts_tmp_file   = File::Temp->new( TEMPLATE => 'tempXXXXX', SUFFIX => '.txt'  );

    print $counts_tmp_file uc(join("\t", 'chr', 'start', 'stop', 'gene', 'target', 'regioncov'."\n"));
    foreach my $line (@bedfile){
        $line =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//; #Remove Unix and Dos style line endings
        chomp $line;
        if ($line =~ m/.+\t[0-9]{1,}\t[0-9]{1,}\t[A-Za-z0-9]{1,}.+/gs){ #Check if line corresponds to chr, start, stop (Should we check if there is a genename, or do we assume this? Otherwise we could check regions and autoincrement them)
            my @array = split("\t", $line); #read line, split by tab
            my $chr=$array[0];
            my $start=$array[1];
            my $stop=$array[2];
            my $gene=$array[3];
            my $target= $array[4] || "";
            my $extractcov = join " ", ( "samtools", 
                                        "depth",
                                        "-d",$params -> {samtoolsdepthmaxcov},
                                        "-r",
                                        $chr.":".$start."-".$stop,
                                        "-a",
                                        "-q",0,
                                        "-Q",0,
                                        $bam,
       	       	       	       	       	'| awk \'BEGIN {
                                                    sum = 0
                                                }{  
       	       	       	       	       	       	       	if($3 == '.$params -> {samtoolsdepthmaxcov}. '){
       	       	       	       	       	       	       	       	print "ERROR: Depth is equal to maxcov, please set the option -samtoolsdepthmaxcov '.$params -> {samtoolsdepthmaxcov}.' to a higher value. (Stacktrace is safe to ignore)" >"/dev/stderr";
       	       	       	       	       	       	       	       	exit 1;
       	       	       	       	       	       	       	};
       	       	       	       	       	       	       	sum+=$3
       	       	       	       	       	       	} END {
       	       	       	       	       	       	      	if(NR > 0 && $3 < '.$params -> {samtoolsdepthmaxcov}.' ){
                                                                print sum/NR
                                                        }else if (NR == 0 && $3 < '.$params -> {samtoolsdepthmaxcov}.' ) {
                                                                print 0
                                                        }else{} 
       	       	       	       	       	        }\'');
            my $regioncov = join("\n",CmdRunner($extractcov));
            chomp $regioncov;
            unless (defined $regioncov) { #Check for empty variable, if true set coverage to 0
                $regioncov = 0;
            }
            $regioncov =~ s/-nan/0/gs;
            print $counts_tmp_file join "\t", $chr, $start, $stop, $gene, $target, $regioncov."\n";
        }else{
            print STDERR "Incorrect BED file format, please check your BED file before processing.\n";
        }
    }
    
    $counts_tmp_file->seek( 0, 0 );
    my @file_data= <$counts_tmp_file>;
    my $args;
    $args->{chrIdx}    = 0;
    $args->{startIdx}  = 1;
    $args->{stopIdx}   = 2;
    $args->{geneIdx}   = 3;
    $args->{targetIdx} = 4;
    $args->{regcovIdx} = 5;
    $args->{lastFileIdx} =$#file_data;
    #Calculate coverage on gene
    my ($coverage, $genehash, $genes, $covchrauto, $covchrsex) = calcGeneCov($args, @file_data);
    #Calculate coverage including sex chromomsomes
    my ($covchrautoval, $covchrall, $counts) = calcCovAutoSex($genes, $covchrauto, $covchrsex);
    
    #Foreach line in input file write away all calculated stats/values
    $outputToWrite .= writeCountFile($args, $covchrautoval, $covchrall, $genehash, $counts, $coverage, @file_data);
    
    my $outputfile = $params->{outputdir}."/".$file.".normalized.coverage.txt"; #Output filename
    print STDERR "Writing normalized coverage counts to: $outputfile\n\n\n";
    writeOutput($outputfile, $outputToWrite); #Write output to above specified file
    print STDERR "Finished processing file: $file.bam\n";
    print STDERR "Saved normalised counts in $outputfile\n";
    print STDERR "\n"; 
    return ($outputfile);

}

#Calculate coverage per region
sub calcGeneCov {
    my $args = shift;
    my @file = @_;
    my (%coverage, %genehash, @genes, @covchrauto, @covchrsex);
    
    for (my $i=1; $i<=$args->{lastFileIdx}; $i++){
        my $line=$file[$i];
        chomp $line;
        my @lines=split("\t",$line);
        my $chr   =$lines[$args->{chrIdx}];
        $chr =~ s/chr//gs; #Remove "chr" from chromosomes

        my $start =$lines[$args->{startIdx}];
        my $stop  =$lines[$args->{stopIdx}];
        my $gene  =$lines[$args->{geneIdx}];
        my $target=(defined $args->{targetIdx} ? $lines[$args->{targetIdx}] : "");
        my $regcov=$lines[$args->{regcovIdx}];
        my $key = $line;
        #Calculate coverage from regions
        push(@genes, $gene); #Save gene
        if ($chr =~ m/[0-9]{1,2}/gs){ #autosomal coverage
            push(@covchrauto, $regcov);
        }elsif ($chr =~ m/[XYxy]{1,2}/gs){ #sex chromosomes
            push(@covchrsex, $regcov);
        }else{ #throw warning
            warn "Chromosome $chr is not supported. any target on this chromossome will be skipped\n";
        }
        #push coverage and genename into corresponding hashes
        $coverage{ $line } = $regcov;
        if (exists $genehash{ $gene }) { #check if gene exist, if yes extract current coverage and add new coverage
            my $currentcov = $genehash{ $gene };
            my $newcov = ($currentcov + $regcov);
            $genehash{ $gene } = $newcov;
        }else {#Gene does not exist in hash yet, add it to hash
            $genehash{ $gene } = $regcov;
        }    
    }

    #Return results to main program
    return (\%coverage,\%genehash, \@genes, \@covchrauto, \@covchrsex);
}


#Calculate coverage for autosomal and sex chromosomes
sub calcCovAutoSex {
    my ($genes, $covchrauto, $covchrsex) = @_;
    my $covchrautosum = 0;
    my $covchrsexsum = 0;
    my $covchrautoval;
    my $covchrall;
    my %counts;
    
    foreach my $num (@$covchrauto){
        $covchrautosum = ($covchrautosum + $num); #total autosomal coverage
    }
    foreach my $num (@$covchrsex){
        $covchrsexsum = ($covchrsexsum + $num); #total sex chromosomes coverage
    }
    my $covchrallsum = ($covchrautosum + $covchrsexsum); #total coverage all chromosomes
    #count number of regions
    my $covchrautolength=scalar(@$covchrauto);
    my $covchrsexlength=scalar(@$covchrsex);
    my $covchralllength = ($covchrautolength + $covchrsexlength);
    #calculate average coverage for all chromosomes and autosomal chromosomes only
    $covchrall = ($covchrallsum/$covchralllength);
    #if sample does not have targets in sex chromossomes $covchrautolength is 0 and calculation fails
    if ($covchrautolength == 0){
        $covchrautoval = 0;
    }else{
        $covchrautoval = ($covchrautosum/$covchrautolength);
    }
    #Count occurences of genenames
    $counts{$_}++ for @$genes;
    
    return($covchrautoval, $covchrall, \%counts);
}

#Calculate coverages
sub writeCountFile {
    my $args = shift;
    my $covchrautoval = shift;
    my $covchrall = shift;
    my $geneh = shift;
    my $countsh = shift;
    my $coverageh  = shift;
    my @file = @_;
    
    my $outputToWrite;
    
    for (my $i=1; $i<=$args->{lastFileIdx}; $i++){
        my $line=$file[$i];
        chomp $line;
        my @lines=split("\t",$line);
        my $chr   =$lines[$args->{chrIdx}];
        $chr =~ s/chr//gs; #Remove "chr" from chromosomes
        my $start = $lines[$args->{startIdx}];
        my $stop  = $lines[$args->{stopIdx}];
        my $gene  = $lines[$args->{geneIdx}];
        my $target= (defined $args->{targetIdx} ? $lines[$args->{targetIdx}] : "-");
        my $regcov= $lines[$args->{regcovIdx}];
        my $key = $line;
        
        my $genename = $geneh->{ $gene }; #total coverage for gene
        my $genecount = $countsh->{ $gene }; #number of regions on gene
        my $genecov = ($genename/$genecount); #avg coverage per gene
        #Calculate normalized coverages
	my $normAuto;
	if($covchrautoval !=0){
            $normAuto = (($coverageh->{ $key })/$covchrautoval);
        }else{$normAuto=0;}
        my $normTotal;
        if($covchrall != 0){
            $normTotal = (($coverageh->{ $key })/$covchrall);
        }else{$normTotal=0;}
        my $normGene;
        if ($genecov != 0) { #Check if coverage for complete gene is not null, if it is, don't calculate
            $normGene = (($coverageh->{ $key })/$genecov);
        }else{ #It's 0, so 
            $normGene = "0";
        }
        #Write all values
        my $lin= join "\t", $chr,
                            $start,
                            $stop,
                            $gene,
                            $target,
                            $coverageh->{ $key },
                            $covchrautoval,
                            $covchrall,
                            $genecov,
                            $normAuto,
                            $normTotal,
                            $normGene."\n";
     
        
        $outputToWrite .= $lin;
    }

    return($outputToWrite);
}

#Do duplicate removal on all BAM files
sub rmDupBam {
    my $bam=shift;
    my $rmDup_file = shift;
    my $tmp_dir = shift;
    #my $filename;
    my ($filename,$dir,$ext) = fileparse($bam, qr/\.[^.]*/);
    #print "$file\t\t$dir\t\t$ext\n";
    print STDERR "Processing file: $bam\n";

    #creating additional temp files
    my $rmdup_bam   = File::Temp->new( TEMPLATE => 'tempXXXXX',DIR => $tmp_dir, SUFFIX => '.rmdup.bam'  );
    my $aligned_sam = File::Temp->new( TEMPLATE => 'tempXXXXX',DIR => $tmp_dir, SUFFIX => '.aligned.only.sam');

    #Clean this up should be something like  samtools rmdup input.bam | samtools view -Sb -h -F 0x400 >  rmdup.bam && samtools index rmdup.bam

    #Mark duplicates command
    my $rmdup = join " ",   "samtools",
                            "rmdup",
                            $params->{inputdir}."/".$bam,
                            $rmdup_bam->filename;
    #Create rmdup index command
    my $rmdupIdx = join " ",    "samtools",
                                "index",
                                $rmdup_bam->filename;
    #Remove all duplicate reads and create SAM file command
    my $sam = join " ", "samtools",
                        "view",
                        "-F",
                        "0x400",
                        "-h",
                        $rmdup_bam->filename,
                        ">",
                        $aligned_sam->filename;
    #Convert SAM to BAM file command
    my $sam2bam = join " ", "samtools",
                            "view",
                            "-hSb",
                            $aligned_sam->filename,
                            ">",
                            $rmDup_file->filename;
    #Create BAM file index command
    my $bamIdx = join " ",  "samtools",
                            "index",
                            $rmDup_file->filename;

    #Execute the above defined steps
    warn "Executed rmdup mark dups on input bam\n". join( "\n" , CmdRunner($rmdup ));
    warn "Executed rmdupIdx: mark dups bam indexing\n". join( "\n" , CmdRunner($rmdupIdx ));
    warn "Executed sam: hard remove dup reads and store to sam\n". join( "\n" , CmdRunner($sam ));
    warn "Executed sam2bam: convert sam to bam\n". join( "\n" , CmdRunner($sam2bam ));
    warn "Executed bamIdx: bam index command\n". join( "\n" , CmdRunner($bamIdx ));
    
    
    return ($rmDup_file);
}

#Do check to see if bam and targeted bed files have same CHR naming scheme
sub getChrMatch{
    my $bed = shift;
    my $bam = shift;
    my $firstLine = $bedfile[0];
    my $hasbedStyleUCSC = 0;
    my $nobedStyleUCSC = 0;;
    my %chrs_on_bed;
    my %chrs_on_bam;
    foreach my $line (@$bed){
        $line =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//; #Remove Unix and Dos style line endings
        my @fields = split /\s/, $line; # split on space because bed files can be split on empty space (this can be regulat space characters or tabs)
        if ($fields[0] =~ m/^chr.+/gs) { #Line starts with chr, so UCSC style
            $hasbedStyleUCSC = 1;
        }else{
            $nobedStyleUCSC = 1;
        }
         $chrs_on_bed{$fields[0]} = 1;
    }
    
    if ($hasbedStyleUCSC == $nobedStyleUCSC && $nobedStyleUCSC == 1){
        die "Chromosome name style in BED file has a mixture of chrs prefixed with chr and without prefix.\n".
            "Please fix your BED file chromosome naming.\n";
    }

    #Extract header from BAM file and check chromosomename style
    my $bamHasStyleUCSC = 0;
    my $bamNostyleUCSC = 0;
    my $retrieveBamHeader = "samtools view -H $bam";
    my $bamHeader = `$retrieveBamHeader`;
    my @bamHeaderLines = split("\n", $bamHeader);
    
    foreach my $line (@bamHeaderLines){
        if ($line =~ m/^\@SQ\tSN:(.+)\tLN:.+/gs) {
            my $chr = $1;
            if ($chr =~ m/chr.+/gs) {
                $bamHasStyleUCSC = 1;
            }else{
                $bamNostyleUCSC = 1;
            }
            $chrs_on_bam{$chr} = 1;
        }
    }
    
    if ($bamHasStyleUCSC == $bamNostyleUCSC){
        if ($bamNostyleUCSC == 1){
            die "Chromosome name style in BAM file has a mixture of chrs prefixed with chr and without prefix.\n".
                "Was your bam aligned to a artificial reference (mixture of chrs from different assemblies?\n".
                "please fix error before continue";
        }else{
            warn "Could not detect any chromossomes on your bam header. Check for chromossome naming will be skipped\n".
                 "Please make sure both your bed file and your bam file use the same chr naming (with or without prefix";
            return 1;
        }
    }
    
    #Check if bam chromosomenames do correspond to BED file chromosomenames
    if ($hasbedStyleUCSC ne $bamHasStyleUCSC) {
        die "Chromosome name style in BED file does not correspond to naming style in BAM file. \n".
            "This is probably caused by using UCSC naming style in one file, and other naming style in the other file. \n".
            "Please fix your BED or BAM file chromosome naming.\n";
    }
    
    foreach my $c (keys %chrs_on_bed){
        unless (exists $chrs_on_bed{$c}){
            warn "Bed file containes target(s) on chr $c but there are not reads mapped to chr $c in the input bam file\n".
                 "This will result in 0 counts for these targets on this sample. Are you sure the bed file is correct for the \n".
                 "For the panel run that generated the input bam file(s) (Y/n)?\n";
                 my $answer = <STDIN>;
                 chomp $answer;
                 unless ($answer =~ /y/i){
                    die "Please check BED - sample correspondence and try again\n";
                 }
                 
        }
    }
    
    return 1;
}

sub undeff {
    @_[0..@_-1] = ();
}

#Usage of software
sub usage {
        print STDERR <<EOF;

#########################################################################################################
   _____      _   ___      __   _____ _____ _   _  _____ 
  / ____|    | \\ | \\ \\    / /  |  __ \\_   _| \\ | |/ ____|
 | |     ___ |  \\| |\\ \\  / /_ _| |  | || | |  \\| | |  __ 
 | |    / _ \\| . ` | \\ \\/ / _` | |  | || | | . ` | | |_ |
 | |___| (_) | |\\  |  \\  / (_| | |__| || |_| |\\  | |__| |
  \\_____\\___/|_| \\_|   \\/ \\__,_|_____/_____|_| \\_|\\_____|
  
#########################################################################################################
\t\tCoNVaDING  Copyright (C) 2015  Freerk van Dijk & Lennart Johansson
\t\t\tAvailable under the GNU LGPLv3 license
#########################################################################################################
This software detects Copy Number Variants (CNVs) in sequencing data using targeted gene panels.

NOTE: This software makes use of samtools (http://www.htslib.org/)
Please make sure you have the samtools executable added to your local environment using the \$PATH variable.

For questions please e-mail: f.van.dijk02\@umcg.nl or l.johansson\@umcg.nl
#########################################################################################################

CoNVaDING software version $version
Usage: $0 <mode> <parameters>
-h\t\t\tThis manual.

-mode\t\t\tMode to run in, one of the following required:
\t\t\tPipelineFromBams :
\t\t\t\tStart with BAM files as input, to enable duplicate
\t\t\t\tremoval use the rmdup variable.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -bed, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-rmDup, -sexChr, controlSamples]
\t\t\t\t[-regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]
\t\t\t\t[-percentageLessReliableTargets]


\t\t\tPipelineFromCounts :
\t\t\t\tStart with BAM files as input, to enable duplicate
\t\t\t\tremoval use the rmdup variable.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -bed, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-rmDup, -sexChr, controlSamples]
\t\t\t\t[-regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]
\t\t\t\t[-percentageLessReliableTargets]

\t\t\taddToControls :
\t\t\t\tStart with BAM files as input, to enable duplicate
\t\t\t\tremoval use the rmdup variable.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -bed]
\t\t\t\tOverWritten:
\t\t\t\t[-useSampleAsControl} this is necessarity true on this mode
\t\t\t\t[-controlsDir] This is set to the same as -outputDir
\t\t\t\tOPTIONAL:
\t\t\t\t[-rmDup, ]

\t\t\tStartWithBam :
\t\t\t\tStart with BAM files as input, to enable duplicate
\t\t\t\tremoval use the rmdup variable.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -bed, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-rmDup, -useSampleAsControl]

\t\t\tStartWithAvgCount :
\t\t\t\tStart with Average Count files as input. This is a five column text file
\t\t\t\twith predefined column names. Please read the manual for instructions.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -bed, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-useSampleAsControl]

\t\t\tStartWithMatchScore :
\t\t\t\tStart with Normalized Coverage files as input.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-controlSamples, -sexChr]

\t\t\tStartWithBestScore :
\t\t\t\tBest score analysis using Match score files as input.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -outputDir, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh]

\t\t\tGenerateTargetQcList :
\t\t\t\tGenerate a target QC list to use as input for finallist creation.
\t\t\t\tREQUIRED:
\t\t\t\t[-outputDir, -controlsDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-controlSamples, -regionThreshold, -ratioCutOffLow, -ratioCutOffHigh, -zScoreCutOffLow, -zScoreCutOffHigh, -sampleRatioScore]

\t\t\tCreateFinalList :
\t\t\t\tCreates the final list using the target QC list for filtering.
\t\t\t\tREQUIRED:
\t\t\t\t[-inputDir, -targetQcList, -outputDir]
\t\t\t\tOPTIONAL:
\t\t\t\t[-percentageLessReliableTargets]


PARAMETERS:
-inputDir\t\tInput directory, depending on the analysis mode this contains
\t\t\tBAM, AvgCount, normalized coverage or match score files.

-bed\t\t\tInput file specifying regions to analyze in BED format.

-outputDir\t\tOutput directory to write results to.

-controlsDir\t\tDirectory containing control samples.

-targetQcList\t\tPath to file containing target QC values.

-controlSamples\t\tNumber of samples to use in Match score analysis. DEFAULT: 30

-regionThreshold\tPercentage of all control samples differing more than 3
\t\t\tstandard deviations from mean coverage of a region in the specified
\t\t\tBED file to exlude from sample CV calculation. DEFAULT: 20

-rmDup\t\t\tSwitch to enable duplicate removal when using BAM files as input.

-sexChr\t\t\tSwitch to include sex chromosomes in analysis.

-useSampleAsControl\tSwitch to use samples as control. Example: when using BAM
\t\t\tfiles to create count files and subsequentially use the
\t\t\tgenerated count files as controls.

-ratioCutOffLow\t\tLower ratio cutoff value. Region ratio values below this
\t\t\tthreshold are marked as deletion. DEFAULT: 0.65

-ratioCutOffHigh\tHigher ratio cutoff value. Region ratio values above this
\t\t\tthreshold are marked as duplication. DEFAULT: 1.4

-zScoreCutOffLow\tLower Z-score cutoff value. Regions with a Z-score below
\t\t\tthis threshold are marked as deletion. DEFAULT: -3

-zScoreCutOffHigh\tHigher Z-score cutoff value. Regions with a Z-score above
\t\t\tthis threshold are marked as duplication. DEFAULT: 3

-sampleRatioScore\tSample CV z-score cutoff value. Sample with a ratio
\t\t\tscore below this value are excluded from analysis. DEFAULT: 0.09
\t\t\tNOTE: this variable now is named "Sample CV" in all outputs!

-percentageLessReliableTargets\tTarget labelled as less reliable in percentage
\t\t\tof control samples. DEFAULT: 20

-samtoolsdepthmaxcov\tConfigure the max coverage of 'samtools depth' tool with this
\t\t\tswitch. The count is capped at this value for an interval. DEFAULT=8000

#########################################################################################################

EOF
 
}

sub CmdRunner {
        my $ret;
        my $cmd = join(" ",@_);

        warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";

	#safety for everything
        @{$ret} = `set -e -o pipefail && ($cmd )`;
        if ($? == -1) {
                die localtime( time() ). " [ERROR] failed to execute: $!\n";
        }elsif ($? & 127) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($? & 127),  ($? & 128) ? 'with' : 'without';
        }elsif ($? != 0) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($? & 127),  ($? & 128) ? 'with' : 'without';
        }else {
               	warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $? >> 8;
        }
	return @{$ret};
}
