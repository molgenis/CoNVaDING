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

print "\n#######################################\n";
print "COMMANDLINE OPTIONS IN AFFECT:\n";
foreach my $par (@ARGV){
    print "$par\n";
}
print "#######################################\n";

######CHANGE VERSION PARAMETER IF VERSION IS UPDATED#####
my $version = "1.2.1";

##############################################################################################
##############################################################################################
##   CoNVaDING, copy number variation detecting in next-generation sequencing gene panels   ##
##   Copyright (C) 2015  Freerk van Dijk & Lennart Johansson                                ##
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
my ($help, $mode, $inputdir, $controlsdir, $outputdir, $bedfile, $rmdup, $sexchr, $sampleAsControl, $numBestMatchSamplesCmdL, $regionThreshold, $ratioCutOffLow , $ratioCutOffHigh, $zScoreCutOffLow, $zScoreCutOffHigh, $sampleRatioScore, $targetQcList, $percentageLessReliableTargets);

#### get options
GetOptions(
                "h"                     => \$help,
                "mode=s"                => \$mode, #For options, see list below (default="StartWithBam")
                "inputDir:s"            => \$inputdir, #optional
                "controlsDir:s"         => \$controlsdir, #optional
                "outputDir=s"           => \$outputdir,
                "bed:s"                 => \$bedfile, #optional
                "controlSamples:s"      => \$numBestMatchSamplesCmdL, #optional
                "regionThreshold:s"     => \$regionThreshold, #optional
                "rmDup:s"               => \$rmdup, #optional
                "sexChr:s"              => \$sexchr, #optional
                "useSampleAsControl:s"  => \$sampleAsControl, #optional
                "ratioCutOffLow:s"      => \$ratioCutOffLow, #optional
                "ratioCutOffHigh:s"     => \$ratioCutOffHigh, #optional
                "zScoreCutOffLow:s"     => \$zScoreCutOffLow, #optional
                "zScoreCutOffHigh:s"    => \$zScoreCutOffHigh, #optional
                "sampleRatioScore:s"    => \$sampleRatioScore, #optional
                "targetQcList:s"        => \$targetQcList, #optional
                "percentageLessReliableTargets:s" => \$percentageLessReliableTargets #optional
          );
usage() and exit(1) if $help;
#Obligatory args
usage() and exit(1) unless $mode;
usage() and exit(1) unless $outputdir;
#Add more parameters later

#Check input parameters
if (not defined $mode){ #If $mode is not defined, assign default value "StartWithBam"
    $mode = "StartWithBam";
}
if ($mode eq "StartWithBam" || $mode eq "StartWithAvgCount" || $mode eq "StartWithMatchScore" || $mode eq "StartWithBestScore" || $mode eq "GenerateTargetQcList" || $mode eq "CreateFinalList") {
    #continue
}else{ #Throw error
    die "Mode $mode is not supported. Please read the manual.\n";
}

#If controlsdir is specified while running best score mode, throw warning
#Set default threshold values
if ($mode eq "StartWithBestScore"){
    if (defined $controlsdir){
        print "\n##### WARNING ##### WARNING #####\n";
        print "User specified parameter \"controlsDir\" not used in best score analysis.";
        print "\n##### WARNING ##### WARNING #####\n";
    }
    if (not defined $regionThreshold) {
        $regionThreshold = 20;
    }
    if (not defined $ratioCutOffLow) {
        $ratioCutOffLow = 0.65;
    }
    if (not defined $ratioCutOffHigh) {
        $ratioCutOffHigh = 1.4;
    }
    if (not defined $zScoreCutOffLow) {
        $zScoreCutOffLow = -3;
    }
    if (not defined $zScoreCutOffHigh) {
        $zScoreCutOffHigh = 3;
    }
}

#Set default threshold values when running in GenerateTargetQcList mode
if ($mode eq "GenerateTargetQcList"){
    if (not defined $controlsdir) { #Throw error when controlsdir is not specified
        die "Directory for controlsamples (-controlsDir) is not specified, please specify to continue analysis.\n";
    }
    if (not defined $regionThreshold) {
        $regionThreshold = 20;
    }
    if (not defined $ratioCutOffLow) {
        $ratioCutOffLow = 0.65;
    }
    if (not defined $ratioCutOffHigh) {
        $ratioCutOffHigh = 1.4;
    }
    if (not defined $zScoreCutOffLow) {
        $zScoreCutOffLow = -3;
    }
    if (not defined $zScoreCutOffHigh) {
        $zScoreCutOffHigh = 3;
    }
    if (not defined $sampleRatioScore) {
        $sampleRatioScore = 0.09;
    }
    #Set inputdir to same value as controlsdir
    $inputdir = $controlsdir;
}

#Check if input directory exists, otherwise error and die
if ($mode eq "StartWithBam" || $mode eq "StartWithAvgCount" || $mode eq "StartWithMatchScore" || $mode eq "StartWithBestScore" || $mode eq "CreateFinalList"){
    if (not defined $inputdir) { #Throw error when inputdir is not specified
        die "Directory for input samples (-inputDir) is not specified, please specify to continue analysis.\n";
    }
}

#Check if output directory exists, otherwise create it
if (!-d $outputdir) { #Output directory does not exist, create it
    `mkdir -p $outputdir`;
}

my @bedfile;
#Check if BED file is specified for StartWithBam and StartWithAvgCount, else die and throw error message
if ($mode eq "StartWithBam" || $mode eq "StartWithAvgCount"){
    if (not defined $bedfile){
        die "Required BED file not specified, please specify to continue analysis.\n";
    }
}elsif ($mode eq "StartWithMatchScore" || $mode eq "StartWithBestScore"){ #Throw warning saying BED file is not used in this analysis.
    if (defined $bedfile){
        print "\n##### WARNING ##### WARNING #####\n";
        print "User specified parameter \"bed\" not used in this analysis.";
        print "\n##### WARNING ##### WARNING #####\n";
    }
}

#Check if number of control samples to use for best match is defined, else use default value
if (not defined $numBestMatchSamplesCmdL){
    $numBestMatchSamplesCmdL = 30; #Default value for number of control samples to use
}

#Check if controlsdir exists for four param options
#if ($mode eq "StartWithBam" || $mode eq "StartWithAvgCount" || $mode eq "StartWithMatchScore" || $mode eq "GenerateTargetQcList"){ #Check if controlsdir is specified
if ($mode eq "StartWithBam" || $mode eq "StartWithAvgCount" || $mode eq "StartWithMatchScore"){ #Check if controlsdir is specified
    if (not defined $controlsdir) { #Throw error when controlsdir is not specified
        die "Directory for controlsamples (-controlsDir) is not specified, please specify to continue analysis.\n";
    }else{
        #Check if controls directory exists, otherwise create it
        if (!-d $controlsdir) { #Controls directory does not exist, create it
            `mkdir -p $controlsdir`;
        }
    }
    #Check if controldir contains avg count files when useSampleAsControl paramater is not specified
    if (not defined $sampleAsControl){
        opendir(DIRHANDLE, $controlsdir) || die "Couldn't open directory $controlsdir: $!\n";
        my @txt = grep {  /\.normalized.coverage.txt$/ && -f "$controlsdir/$_" } readdir DIRHANDLE;
        closedir DIRHANDLE;
        if (@txt) { #not empty
            #continue
            my $numFiles = scalar(@txt); #If match score or best score mode is used, check if number of control sample files is greater than or equal to what user requested via parameter.
            if ($mode eq "StartWithMatchScore") {
                if ( $numFiles < $numBestMatchSamplesCmdL){
                    die "The number of controlsamples in $controlsdir is less than specified number of $numBestMatchSamplesCmdL controlsamples by \"controlSamples\" parameter.\n";
                }
            }
        }else{
            die "Directory specified by \"controlsDir\" parameter does not contain normalized coverage files.\n";
        }
    }    
}

#Checks when running in CreateFinalList mode
if ($mode eq "CreateFinalList"){ #Check if targetQcList is specified
    if (not defined $targetQcList) { #Throw error when targetQcList is not specified
        die "Target QC list (-targetQcList) is not specified, please specify to continue analysis.\n";
    }
    if (not defined $percentageLessReliableTargets) {
        $percentageLessReliableTargets = 20;
    }
}

#Global vars
my $extension;
my $rmdupfile;
my $file;
my $filename;
my $bam;
my $chr;
my $start;
my $stop;
my $gene;
my $regioncov;
my $line;
my $covchrall;
my $covchrauto;
my $key;
my $covchrautosum;
my $covchrsexsum;
my $covchrautoval;
my $outputfile;
my $outputToWrite;
my $columnName;
my $mean;
my $sd;
my $fwdbam;
my $rvrsbam;
my $outputExtnsn;
my $colsToExtract;
my $perfectMatch=0;
my $TAi;
my $autoMean;
my $autoSD;
my $autoRatio;
my $autoZscore;
my $autoVc;
my $outputdirOriginal;
#Global arrays
my @inputfiles;
my @bamstocount;
my @filestodel;
my @covchrauto;
my @covchrsex;
my @genes;
my @headerarray;
my @indices;
my @values;
my @normCountFiles;
my @fwdControlLineArray;
my @rvrsControlLineArray;
my @fwds;
my @revrs;
my @sampleRatio;
my @TNciArray;
my @passSampleRatioSamples;
#Global hashes
my %coverage;
my %genehash;
my %counts;
my %autodiff;
my %sexdiff;
my %DFfile;
my %resultAuditTtest;

#Retrieve and print starttime
my $starttime = localtime();
print "\nStarting analysis $starttime\n";

##################################################################
##################################################################
#Start analysis from BAM file
if ($mode eq "StartWithBam"){   
    #continue
    if (defined $rmdup) { #Remove duplicate switch added in cmdline
        print "\n############\nrmdup switch detected, duplicate removal included in analysis\n############\n\n";
        print "Starting removing duplicates and creating new BAM files..\n";
    }
    #Read BAM files
    print "Reading BAM files to process..\n";
    #Set file extension
    $extension = ".bam";
    readFile($inputdir, $extension);
    
    #Start analysis from BAM file
    startWithBam(\@inputfiles);

##################################################################
##################################################################
#Start analysis from average count files
}elsif ($mode eq "StartWithAvgCount"){
    #Read count TXT files
    print "Reading count files..\n";
    $extension = ".txt";
    readFile($inputdir, $extension);
    print "Starting counts analysis..\n";
    
    startWithAvgCount(\@inputfiles);

##################################################################
##################################################################
#Start analysis from match score
}elsif ($mode eq "StartWithMatchScore"){
    #continue
    if (defined $sexchr) { #Use sex chromosomes switch added in cmdline
        print "\n\n############\nsexchr switch detected, sex chromosomes included in analysis\n############\n\n";
    }
    #Read count TXT files
    print "Starting search for best match scores..\n";
    print "Reading count files..\n";
    $extension = ".normalized.coverage.txt";
    readFile($controlsdir, $extension); #Read files in controls directory
    my @controlfiles = @inputfiles;
    undeff(@inputfiles);
    readFile($inputdir, $extension); #Read files in input directory
    print "Starting match score analysis..\n";
    
    #Start analysis from match score file
    startWithMatchScore($extension, \@inputfiles, \@controlfiles);
    
##################################################################
##################################################################
#Start analysis from best score
}elsif ($mode eq "StartWithBestScore"){
    #continue
    if (defined $sexchr) { #Use sex chromosomes switch added in cmdline
        print "\n\n############\nsexchr switch detected, sex chromosomes included in analysis\n############\n\n";
    }
    #Read count TXT files
    print "Starting search for best scores..\n";
    #Read all normalized autosomal coverage control files into array
    $extension = ".normalized.autosomal.coverage.all.controls.txt";
    readFile($inputdir, $extension); #Read files in input directory
    my @normAutoControls = @inputfiles;
    undeff(@inputfiles);
    #Read best match score files
    print "Reading best match score files..\n";
    $extension = ".best.match.score.txt";
    readFile($inputdir, $extension); #Read files in input directory
    print "Starting best score analysis..\n";
    
    #Start CNV detection analysis
    startWithBestScore(\@inputfiles, \@normAutoControls);
    
    #empty output variable
    undeff($outputToWrite);
    
##################################################################
##################################################################
#Generate target QC list from all controlsamples
}elsif ($mode eq "GenerateTargetQcList"){
    #Read count TXT files
    print "Reading controls directory..\n";
    #Read all *.txt files in control directory into array
    $extension = ".normalized.coverage.txt";
    readFile($controlsdir, $extension); #Read files in controls directory
    my @controlfiles = @inputfiles;
    undeff(@inputfiles);
    readFile($inputdir, $extension); #Read files in input directory
    #Set temporary outputdir to do intermediate work
    $outputdirOriginal = $outputdir;
    $outputdir="$outputdirOriginal/tmpStartWithMatchScore/";
    `mkdir $outputdir`;
    
    #Start analysis from match score file
    startWithMatchScore($extension, \@inputfiles, \@controlfiles);
    
    #Read count TXT files
    print "Starting search for best scores..\n";
    #Read all normalized autosomal coverage control files into array
    $extension = ".normalized.autosomal.coverage.all.controls.txt";
    readFile("$outputdirOriginal/tmpStartWithMatchScore/", $extension); #Read files in input directory
    my @normAutoControls = @inputfiles;
    undeff(@inputfiles);
    #Read best match score files
    print "Reading best match score files..\n";
    $extension = ".best.match.score.txt";
    readFile("$outputdirOriginal/tmpStartWithMatchScore/", $extension); #Read files in input directory
    print "Starting best score analysis..\n";
    
    #Set temporary outputdir to do intermediate work
    $outputdir="$outputdirOriginal/tmpStartWithBestScore/";
    `mkdir $outputdir`;
    
    $inputdir="$outputdirOriginal/tmpStartWithMatchScore/";

    #Start CNV detection analysis
    startWithBestScore(\@inputfiles, \@normAutoControls);
    
    #empty output variable
    undeff($outputToWrite);
    
    #Remove temporary matchscore directory
    `rm -r $outputdirOriginal/tmpStartWithMatchScore/`;
    
    #Read all *.log files in control directory into array
    $extension = ".log";
    $inputdir = "$outputdirOriginal/tmpStartWithBestScore/";
    readFile($inputdir, $extension); #Read files in controls directory
    my @logfiles = @inputfiles;
    undeff(@inputfiles);
    
    #Extract sampleratio from *.log file

    print "\n#######################################\n";
    print "\n\nSamples failing sample CV threshold of $sampleRatioScore:\n\n";
    
    
    foreach my $logfile (@logfiles){
        my $grep = `grep SAMPLE_CV: $inputdir/$logfile`;
        chomp $grep;
        if ($grep =~ m/SAMPLE_CV: ([0-9].+)/gs){ #Check if sampleratio is below threshold(default 0.09), otherwise don't use it in targetlist analysis
            my $sampleRatio = $1;
            if ($sampleRatio <= $sampleRatioScore) {
                $logfile =~ s/.log/.totallist.txt/gs; #Change *.log extension to *.totallist.txt
                push(@passSampleRatioSamples, "$inputdir/$logfile");
            }else{
                #Print samples failing sample CV score to stdout
                print "$logfile\n";
            }
        }
    }
    
    print "\n\n#######################################\n\n";
    print "\nGenerating target QC list..\n";
    
    #Calculate target autoVC
    generateTargetQcList(\@passSampleRatioSamples);
    
    print "\nDone generating target QC list\n";

##################################################################
##################################################################
#Create final list based on target QC file
}elsif ($mode eq "CreateFinalList"){ #Apply the target filtering using the file created in previous step
    #Read shortlist TXT files
    print "Reading input directory..\n";
    #Read all *.txt files in input directory into array
    $extension = ".shortlist.txt";
    readFile($inputdir, $extension); #Read files in controls directory
    
    #Generate the final list
    createFinalList(\@inputfiles);
    
}

#Retrieve and print end time
my $endtime = localtime();
print "\nFinished analysis $endtime\n";


########################################################################################################
########## SUBS #################### SUBS #################### SUBS #################### SUBS ##########
########################################################################################################

###############################################
## Main code to generate final list using    ##
## the target QC list                        ##
sub createFinalList{
    my ($inputfiles) = @_;
    
    #Read target QC list into array;
    open(FILE, "$targetQcList") or die("Unable to open file: $!"); #Read count file
    my @file= <FILE>;
    close(FILE);
    #Retrieve chr, start, stop and genename from file by searching column indices
    my $header = uc($file[0]); #Header in uppercase
    chomp $header;
    #Extract columns from header and use as index
    my @colNames = qw(CHR START STOP GENE);
    getColumnIdx($header, \@colNames);
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
        $line=$file[$i];
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
        if ($perc > $percentageLessReliableTargets){ #More than X percentage (default 20%) of targets are low quality, so FAIL
            $targetQC{ $key } = "FAIL";
        }else{
            $targetQC{ $key } = "PASS";
        }
    }

    #Iterate over input files, asses if failing targets are within the calls in the shortlist file
    foreach my $shortlist (@$inputfiles){
        print "\n\n#######################################\n";
        print "Analyzing sample: $shortlist..\n";
        #Read file into array;
        my $dir;
        my $ext;
        ($file,$dir,$ext) = fileparse($shortlist, qr/\.[^.]*/);
        $file =~ s/.shortlist//gs;
        open(FILE, "$inputdir/$shortlist") or die("Unable to open file: $!"); #Read count file
        my @file= <FILE>;
        close(FILE);
        #Retrieve chr, start, stop, genename and region coverage from file by searching column indices
        my $header = uc($file[0]); #Header in uppercase
        chomp $header;
        $outputToWrite .= "$header\n";
        my @colNames = qw(CHR START STOP GENE);
        getColumnIdx($header, \@colNames);
        my $chrIdx = $indices[0];
        my $startIdx = $indices[1];
        my $stopIdx = $indices[2];
        my $geneIdx = $indices[3];
        my $lastFileIdx=$#file;
        #Retrieve chr, start, stop, gene and regcov for each line in avg count file
        for (my $k=1; $k<=$lastFileIdx; $k++){
            $line=$file[$k];
            chomp $line;
            my @lines=split("\t",$line);
            $chr=$lines[$chrIdx];
            $chr =~ s/X/23/g;
            $chr =~ s/Y/24/g;
            $chr =~ s/Z/25/g;
            $start=$lines[$startIdx];
            $stop=$lines[$stopIdx];
            $gene=$lines[$geneIdx];
            
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
                print "\nEvent failing target QC: $line\n";
            }else{
                $outputToWrite .= "$line\n";
            }
        }
    #Write output to *.finallist.txt file
    $outputfile = "$outputdir/$file.finallist.txt"; #Output filename
    writeOutput($outputfile, $outputToWrite); #Write output to above specified file
    print "#######################################\n\n";
    undeff($outputToWrite);
    }
}

###############################################
## Main code to generate target QC list from ##
## files in controls directory               ##
sub generateTargetQcList{
    my ($passSampleRatioSamples) = @_;
    
    #Set header for output file
    $outputToWrite = "CHR\tSTART\tSTOP\tGENE"; #Instantiate header to write
    my %autoVcHash;
    my @keyFiles;
    my $lastInputfilesIdx = $#passSampleRatioSamples; #Idx of last inputfile
    my @autoVCsToOutput;
    for (my $j=0; $j <= $lastInputfilesIdx; $j++){ #Iterate over inputfiles
        my $totallist = @$passSampleRatioSamples[$j];
        #Read file into array;
        my $dir;
        my $ext;
        ($file,$dir,$ext) = fileparse($totallist, qr/\.[^.]*/);
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
        getColumnIdx($header, \@colNames);
        my $chrIdx = $indices[0];
        my $startIdx = $indices[1];
        my $stopIdx = $indices[2];
        my $geneIdx = $indices[3];
        my $autoVcIdx = $indices[4];
        my $lastFileIdx=$#file;
        #Retrieve chr, start, stop, gene and regcov for each line in avg count file
        my @autoVCs;
        for (my $i=1; $i<=$lastFileIdx; $i++){ #Iterate over lines in file
            $line=$file[$i];
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
    
    #Write output to file
    $outputfile = "$outputdirOriginal/targetQcList.txt"; #Output filename
    writeOutput($outputfile, $outputToWrite); #Write output to above specified file

    #Remove temporary directory
    `rm -r $outputdir`;
    
}

###############################################
## Main code to extract region coverage from ##
## BAM file(s)                               ##
sub startWithBam{
    my ($inputfiles) = @_;
    foreach my $bam (@$inputfiles){
        #Check if *.bam.bai or *.bai file exist, otherwise skip this bam file
        my $dir;
        my $ext;
        ($file,$dir,$ext) = fileparse($bam, qr/\.[^.]*/);
        if (-e "$inputdir/$file.bai" || -e "$inputdir/$bam.bai") {
            #Check if duplicates need to be removed
            if (defined $rmdup){
                #Process BAM files generating duplicate removed BAM files
                rmDupBam("$bam");
                print "Starting counts analysis..\n"; #Start to count regions
                #rmdupfile is returned from rmDupBam function
                countFromBam($rmdupfile);
                print "\n";
                rmTmpBAMs($outputdir, $filename); #Remove tmp generate files while running remove duplicates
            }else{ #BAM files are already rmdupped, add them to list of files to process (retrieve them from inputdir cmdline)
                print "Starting counts analysis..\n";
                countFromBam("$inputdir/$bam");
                print "\n";
            }
            print "Writing normalized coverage counts to: $outputfile\n\n\n";
            undeff($outputToWrite);
            #If $sampleAsControl variable is specified, write coverage.txt files to controlsdir too
            #if (defined $sampleAsControl){
            #    `cp $outputfile $controlsdir`;
            #}
        }else{ #Don't process this file, because it doesn't have an indexfile
            print "##### WARNING ##### WARNING #####\nCannot find an index file for file: $inputdir/$bam, skipping this file from analysis\n##### WARNING ##### WARNING #####\n\n";
        }
    }
}

sub startWithAvgCount{
    my ($inputfiles) = @_;
    foreach my $txt (@$inputfiles){
        #Set header for output file
        $outputToWrite = "CHR\tSTART\tSTOP\tGENE\tREGION_COV\tAVG_AUTOSOMAL_COV\tAVG_TOTAL_COV\tAVG_GENE_COV\tNORMALIZED_AUTOSOMAL\tNORMALIZED_TOTAL\tNORMALIZED_GENE\n";
        #Read file into array;
        my $dir;
        my $ext;
        ($file,$dir,$ext) = fileparse($txt, qr/\.[^.]*/);
        open(FILE, "$inputdir/$txt") or die("Unable to open file: $!"); #Read count file
        my @file= <FILE>;
        close(FILE);
        #Retrieve chr, start, stop, genename and region coverage from file by searching column indices
        my $header = uc($file[0]); #Header in uppercase
        chomp $header;
        my @colNames = qw(CHR START STOP GENE REGION_COV);
        getColumnIdx($header, \@colNames);
        my $chrIdx = $indices[0];
        my $startIdx = $indices[1];
        my $stopIdx = $indices[2];
        my $geneIdx = $indices[3];
        my $regcovIdx = $indices[4];
        my $lastFileIdx=$#file;
        #Retrieve chr, start, stop, gene and regcov for each line in avg count file
        for (my $i=1; $i<=$lastFileIdx; $i++){
            $line=$file[$i];
            chomp $line;
            my @lines=split("\t",$line);
            my $chr=$lines[$chrIdx];
            my $start=$lines[$startIdx];
            my $stop=$lines[$stopIdx];
            my $gene=$lines[$geneIdx];
            my $regcov=$lines[$regcovIdx];
            $key = $line;
            #Calculate coverage from regions
            calcGeneCov($chr, $start, $stop, $gene, $regcov, $line);
        }
        #Calculate coverage including sex chromomsomes
        calcCovAutoSex(\@genes, \@covchrauto, \@covchrsex);
        
        #Foreach line in input file write away all calculated stats/values
        for (my $i=1; $i<=$lastFileIdx; $i++){
            $line=$file[$i];
            chomp $line;
            my @lines=split("\t",$line);
            my $chr=$lines[$chrIdx];
            my $start=$lines[$startIdx];
            my $stop=$lines[$stopIdx];
            my $gene=$lines[$geneIdx];
            my $regcov=$lines[$regcovIdx];
            $key = $line;
            #Write all values
            $line = "$chr\t$start\t$stop\t$gene\t$regcov"; #Produce line to print in correct order for output
            writeCountFile($line, $key, $gene, $covchrautoval, $covchrall, \%genehash, \%counts, \%coverage);
            $covchrautosum = 0;
            $covchrsexsum = 0;
        }
        $outputfile = "$outputdir/$file.normalized.coverage.txt"; #Output filename
        writeOutput($outputfile, $outputToWrite, $controlsdir); #Write output to above specified file
        #If $sampleAsControl variable is specified, write coverage.txt files to controlsdir too
        #if (defined $sampleAsControl){
        #    `cp $outputfile $controlsdir`;
        #}
        #Empty all arrays and hashes
        undef(@genes); undef(@covchrauto); undef(@covchrsex);
        undef(%counts);
        undef(%coverage);
        undef(%genehash);
        print "Finished processing file: $txt\n\n\n";
    }
}

###############################################
## Main code to select most informative con- ##
## trol samples                              ##
sub startWithMatchScore{
    my $extension = shift;
    my ($inputfiles) = $_[0];
    my ($controlfiles) = $_[1];
    foreach my $inputfile (@$inputfiles){ #Open sample file
        $outputExtnsn = "normalized.autosomal.coverage.all.controls.txt"; #Specify extension for output normalized coverage file
        $colsToExtract = "CHR START STOP GENE REGION_COV NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL"; #Specify Columns to extract
        createNormalizedCoverageFiles($inputfile, $extension, $outputExtnsn, $colsToExtract, \@$controlfiles); #Give input file to analyze, extension, columns to extract and list of controlfiles to use to function to retrieve normalized coverage
        
        #Open output bestmatch file
        my $outputPostfixRemoved = $inputfile;
        $outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
        $outputfile = "$outputdir/$outputPostfixRemoved.best.match.score.txt"; #Output filename
        #Do additional checks to detect perfect match samples and check when this occurs if there still are enough samples available in controlsdirectory to continues analysis
        my $numBestMatchSamplesToProcess = keys %autodiff; #Count number of samples to process (all - perfectMatch)
        my $numBestMatchSamples=$numBestMatchSamplesCmdL;
        if ( $perfectMatch > 0){
            if ($numBestMatchSamplesToProcess == $numBestMatchSamples){ #Detected perfect matches, but number of samples to process is the same as number of samples requested on cmdline
                print "\n##### WARNING ##### WARNING #####\n";
                print "Detected $perfectMatch perfect match(es) between sample and controlsamples, excluding these samples from the analysis.";
                print "\n##### WARNING ##### WARNING #####\n";
            }elsif ($numBestMatchSamplesToProcess >= ($numBestMatchSamples+$perfectMatch)){ #If there are more or equal #samples in controlsdir then specified on cmdline, use best #cmdline samples
                $numBestMatchSamples=$numBestMatchSamples;
            }else{ #Less samples than cmdline specified available, continue analysis but throw warning
                print "\n##### WARNING ##### WARNING #####\n";
                print "Detected $perfectMatch perfect match(es) between sample and controlsamples, which means only $numBestMatchSamplesToProcess instead of $numBestMatchSamples samples from controls directory are used for Match score analysis.";
                print "\n##### WARNING ##### WARNING #####\n";
                $numBestMatchSamples = $numBestMatchSamplesToProcess;
            }
        }
        $perfectMatch=0; #Reset perfectMatch variable
        print "\n#######################################\n";
        print "Selecting best $numBestMatchSamples control samples for analysis..\n";
        print "#######################################\n";
        my @keys;
        my @vals;
        if (defined $sexchr) { #Use sex chromosomes switch, print all abs diffs including sex chromosomes
            @keys = sort { $sexdiff{$a} <=> $sexdiff{$b} } keys(%sexdiff);
            @vals = @sexdiff{@keys};
        }else { #Only use autosomal chrs
            @keys = sort { $autodiff{$a} <=> $autodiff{$b} } keys(%autodiff);
            @vals = @autodiff{@keys};
        }
        undeff($outputToWrite);
        $outputToWrite= "SAMPLE\tSAMPLE_PATH\tCONTROL_SAMPLE\tCONTROL_SAMPLE_PATH\tAVERAGE_BEST_MATCH_SCORE\n"; #Assign header to output best match file
        
        for (my $k=0; $k < $numBestMatchSamples; $k++){
            print "Control: " . $keys[$k] . "\t\t\tAvg abs diff: " . $vals[$k] . "\n";
            my $lin = $inputfile . "\t$inputdir/$inputfile\t" . $keys[$k] . "\t$controlsdir/" . $keys[$k] . "\t" . $vals[$k] . "\n";
            $outputToWrite .= $lin; #concatenate full generated line to files
        }
        print "#######################################\n\n";
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        undef(%autodiff);
        undef(%sexdiff);
    }
}

###############################################
## Main code to detect CNVs                  ##
sub startWithBestScore{
    my ($inputfiles) = $_[0];
    my ($normAutoControls) = $_[1];
    my $lastFileIdx = $#inputfiles;
    for (my $m=0; $m<=$lastFileIdx; $m++){
        my $inputfile = $inputfiles[$m];
        my @sampleFile;
        my @controlFile;
        my @controlChr;
        my @controlStart;
        my @controlStop;
        my @controlGene;
        print "\nAnalyzing sample: $inputfile..\n";
        
        $outputToWrite= "CHR\tSTART\tSTOP\tGENE\t$inputfile\t"; #Assign header to output best match file
        open(INPUTFILE, "$inputdir/$inputfile") or die("Unable to open file: $!"); #Read best match file
        my @inputfile= <INPUTFILE>;
        close(INPUTFILE);
        my $header = uc($inputfile[0]); #Header in uppercase
        chomp $header;
        
        my @colNames = qw(SAMPLE SAMPLE_PATH CONTROL_SAMPLE CONTROL_SAMPLE_PATH);
        getColumnIdx($header, \@colNames);
        my $sampleIdx = $indices[0];
        my $samplePathIdx = $indices[1];
        my $controlIdx = $indices[2];
        my $controlPathIdx = $indices[3];
        my @array = split("\t", $header);
        my $samplename = $array[$sampleIdx];
        my $samplefile = $array[$samplePathIdx];
        my $controlname = $array[$controlIdx];
        my $controlfile = $array[$controlPathIdx];
        my $lastLine=$#inputfile;
        
        #Initialize arrays for arrayreferences
        my @NORMAUTOSOMAL;
        my @NORMTOTAL;
        my @NORMGENE;
        
        # Read samplefile into array (only do this once)
        my $currentline = $inputfile[1]; #Open first samplefile from inputfile
        my @crntlnArray = split("\t", $currentline);
        # Read input samplefile into array
        open(SAMPLEFILE, "$crntlnArray[$samplePathIdx]") or die("Unable to open file: $!"); #Read SAMPLE best match file
        @sampleFile= <SAMPLEFILE>;
        close(SAMPLEFILE);
        my $headerSample = uc($sampleFile[0]); #Header in uppercase
        chomp $headerSample;
        
        # Read normalized_autosomal values file
        #Open norm auto file, further on calculate if sample is within 3SD from mean, otherwise exclude in sampleratio calculation
        my @normAutoCon = @$normAutoControls;
        my $normautofile = $normAutoCon[$m];
        open(NORMAUTOFILE, "$inputdir/$normautofile") or die("Unable to open file: $!"); #Read count file
        my @normautofile= <NORMAUTOFILE>;
        close(NORMAUTOFILE);
        
        #Extract sample file header indices
        @colNames = qw(CHR START STOP GENE NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL NORMALIZED_GENE);
        getColumnIdx($headerSample, \@colNames);
        my $chrIdxControl = $indices[0];
        my $startIdxControl = $indices[1];
        my $stopIdxControl = $indices[2];
        my $geneIdxControl = $indices[3];
        my $normAutoIdxSample = $indices[4];
        my $normSexIdxSample = $indices[5];
        my $normGeneIdxSample = $indices[6];
        #Read through samplefile and push all three values into array of arrays;
        my $lastLineSampleFile = $#sampleFile;
        my @sampleNormAuto;
        my @sampleNormTotal;
        my @sampleNormGene;

        for (my $i=1; $i<=$lastLineSampleFile; $i++){
            my $currentline = $sampleFile[$i];
            my @crntlnArray = split("\t", $currentline);
            push(@controlChr, $crntlnArray[$chrIdxControl]);
            push(@controlStart, $crntlnArray[$startIdxControl]);
            push(@controlStop, $crntlnArray[$stopIdxControl]);
            push(@controlGene, $crntlnArray[$geneIdxControl]);
            push(@sampleNormAuto, $crntlnArray[$normAutoIdxSample]);
            push(@sampleNormTotal, $crntlnArray[$normSexIdxSample]);
            push(@sampleNormGene, $crntlnArray[$normGeneIdxSample]);
        }
        my $refNormAuto = \@sampleNormAuto; #Create array references
        my $refNormTotal = \@sampleNormTotal;
        my $refNormGene = \@sampleNormGene;
        push(@NORMAUTOSOMAL, $refNormAuto);
        push(@NORMTOTAL, $refNormTotal);
        push(@NORMGENE, $refNormGene);
        #Finalize header line for outputfile
        my $lin = "AUTO_RATIO\tAUTO_ZSCORE\tAUTO_VC\tGENE_RATIO\tGENE_ZSCORE\tGENE_VC\tSHAPIRO-WILK\n";
        $outputToWrite .= $lin; #concatenate full generated line to files
        
        #Iterate over controlfiles and put into array of arrays
        for (my $i=1; $i<=$lastLine; $i++){
            # Read line from inputfile
            my $currentline = $inputfile[$i];
            my @crntlnArray = split("\t", $currentline);
            #Read control file into array
            open(CONTROLFILE, "$crntlnArray[$controlPathIdx]") or die("Unable to open file: $!"); #Read CONTROL best match file
            @controlFile= <CONTROLFILE>;
            close(CONTROLFILE);
            #Control file is now in an array
            my $headerControl = uc($controlFile[0]); #Header in uppercase
            chomp $headerControl;
            #Extract control header columnindices
            my @colNames = qw(NORMALIZED_AUTOSOMAL NORMALIZED_TOTAL NORMALIZED_GENE);
            getColumnIdx($headerControl, \@colNames);
            my $normAutoIdxControl = $indices[0];
            my $normSexIdxControl = $indices[1];
            my $normGeneIdxControl = $indices[2];
            my $lastLineControlFile = $#controlFile;
            #Walk through contents of control file
            my @controlNormAuto;
            my @controlNormTotal;
            my @controlNormGene;
            for (my $i=1; $i<=$lastLineControlFile; $i++){
                my $currentline = $controlFile[$i];
                my @crntlnArray = split("\t", $currentline);
                push(@controlNormAuto, $crntlnArray[$normAutoIdxControl]);
                push(@controlNormTotal, $crntlnArray[$normSexIdxControl]);
                push(@controlNormGene, $crntlnArray[$normGeneIdxControl]);
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
            if (defined $sexchr) { #Use sex chromosomes switch, so NORMTOTAL needs to be used
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
            my $lin = $controlChr[$i] . ":" . $controlStart[$i] . "-" . $controlStop[$i] . "\t" . $controlGene[$i] . "\n";
            
            calcAutoRatioZscoreVc($sampleValue, $lin, \@controlAutoArray);
            
            $lin = $controlChr[$i] . "\t" . $controlStart[$i] . "\t" . $controlStop[$i] . "\t" . $controlGene[$i] . "\t$sampleValue\t";
            $outputToWrite .= $lin; #concatenate full generated line to files
            $lin = "$autoRatio\t$autoZscore\t$autoVc\t";
            $outputToWrite .= $lin; #concatenate full generated line to files
            
            my $TNsi = $NORMGENE[0]->[$i];
            calcMeanSD(\@TNciArray);
            undef(@TNciArray);
            my $TNciMean = $mean;
            my $TNciSD=  $sd;
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
            calcMeanSD(\@controlGeneArray);
            #mean and sd are returned by function, calculate ratio, z-score and variation coefficient
            #ratio, observed devided by mean
            my $geneMean = $mean;
            my $geneSD = $sd;
            my $geneRatio;
            my $geneZscore;
            my $geneVc;
            if ($geneMean == 0 || $geneSD == 0) {
                $geneRatio = "NA";
                $geneZscore = "NA";
                $geneVc = "NA";
                print "\n##### WARNING ##### WARNING #####\n";
                print $controlChr[$i] . ":" . $controlStart[$i] . "-" . $controlStop[$i] . "\t" . $controlGene[$i] . "\n";
                print "Gene normalization route not available.\nMean or Standard Deviation is 0.\nCan not calculate ratio, zscore and variation coefficient on this region\n";
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
            $lin = "$geneRatio\t$geneZscore\t$geneVc\t$pVal\n";
            $outputToWrite .= $lin; #concatenate full generated line to files
            
            undef($autoMean); undef($autoSD); undef($autoRatio); undef($autoZscore); undef($autoVc);
        }
        #Open output bestmatch file
        my $outputPostfixRemoved = $inputfile;
        $outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
        $outputfile = "$outputdir/$outputPostfixRemoved.best.score.txt"; #Output filename
        print "#######################################\n\n";
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        undeff($outputToWrite);
        
        #Write sample CV score to log file in output directory
        $outputfile = "$outputdir/$outputPostfixRemoved.best.score.log"; #Output filename
        my $lastLineIdx = $#inputfile;
        $header = uc($inputfile[0]); #Header in uppercase
        chomp $header;
        #Extract control header columnindices
        @colNames = qw(AVERAGE_BEST_MATCH_SCORE);
        getColumnIdx($header, \@colNames);
        my $avgBestMatchScoreValIdx = $indices[0];
        my @avgBestMatchScores;
        for (my $i=1; $i<=$lastLineIdx; $i++){
            my $line = $inputfile[$i];
            chomp $line;
            my @array = split("\t", $line);
            my $avgBestMatchScoreVal = $array[$avgBestMatchScoreValIdx];
            push(@avgBestMatchScores, $avgBestMatchScoreVal);
            $outputToWrite .= "$line\n"; #concatenate full generated line to files
        }
        
        ######
        my $failedRegionsToWrite = "\n\n###REGIONS FAILING USER SPECIFIED QUALITY THRESHOLD OF $regionThreshold PERCENT###\n###THESE REGIONS ARE OMMITTED FROM SAMPLE_CV CALCULATION###\n";
        my $lastLin=$#normautofile; #lastlineIdx of file
        my @idxToKeep;
        for (my $i=1; $i<=$lastLin; $i++){ #Loop through norm auto file
            my $currentline = $normautofile[$i];
            chomp $currentline;
            my @crntlnArray = split("\t", $currentline); #Complete line containing region, sample and control values
            my $region = $crntlnArray[0];
            shift(@crntlnArray); #Remove first element (in this case region) from array
            calcMeanSD(\@crntlnArray); #returns mean and sd
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
            if ($percentageFail >= $regionThreshold) {
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
            push(@calcSampleRatio, $sampleRatio[$idxKeep]); #Push sample CVs to calculate mean and sd on into new array
        }
        undef(@sampleRatio);
        #Select 95% samples (exclude low and high) for sampleRatio calculation
        my @sortedNormVal = sort { $a <=> $b } @calcSampleRatio;
        my $numValues = scalar(@sortedNormVal);
        my $low25perc = ceil((0.025*$numValues));
        my $high25perc = floor((0.975*$numValues));
        my @sliceNormVal = @sortedNormVal[($low25perc) .. ($high25perc-1)]; #Slice values out of array (why not $low25perc-1??)
        
        #Calculate mean average best match score over all control samples
        calcMeanSD(\@avgBestMatchScores);
        my $meanAvgBestMatchScore = $mean;
        #Calculate sample CV
        calcMeanSD(\@sliceNormVal); #Calculate mean and sd
        my $sampleRatio = ($sd/$mean);
        $lin = "\n\nSAMPLE_CV: $sampleRatio\nMEAN_AVERAGE_BEST_MATCHSCORE: $meanAvgBestMatchScore\n"; #Add mean average best match score and sample CV to output logfile
        $outputToWrite .= $lin;
        $outputToWrite .= $failedRegionsToWrite; #Add failed regions to output logfile
        writeOutput($outputfile, $outputToWrite); #Write output to above specified file
        print "Sample CV: $sampleRatio\n";
        print "Mean average best match score of all control samples: $meanAvgBestMatchScore\n";
        print "#######################################\n\n";
        
        createOutputLists($outputdir, $inputfile);
    }
    
}

sub calcAutoRatioZscoreVc{
    my $sampleValue = $_[0];
    my $lin = $_[1];
    my ($controlAutoArray) = $_[2];
    #Calculate mean and SD for Auto array
    calcMeanSD(\@$controlAutoArray);
    #mean and sd are returned by function, calculate ratio, z-score and variation coefficient
    #ratio, observed devided by mean
    $autoMean = $mean;
    $autoSD = $sd;
    push (@TNciArray, "$autoMean");
    if ($autoMean == 0 || $autoSD == 0) {
        $autoRatio = "NA";
        $autoZscore = "NA";
        $autoVc = "NA";
        print "\n##### WARNING ##### WARNING #####\n";
        print $lin;
        print "Mean or Standard Deviation is 0.\nCan not calculate ratio, zscore and variation coefficient on this region\n";
    }else{
        $autoRatio = ($sampleValue/$autoMean);
        #z-score, observed minus mean devided by sd
        $autoZscore = (($sampleValue-$autoMean)/$autoSD);
        #variation coefficient, sd devided by mean
        $autoVc = ($autoSD/$autoMean);
        #push autoRatio in array, to calculate the VC ratio for the complete sample
        push(@sampleRatio, $autoRatio);
    }
    return(@TNciArray, @sampleRatio, $autoRatio, $autoZscore, $autoVc);
}

sub createOutputLists{
    my $outputdir = shift;
    my $inputfile = shift;
    my $outputPostfixRemoved = $inputfile;
    $outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
    #Open best.score.txt file to extract all targets and mark them
    open(BESTSCORE, "$outputdir/$outputPostfixRemoved.best.score.txt") or die("Unable to open file: $!"); #Read best match file
    my @bestScoreFile= <BESTSCORE>;
    close(BESTSCORE);
    `rm $outputdir/$outputPostfixRemoved.best.score.txt`; #Remove *.best.score.txt file, since it is almost the same as other files produced
    my @targets;
    my @genes;
    my @values;
    my @autoVCs;
    my @shaps;
    my $header = uc($bestScoreFile[0]); #Header in uppercase
    chomp $header;
    my @colNames = qw(CHR START STOP GENE AUTO_RATIO AUTO_ZSCORE AUTO_VC GENE_RATIO GENE_ZSCORE GENE_VC SHAPIRO-WILK);
    getColumnIdx($header, \@colNames);
    my $chrIdx = $indices[0];
    my $startIdx = $indices[1];
    my $stopIdx = $indices[2];
    my $geneIdx = $indices[3];
    my $autoRatioIdx = $indices[4];
    my $autoZscoreIdx = $indices[5];
    my $autoVcIdx = $indices[6];
    my $geneRatioIdx = $indices[7];
    my $geneZscoreIdx = $indices[8];
    my $geneVcIdx = $indices[9];
    my $shapIdx = $indices[10];
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
        my $autoRatio=$line[$autoRatioIdx];
        my $autoZscore=$line[$autoZscoreIdx];
        my $autoVc=$line[$autoVcIdx];
        my $geneRatio=$line[$geneRatioIdx];
        my $geneZscore=$line[$geneZscoreIdx];
        my $geneVc=$line[$geneVcIdx];
        my $shap=$line[$shapIdx];
        my $target = "$chr:$start-$stop";
        my $value = "$autoRatio\t$autoZscore\t$autoVc\t$geneRatio\t$geneZscore\t$geneVc";
        chomp $value;
        chomp $shap;
        push(@targets, $target);
        push(@genes, $gene);
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
    my @arrayRefs = ($ref1, $ref2, $ref3, $ref4, $ref5, $ref6, $ref7);
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
                }elsif ($autoratio < $ratioCutOffLow && $autozscore < $zScoreCutOffLow && $genezscore < $zScoreCutOffLow) { #Target labelled as low
                    $abberationValue = "DEL";
                }elsif ($autoratio > $ratioCutOffHigh && $autozscore > $zScoreCutOffHigh && $genezscore > $zScoreCutOffHigh) { #Target labelled as high
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
                                        if ($autoratio < $ratioCutOffLow || $autozscore < $zScoreCutOffLow || $genezscore < $zScoreCutOffLow) {
                                            $abberationValue = "DEL";
                                        }else{
                                            $abberationValue = ".";
                                            last;
                                        }
                                    }
                                }
                                if ($currentAbberation eq "DUP") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio > $ratioCutOffHigh || $autozscore > $zScoreCutOffHigh || $genezscore > $zScoreCutOffHigh) {
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
                                        if ($autoratio < $ratioCutOffLow || $autozscore < $zScoreCutOffLow || $genezscore < $zScoreCutOffLow) {
                                            $abberationValue = "DEL";
                                        }else{
                                            $abberationValue = ".";
                                            $d=($idxs[$lastIdx]+1);
                                        }
                                    }
                                }
                                if ($currentAbberation eq "DUP") {
                                    if ($autoratio ne "NA" && $autozscore ne "NA" && $genezscore ne "NA") {
                                        if ($autoratio > $ratioCutOffHigh || $autozscore > $zScoreCutOffHigh || $genezscore > $zScoreCutOffHigh) {
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
            if ($autoratio < $ratioCutOffLow || $autoratio > $ratioCutOffHigh){ #Count if auto ratio is either LOW or HIGH
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
                    if ($autoratio < $ratioCutOffLow && $autozscore < $zScoreCutOffLow){
                        $abberationValue = "DEL";
                        $arrayRefs[3]->[$g] = $abberationValue;
                        $highQualityTargets++;
                    }elsif ($autoratio > $ratioCutOffHigh && $autozscore > $zScoreCutOffHigh){
                        $abberationValue = "DUP";
                        $arrayRefs[3]->[$g] = $abberationValue;
                        $highQualityTargets++;
                    }else {
                        #No updates to be done;
                    }
                }else{ #Low quality target, check if more high quality targets are marked as DEL or DUP
                    if ($autoratio < $ratioCutOffLow || $autoratio > $ratioCutOffHigh && $highQualityTargets > 0) { #if autoratio is either LOW or HIGH and other high quality targets passed ratio and zscore
                        if ($autoratio < $ratioCutOffLow){
                            $abberationValue = "DEL";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }elsif ($autoratio > $ratioCutOffHigh){
                            $abberationValue = "DUP";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }
                    }else{#Just use normal thresholds, autoratio and autozscore
                        if ($autoratio < $ratioCutOffLow && $autozscore < $zScoreCutOffLow){
                            $abberationValue = "DEL";
                            $arrayRefs[3]->[$g] = $abberationValue;
                        }elsif ($autoratio > $ratioCutOffHigh && $autozscore > $zScoreCutOffHigh){
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
    
    my $outputfileTotal = "$outputdir/$outputPostfixRemoved.best.score.totallist.txt"; #Output total filename
    my $outputfileLong = "$outputdir/$outputPostfixRemoved.best.score.longlist.txt"; #Output long filename
    my $outputfileShort = "$outputdir/$outputPostfixRemoved.best.score.shortlist.txt"; #Output short filename
    my $outputTotalToWrite;
    my $outputLongToWrite;
    my $outputShortToWrite;
    my $lin = "CHR\tSTART\tSTOP\tGENE\tAUTO_RATIO\tAUTO_ZSCORE\tAUTO_VC\tGENE_RATIO\tGENE_ZSCORE\tGENE_VC\tABBERATION\tQUALITY\tSHAPIRO-WILK\n"; #Set header for output files
    my $longShortLin = "CHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION\n";
    $outputTotalToWrite .= $lin; #concatenate full generated line to files
    $outputLongToWrite .= $longShortLin;
    $outputShortToWrite .= $longShortLin;
    my $abberationCount = 0;
    my $shapPassCount = 0;
    my $highQualAbberationCount = 0;
    my @chrStarts;
    my @chrStartsHOM;
    my $abberationCountHOM = 0;
    my $shapPassCountHOM = 0;
    my $highQualAbberationCountHOM = 0;
    for (my $m=0; $m <= $lastTargetIdx; $m++){ #Iterate over targets
        my $chr;
        my $start;
        my $stop;
        my $target = $arrayRefs[0][$m]; #Extract current target
        $target =~ s/:/\t/g;
        $target =~ s/-/\t/g;
        $target =~ s/^23/X/gs;
        $target =~ s/^24/Y/gs;
        $target =~ s/^25/MT/gs;
        my $gene = $arrayRefs[1][$m];
        my $vals = $arrayRefs[2][$m];
        $vals =~ s/&&/\t/g;
        my $abberation = $arrayRefs[3][$m];
        my $quality = $arrayRefs[4][$m];
        my $shap = $arrayRefs[6][$m];
        my @array = split("\t", $target);
        $outputTotalToWrite .= "$target\t$gene\t$vals\t$abberation\t$quality\t$shap\n"; #Add failed regions to output total file
        if ($abberation eq "DEL" || $abberation eq "DUP" || $abberation eq "HOM_DEL" || $abberation eq "HOM_DUP") { #If abberation detected, check next abberation, afterwards write to longlist file
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
                
                #Write the previous chr and startPos to output
                $outputLongToWrite .= "$chrPrev\t$startPrev\t$stop\t$gene\t$abberationCountToPrint\t$shapPassCount\t$abberation\n"; #Write event end and details away
                if ($highQualAbberationCount > 0) { #If total abberation counts is equal to high quality calls all target of an abberation are PASS, so the event can be written to the shortlist
                    if ($abberationCount == 0 ) {
                        $outputShortToWrite .= "$chr\t$start"; #Write event chr and start
                    }else{
                        $outputShortToWrite .= $chrStarts[0] . "\t" . $chrStarts[1];
                    }
                    $outputShortToWrite .= "\t$stop\t$gene\t$abberationCountToPrint\t$shapPassCount\t$abberation\n"; #Write high quality events
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
                #$outputLongToWrite .= "$chr\t$start\t$stop\t$gene\t$abberationCountToPrint\t$shapPassCountHOM\t$abberation\n"; #Write event end and details away
                if ($highQualAbberationCountHOM > 0) { #If total abberation counts is equal to high quality calls all target of an abberation are PASS, so the event can be written to the shortlist
                    if ($abberationCountHOM == 0 ) {
                        #$outputShortToWrite .= "$chr\t$start"; #Write event chr and start
                    }else{
                        #push(@chrStarts, $chr);
                        #push(@chrStarts, $start);
                        #$outputShortToWrite .= $chrStarts[0] . "\t" . $chrStarts[1];
                    }
                    $outputShortToWrite .= "$chrP\t$startP\t$stop\t$gene\t$abberationCountToPrint\t$shapPassCountHOM\t$abberation\n"; #Write high quality events
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
        @fwdControlLineArray = split("\t", $fwdControlLine);
        my $target = $fwdControlLineArray[0];
        if(exists($targetsToSlct{$target})) { #If target exists continue further analysis
            #Only select columns from samples we want to use
            my @fwdControlLineArrayValues;
            foreach my $idx (@headerFwdControlsArrayIdx){ #For every index value obtained before, extract value
                push(@fwdControlLineArrayValues, $fwdControlLineArray[$idx]);
            }
            my $rvrsControlLine = $normRvrsControls[$i];
            @rvrsControlLineArray = split("\t", $rvrsControlLine);
            #Only select columns from samples we want to use
            my @rvrsControlLineArrayValues;
            foreach my $idx (@headerRvrsControlsArrayIdx){ #For every index value obtained before, extract value
                push(@rvrsControlLineArrayValues, $rvrsControlLineArray[$idx]);
            }
            
            calcMeanSD(\@fwdControlLineArrayValues); #Calculate mean and SD for forward reads in target
            my $fwdmean = $mean;
            my $fwdSD = $sd;
            calcMeanSD(\@rvrsControlLineArrayValues); #Calculate mean and SD for reverse reads in target
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
    print "Finished target normalization\n";
    print "#######################################\n\n";
}

sub targetAudit {
    my $file = shift;
    my $choose = shift;
    my ($samplesToSlct) = @_;
    
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
        @fwdControlLineArray = split("\t", $fwdControlLine);
        my $target = $fwdControlLineArray[0];
        #Only select columns from samples we want to use
        my @fwdControlLineArrayValues;
        foreach my $idx (@headerFwdControlsArrayIdx){ #For every index value obtained before, extract value
            push(@fwdControlLineArrayValues, $fwdControlLineArray[$idx]);
        }
        my $rvrsControlLine = $normRvrsControls[$i];
        @rvrsControlLineArray = split("\t", $rvrsControlLine);
        #Only select columns from samples we want to use
        my @rvrsControlLineArrayValues;
        foreach my $idx (@headerRvrsControlsArrayIdx){ #For every index value obtained before, extract value
            push(@rvrsControlLineArrayValues, $rvrsControlLineArray[$idx]);
        }
        
        calcMeanSD(\@fwdControlLineArrayValues); #Calculate mean and SD for forward reads in target
        my $fwdmean = $mean;
        my $fwdSD = $sd;
        calcMeanSD(\@rvrsControlLineArrayValues); #Calculate mean and SD for reverse reads in target
        my $rvrsmean = $mean;
        my $rvrsSD = $sd;
        my $numSamples = scalar(@fwdControlLineArrayValues);
        #Run audit for T-test to determine which targets can be used for further downstream analysis
        auditTtest($fwdmean, $fwdSD, $rvrsmean, $rvrsSD, $numSamples, \@fwdControlLineArrayValues, \@rvrsControlLineArrayValues);
        my $resultTAi = $TAi;
        $resultAuditTtest{ $target } = $resultTAi;
        undef($TAi);
    }
    return(%resultAuditTtest);
}

sub auditTtest {
    my $forwardMean = shift;
    my $forwardSD = shift;
    my $reverseMean = shift;
    my $reverseSD = shift;
    my $N = shift; #number of samples
    my (@fwd) = @_;
    my (@revrs) = @_;
    #first part
    my $SDFRN = sqrt(((($forwardSD * $forwardSD)/$N) + (($reverseSD * $reverseSD)/$N)));
    
    #second part, calc TAi
    $TAi = (abs(($forwardMean-$reverseMean))) / $SDFRN;
    undef(@fwdControlLineArray); undef(@rvrsControlLineArray); undef(@fwds); undef(@revrs);
    return($TAi);
}

sub createNormalizedCoverageFiles {
    my $inputfile = shift;
    my $extension = shift;
    my $outputExtnsn = shift;
    my $colsToExtract = shift;
    my ($controfiles) = @_;
    
    print "\nAnalyzing sample: $inputfile..\n";
    open(INPUTFILE, "$inputdir/$inputfile") or die("Unable to open file: $!"); #Read count file
    my @inputfile= <INPUTFILE>;
    close(INPUTFILE);
    #Set counter for number of perfect matches between sample and controls
    #Retrieve chr, start, stop, genename, region coverage, normalized autosomal coverage and normalized coverage incl. sex chrs from file by searching column indices
    my $header = uc($inputfile[0]); #Header in uppercase
    chomp $header;
    my @colNames = split(/ /, $colsToExtract);
    getColumnIdx($header, \@colNames);
    my $chrIdxSample = $indices[0];
    my $startIdxSample = $indices[1];
    my $stopIdxSample = $indices[2];
    my $geneIdxSample = $indices[3];
    my $regcovIdxSample = $indices[4];
    my $normAutoIdxSample = $indices[5];
    my $normSexIdxSample = $indices[6];
    my $lastFileIdxSample=$#inputfile;
    
    #Check if output match score file containing normalized_autosomal value per region for all control samples already exists
    my $outputPostfixRemoved = $inputfile;
    $outputPostfixRemoved =~ s/$extension//g; #Remove old extension from inputfile
    my $normautofile = "$outputdir/$outputPostfixRemoved.$outputExtnsn"; #Output filename
    my $normautofiletmp = "$outputdir/$outputPostfixRemoved.$outputExtnsn.tmp"; #Output filename

    my $normAutoToWrite= "CHR:START-STOP\tSAMPLE\n"; #Generate outputfile header
    for (my $i=1; $i<=$lastFileIdxSample; $i++){ #Extract information from sample file
        my $lineSample=$inputfile[$i];
        chomp $lineSample;
        my @linesSample=split("\t",$lineSample);
        my $chrSample=$linesSample[$chrIdxSample];
        $chrSample =~ s/^X/23/;
        $chrSample =~ s/^Y/24/;
        $chrSample =~ s/^MT/25/;
        my $startSample=$linesSample[$startIdxSample];
        my $stopSample=$linesSample[$stopIdxSample];
        my $normautoSample=$linesSample[$normAutoIdxSample];
        my $normsexSample=$linesSample[$normSexIdxSample];
        
        #Change parameter name later, since the selection for total or autosomal only is made earlier
        if (defined $sexchr) {
            $normAutoToWrite .= "$chrSample:$startSample\-$stopSample\t$normsexSample\n";
        } else{
           $normAutoToWrite .= "$chrSample:$startSample\-$stopSample\t$normautoSample\n"; 
        }
    }
    writeOutput($normautofile, $normAutoToWrite); #Write output to above specified file
    
    foreach my $ctrlfile (@$controfiles){ #Open control file
        my @absDiffsAuto; #Store absolute differences per control file
        my @absDiffsSex;
        if ($inputfile ne $ctrlfile){ #If inputfile and controlfile have the same name, do skip file
            open(CONTROLFILE, "$controlsdir/$ctrlfile") or die("Unable to open file: $!"); #Read count file
            my @controlfile= <CONTROLFILE>;
            close(CONTROLFILE);
            #Retrieve chr, start, stop, genename, region coverage, normalized autosomal coverage and normalized coverage incl. sex chrs from file by searching column indices
            my $header = uc($controlfile[0]); #Header in uppercase
            chomp $header;
            my @colNames = split(" ", $colsToExtract);
            getColumnIdx($header, \@colNames);
            my $chrIdxControl = $indices[0];
            my $startIdxControl = $indices[1];
            my $stopIdxControl = $indices[2];
            my $geneIdxControl = $indices[3];
            my $regcovIdxControl = $indices[4];
            my $normAutoIdxControl = $indices[5];
            my $normSexIdxControl = $indices[6];
            my $lastFileIdxControl=$#controlfile;
            
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
                    my @linesSample=split("\t",$lineSample);
                    my $chrSample=$linesSample[$chrIdxSample];
                    $chrSample =~ s/^X/23/;
                    $chrSample =~ s/^Y/24/;
                    $chrSample =~ s/^MT/25/;
                    my $startSample=$linesSample[$startIdxSample];
                    my $stopSample=$linesSample[$stopIdxSample];
                    my $geneSample=$linesSample[$geneIdxSample];
                    my $regcovSample=$linesSample[$regcovIdxSample];
                    my $normautoSample=$linesSample[$normAutoIdxSample];
                    my $normsexSample=$linesSample[$normSexIdxSample];
                    my $keySample = $lineSample;
                    #Extract information from control file
                    my $lineControl=$controlfile[$i];
                    chomp $lineControl;
                    my @linesControl=split("\t",$lineControl);
                    my $chrControl=$linesControl[$chrIdxControl];
                    $chrControl =~ s/^X/23/;
                    $chrControl =~ s/^Y/24/;
                    $chrControl =~ s/^MT/25/;
                    my $startControl=$linesControl[$startIdxControl];
                    my $stopControl=$linesControl[$stopIdxControl];
                    my $geneControl=$linesControl[$geneIdxControl];
                    my $regcovControl=$linesControl[$regcovIdxControl];
                    my $normautoControl=$linesControl[$normAutoIdxControl];
                    my $normsexControl=$linesControl[$normSexIdxControl];
                    my $keyControl = $lineControl;
                    if ($chrSample == $chrControl && $startSample == $startControl && $stopSample == $stopControl){ #Check if chr, start and stop match, if not throw error and skip this file from analysis
                        my $absDiffAuto = abs($normautoSample-$normautoControl); #Calculate absolute difference autosomal coverage
                        my $absDiffSex = abs($normsexSample-$normsexControl); #Calculate absolute difference all coverage
                        push(@absDiffsAuto, $absDiffAuto);
                        push(@absDiffsSex, $absDiffSex);
                    }else{
                        #Throw error and continue with next file
                        print "File $controlsdir/$ctrlfile and sample file do not have the same chromosome, start and stop position for regions, therefore skipping this file from the analysis\n";
                        print "Error at sample line: $chrSample\t$startSample\t$stopSample\n\n";
                        last;
                    }
                    
                    #Add original line and new norm_auto value to output line
                    my $existingLine = $normautofile[$i];
                    $existingLine =~ s/^23/X/gs;
                    $existingLine =~ s/^24/Y/gs;
                    $existingLine =~ s/^25/MT/gs;
                    chomp $existingLine;
                    #Change this parameter name later, since distinction between autosomal or total is made earlier
                    if (defined $sexchr) {
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
                print "Control: $ctrlfile\t\t\tAvg abs diff autosomal: $avgDiffAuto\tAvg abs diff incl sex: $avgDiffSex\n";
                if ($avgDiffAuto != 0 || $avgDiffSex != 0) { #If 0 then perfect match with itself, so exclude from analysis
                    $autodiff{ $ctrlfile } = $avgDiffAuto;
                    $sexdiff{ $ctrlfile } = $avgDiffSex;
                }else{ #perfect match between sample and control, throw warning and use 1 sample less in match score output
                    $perfectMatch++;
                }
            }else{
                #Different number of regions in both files, comparison can't be made
                #Throw error and continue with next file
                print "File $controlsdir/$ctrlfile does not contain the same number of regions as the sample file, therefore skipping this file from the analysis\n";
                #Continue with next element in controls array
                next;
            }
            
            #Write all lines to file, mv tmp file to permanent
            writeOutput($normautofiletmp, $normAutoToWrite); #Write output to above specified file
            `mv $normautofiletmp $normautofile`;#mv tmp file to orig file
        }
    }
    return(%autodiff, %sexdiff, $perfectMatch);
}

#Grep column indices and push into array
sub getColumnIdx {
    my $header = shift;
    my ($colnames) = @_;
    my @headerarray = split("\t", $header);
    undef(@indices);
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
    $mean = ($total/$countEle);
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
    $sd = sqrt($variance);
    return($mean, $sd);
    undef(@$values);
    undef(@sqdiff);
}

#Read inputfiles by extension and push into array
sub readFile {
    my $inputdir=shift;
    my $extension=shift;
    opendir(DIR, "$inputdir") or die "Cannot open directory $inputdir";
    @inputfiles = grep(/$extension$/,readdir(DIR));
    return @inputfiles;
}

#Open output file and write contents from string away to file
sub writeOutput {
    #$outputfile = shift;
    #$outputToWrite = shift;
    #$controlsdir = shift;
    my ($outputfile, $outputToWrite, $controlsdir) = @_;
    open (OUTPUT, ">$outputfile") or die "Cannot open outputfile: $outputfile\n";
    print OUTPUT $outputToWrite;
    close(OUTPUT);
    #If $sampleAsControl variable is specified write output files to controlsdir too
    if (defined $sampleAsControl && defined $controlsdir){
        `cp $outputfile $controlsdir`;
    }
}


#Create avg count files from BAM
sub countFromBam {
    $bam = shift;
    #Specify header to write in outputfile
    $outputToWrite = "CHR\tSTART\tSTOP\tGENE\tREGION_COV\tAVG_AUTOSOMAL_COV\tAVG_TOTAL_COV\tAVG_GENE_COV\tNORMALIZED_AUTOSOMAL\tNORMALIZED_TOTAL\tNORMALIZED_GENE\n";
    my $dir;
    my $ext;
    ($file,$dir,$ext) = fileparse($bam, qr/\.[^.]*/);
    open (BED, "<$bedfile") or die "Cannot open file: $bedfile\n";
    @bedfile = <BED>;
    close(BED);
    
    #Check if style of regions in BED file is normal or UCSC (so incl. "chr" in chromosomename)
    my $firstLine = $bedfile[0];
    $firstLine =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//; #Remove Unix and Dos style line endings
    my $bedStyleUCSC = "FALSE";
    if ($firstLine =~ m/^chr.+/gs) { #Line starts with chr, so UCSC style
        $bedStyleUCSC = "TRUE";
    }
    
    #Extract header from BAM file and check chromosomename style
    my $bamStyleUCSC = "FALSE";
    my $retrieveBamHeader = "samtools view -H $bam";
    my $bamHeader = `$retrieveBamHeader`;
    my @bamHeaderLines = split("\n", $bamHeader);
    my $bamFirstChr = $bamHeaderLines[1];
    if ($bamFirstChr =~ m/^\@SQ\tSN:(.+)\tLN:.+/gs) {
        my $chr = $1;
        if ($chr =~ m/chr.+/gs) {
            $bamStyleUCSC = "TRUE";
        }
    }
    
    #Check if bam chromosomenames do correspond to BED file chromosomenames
    if ($bedStyleUCSC ne $bamStyleUCSC) {
        die "Chromosome name style in BED file does not correspond to naming style in BAM file. This is probably caused by using UCSC naming style in one file, and other naming style in the other file. Please fix your BED or BAM file chromosome naming.\n";
    }
    
    foreach $line (@bedfile){
        $line =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//; #Remove Unix and Dos style line endings
        chomp $line;
        if ($line =~ m/.+\t[0-9]{1,}\t[0-9]{1,}\t[A-Za-z0-9]{1,}.+/gs){ #Check if line corresponds to chr, start, stop (Should we check if there is a genename, or do we assume this? Otherwise we could check regions and autoincrement them)
            my @array = split("\t", $line); #read line, split by tab
            $chr=$array[0];
            $start=$array[1];
            $stop=$array[2];
            $gene=$array[3];
            my $extractcov = "samtools depth -r $chr:$start-$stop -a -q 0 -Q 0 $bam | awk \'\{sum+=\$3\} END \{print sum\/NR\}\'";
            $regioncov = `$extractcov`;
            chomp $regioncov;
            if (defined $regioncov and length $regioncov) { #Check for empty variable, if true set coverage to 0
                #
            }else{
                $regioncov = 0;
            }
            $regioncov =~ s/-nan/0/gs;
            calcGeneCov($chr, $start, $stop, $gene, $regioncov, $line); #Calculate coverage per region
        }else{
            print "Incorrect BED file format, please check your BED file before processing.\n";
        }
    }
    #Calculate coverage including sex chromosomes, by iterating over bed file again
    calcCovAutoSex(\@genes, \@covchrauto, \@covchrsex);
    foreach $line (@bedfile){ #Read input file per line (region)
        chomp $line;
        my $chr;
        my $start;
        my $stop;
        my $gene;
        $key = $line;
        my @test = split("\t", $line); #read line BED file and split by tab
        $chr=$test[0];
        $start=$test[1];
        $stop=$test[2];
        $gene=$test[3];
        $gene =~ s/[\r\n]+//g;
        $line =~ s/^23/X/gs;
        $line =~ s/^24/Y/gs;
        $line =~ s/^25/MT/gs;
        #Write all values
        writeCountFile($line, $key, $gene, $covchrautoval, $covchrall, \%genehash, \%counts, \%coverage);
        $covchrautosum = 0;
        $covchrsexsum = 0;
    }
    $outputfile = "$outputdir/$file.normalized.coverage.txt"; #Output filename
    writeOutput($outputfile, $outputToWrite, $controlsdir); #Write output to above specified file
    #Empty all arrays and hashes
    undef(@genes); undef(@covchrauto); undef(@covchrsex);
    undef(%counts);
    undef(%coverage);
    undef(%genehash);
    print "Finished processing file: $bam\n";
    return ($outputfile);
}

#Calculate coverage per region
sub calcGeneCov {
    $chr = shift;
    $chr =~ s/chr//gs; #Remove "chr" from chromosomes
    $start = shift;
    $stop = shift;
    $gene = shift;
    $regioncov = shift;
    $line = shift;
    push(@genes, $gene); #Save gene
    if ($chr =~ m/[0-9]{1,2}/gs){ #autosomal coverage
        push(@covchrauto, $regioncov);
    }elsif ($chr =~ m/[XYxy]{1,2}/gs){ #sex chromosomes
        push(@covchrsex, $regioncov);
    }else{ #throw warning
        die("Chromosome $chr is not supported\n");
    }
    
    #push coverage and genename into corresponding hashes
    $coverage{ $line } = $regioncov;
    if (exists $genehash{ $gene }) { #check if gene exist, if yes extract current coverage and add new coverage
        my $currentcov = $genehash{ $gene };
        my $newcov = ($currentcov + $regioncov);
        $genehash{ $gene } = $newcov;
    }else {#Gene does not exist in hash yet, add it to hash
        $genehash{ $gene } = $regioncov;
    }
    #Return results to main program
    return (%coverage, %genehash, @genes, @covchrauto, @covchrsex);
}

#Calculate coverage for autosomal and sex chromosomes
sub calcCovAutoSex {
    my ($genes, $covchrauto, $covchrsex) = @_;
    my $covchrautosum = 0;
    my $covchrsexsum = 0;
    foreach my $num (@$covchrauto){
        $covchrautosum = ($covchrautosum + $num); #total autosomal coverage
    }
    foreach my $num (@$covchrsex){
        $covchrsexsum = ($covchrsexsum + $num); #total sex chromosomes coverage
    }
    my $covchrallsum = ($covchrautosum + $covchrsexsum); #total coverage all chromosomes
    #count number of regions
    my $covchrautolength=scalar(@covchrauto);
    my $covchrsexlength=scalar(@covchrsex);
    my $covchralllength = ($covchrautolength + $covchrsexlength);
    #calculate average coverage for all chromosomes and autosomal chromosomes only
    $covchrall = ($covchrallsum/$covchralllength);
    $covchrautoval = ($covchrautosum/$covchrautolength);
    #Count occurences of genenames
    $counts{$_}++ for @$genes;
    
    return($covchrautoval, $covchrall, %counts, @genes);
}

#Calculate coverages
sub writeCountFile {
    $line = shift;
    $key = shift;
    $gene = shift;
    $covchrautoval = shift;
    $covchrall = shift;
    
    my($geneh,$countsh,$coverageh)=@_;
    my $genename = $geneh->{ $gene }; #total coverage for gene
    my $genecount = $countsh->{ $gene }; #number of regions on gene
    my $genecov = ($genename/$genecount); #avg coverage per gene
    #Calculate normalized coverages
    my $normAuto = (($coverageh->{ $key })/$covchrautoval);
    my $normTotal = (($coverageh->{ $key })/$covchrall);
    my $normGene;
    if ($genecov != 0) { #Check if coverage for complete gene is not null, if it is, don't calculate
        $normGene = (($coverageh->{ $key })/$genecov);
    }else{ #It's 0, so 
        $normGene = "0";
    }
    chomp $line;
    my @array = split("\t", $line);
    my $chr = $array[0];
    $chr =~ s/23/X/gs;
    $chr =~ s/24/Y/gs;
    $chr =~ s/25/MT/gs;
    my $start = $array[1];
    my $stop = $array[2];
    my $gene = $array[3];
    my $line = "$chr\t$start\t$stop\t$gene";
    my $lin="$line\t" . $coverage{ $key } . "\t" . $covchrautoval . "\t" . $covchrall . "\t" . $genecov . "\t" . $normAuto . "\t" . $normTotal . "\t" . $normGene . "\n";
    $outputToWrite .= $lin; #concatenate full generated line to files
    return($outputToWrite);
}

#Do duplicate removal on all BAM files
sub rmDupBam {
    my $bam=shift;
    #my $filename;
    my $dir;
    my $ext;
    ($filename,$dir,$ext) = fileparse($bam, qr/\.[^.]*/);
    #print "$file\t\t$dir\t\t$ext\n";
    print "Processing file: $bam\n";
    #Mark duplicates command
    my $rmdup="samtools rmdup $inputdir/$bam $outputdir/$filename.rmdup.bam";
    #Create rmdup index command
    my $rmdupIdx="samtools index $outputdir/$filename.rmdup.bam";
    #Remove all duplicate reads and create SAM file command
    my $sam="samtools view -F 0x400 -h $outputdir/$filename.rmdup.bam > $outputdir/$filename.aligned.only.sam";
    #Convert SAM to BAM file command
    my $sam2bam="samtools view -hSb $outputdir/$filename.aligned.only.sam > $outputdir/$filename.aligned.only.bam";
    #Create BAM file index command
    $rmdupfile="$outputdir/$filename.aligned.only.bam";
    my $bamIdx="samtools index $rmdupfile";

    #Execute the above defined steps
    system( $rmdup ) == 0 or die "Removing duplicates from file failed: $!";
    system( $rmdupIdx ) == 0 or die "Removing duplicates from file failed: $!";
    system( $sam ) == 0 or die "Removing duplicates from file failed: $!";
    system( $sam2bam ) == 0 or die "Removing duplicates from file failed: $!";
    system( $bamIdx ) == 0 or die "Removing duplicates from file failed: $!";    
    
    return ($rmdupfile, $filename);
}

#Remove temporary bam and sam files
sub rmTmpBAMs{
    my $outputdir = shift;
    my $file = shift;
    #Remove temporary BAM files when rmdup is defined
    `rm $outputdir/$file.rmdup.bam`;
    `rm $outputdir/$file.rmdup.bam.bai`;
    `rm $outputdir/$file.aligned.only.sam`;
    `rm $outputdir/$file.aligned.only.bam`;
    `rm $outputdir/$file.aligned.only.bam.bai`;
}

sub undeff {
    @_[0..@_-1] = ();
}

#Usage of software
sub usage {
        print <<EOF;

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

Usage: ./CoNVaDING.pl <mode> <parameters>
-h\t\t\tThis manual.
-mode\t\t\tMode to run in, one of the following required:
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
#########################################################################################################

EOF
 
}
