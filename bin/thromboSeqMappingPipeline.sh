#!/bin/bash

# bash-tools and software pipeline for thromboSeq
# Software lists all FASTQ-files in the input directory, 
# creates output folders, performs per sample Trimmomatic,
# STAR-mapping, Picard tools add-read-groups, Samtools sorting,
# and HTSeq intron-spanning read summarization. Currently only
# Single Read datasets ('R1') are analyzed.
# Authors       : Myron G. Best, Edward Post
# Email         : m.best@vumc.nl
# Summary       : thromboSeq FASTQ-file processing
# Date          : 1st of September 2018

# input arguments
DIR_INPUT=$1 # directory in which FASTQ-input files are stored
DIR_OUTPUT=$2 # output directory for the analyzed data
DIR_PROGRAMS=$3 # directory with software modules
DIR_LIBRARY=$4 # directory with library files
THREADS=$5 # number of computational cores available for analysis

# output
# folder with soft symlinks to mapped and summarized HTSeq-files
# also outputs per FASTQ file a output-folder with output files from STAR-mapping,
# and a BAM and BAI file of the summarized intron-spanning reads for visualization in IGV.

echo "Start thromboSeq mapping pipeline"
echo "Check whether programs and library files are present"
# programs; ensure Trimmomatic, STAR, Picardtools, Samtools, and HTSeq are installed in the 
# same directory, e.g. 'programs'. Also ensure the software versions on your system match
# the software versions listed here below.
# the custom script filterBamCigar.pl provided with the Supplemental Data of this manuscript
# should be present in the 'programs'-directory.
export trimJar=$DIR_PROGRAMS/Trimmomatic-0.22/trimmomatic-0.22.jar
[ -e $trimJar ] && echo "trimmomatic-0.22.jar found" || echo "trimmomatic-0.22.jar not found"
export SCRIPT_STAR=$DIR_PROGRAMS/STAR_2.3.0e.Linux_x86_64/STAR
[ -e $SCRIPT_STAR ] && echo "STAR_2.3.0e.Linux_x86_64 found" || echo "STAR_2.3.0e.Linux_x86_64 not found"
export SCRIPT_JAVA=/usr/bin/java
[ -e $SCRIPT_JAVA ] && echo "java found" || echo "java not found"
export SCRIPT_PT_ADDREADGROUP=$DIR_PROGRAMS/picard-tools-1.115/AddOrReplaceReadGroups.jar
[ -e $SCRIPT_PT_ADDREADGROUP ] && echo "AddOrReplaceGroups.jar found" || echo "AddOrReplaceGroups.jar not found"
export SCRIPT_HTSEQ=$DIR_PROGRAMS/HTSeq-0.6.1/scripts/htseq-count
[ -e $SCRIPT_HTSEQ ] && echo "htseq-count found" || echo "htseq-count not found"
export filterBamCigar=$DIR_PROGRAMS/filterBamCigar.pl
[ -e $filterBamCigar ] && echo "filterBamCigar.pl found" || echo "filterBamCigar.pl not found"

# Reference files; ensure the reference files as can be downloaded from GSE107868 are 
# available in a reference file library-directory. 
# the adapter reference file TruSeq3-SE.fa provided with the Supplemental Data of this manuscript
# should be present in the 'library'-directory.
export adapters=$DIR_LIBRARY/TruSeq3-SE.fa
[ -e $adapters ] && echo "TruSeq3-SE.fa found" || echo "TruSeq3-SE.fa not found"
export INDEX_STAR_HG19=$DIR_LIBRARY/STAR-hg19-ensembl75-sjdboverhang99
[ -e $INDEX_STAR_HG19 ] && echo "STAR-hg19-ensembl75-sjdboverhang99 found" || echo "STAR-hg19-ensembl75-sjdboverhang99 not found"
export GTF=$DIR_LIBRARY/hg19/Homo_sapiens.GRCh37.75.chr1-22XY.gtf
[ -e $GTF ] && echo "Homo_sapiens.GRCh37.75.chr1-22XY.gtf found" || echo "Homo_sapiens.GRCh37.75.chr1-22XY.gtf not found"


# create output directory if not exists
if [ ! -d $DIR_OUTPUT ]
then
	mkdir $DIR_OUTPUT
fi


# list (un)zipped samples for processing
echo "" # print empty row
echo "Loop all files in input directory"

# list samples
FILES=$( find $DIR_INPUT -name '*_R1*fastq*' )
echo "Files found for processing:"
echo $FILES

# loop all files for processing
for SAMPLE in $FILES

do

	# select the sample ID
	NAME_SAMPLE=$(basename "$SAMPLE")
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	NAME_SAMPLE=$(echo $NAME_SAMPLE | cut -d _ -f 1)
	echo "" # print empty row
	echo $NAME_SAMPLE
	
# prepare output folder per sample
mkdir $DIR_OUTPUT/$NAME_SAMPLE
export SAMPLE_OUTPUT=$DIR_OUTPUT/$NAME_SAMPLE

echo "" # print empty row
# Run Trimmomatic, settings are available in the default-file
java -classpath $trimJar \
		org.usadellab.trimmomatic.TrimmomaticSE \
		-threads $THREADS \
		-phred33 \
		$SAMPLE \
		$SAMPLE_OUTPUT/$NAME_SAMPLE\_R1_clean.fastq \
		ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "" # print empty row
echo "Trimmomatic " $NAME_SAMPLE " finished"
echo "" # print empty row

# Run STAR splice-aware alignment
INPUT=$NAME_SAMPLE
INPUT_R1=$SAMPLE_OUTPUT

# Print fastq-file that will be processed
echo "Fastq-file used for $NAME_SAMPLE: $STAR_INPUT"
echo "" # print empty row

# Run STAR alignment. 
$SCRIPT_STAR \
	--runMode alignReads \
	--genomeDir $INDEX_STAR_HG19 \
	--readFilesIn $SAMPLE_OUTPUT/$NAME_SAMPLE\_R1_clean.fastq \
	--runThreadN $THREADS \
	--sjdbScore 2 \
	--outFilterMismatchNmax 10 \
	--chimSegmentMin 1 \
	--outFilterIntronMotifs 'None' \
	--outFileNamePrefix $SAMPLE_OUTPUT/$NAME_SAMPLE.

# remove trimmomatic output file
rm $SAMPLE_OUTPUT/$NAME_SAMPLE\_R1_clean.fastq
echo "STAR  " $NAME_SAMPLE " finished"
echo "" # print empty row

# Run PICARD TOOLS ADDREADGROUP
sleep 10s # 10 seconds sleep to close STAR and its virtual memory use 

echo "" # print empty row
$SCRIPT_JAVA -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_ADDREADGROUP \
	INPUT=$SAMPLE_OUTPUT/$NAME_SAMPLE.Aligned.out.sam \
	OUTPUT=$SAMPLE_OUTPUT/$NAME_SAMPLE.rg.sam \
	RGLB=RNAseq \
	RGPL=illumina \
	RGPU=platform1 \
	RGSM=$NAME_SAMPLE \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$SAMPLE_OUTPUT \
	CREATE_INDEX=true

rm $SAMPLE_OUTPUT/$NAME_SAMPLE.Aligned.out.sam # remove temporary files

echo "" # print empty row
echo "Picard Tools " $NAME_SAMPLE " finished"
echo "" # print empty row

# Run SAMTOOLS data sorting
set -e 

# Convert SAM to BAM
samtools view -@$THREADS -bS $SAMPLE_OUTPUT/$NAME_SAMPLE.rg.sam > $SAMPLE_OUTPUT/$NAME_SAMPLE.convert.bam
echo "" # print empty row
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.rg.sam # remove temporary files
echo "" # print empty row

# Sort BAM file by read name
samtools sort $SAMPLE_OUTPUT/$NAME_SAMPLE.convert.bam $SAMPLE_OUTPUT/$NAME_SAMPLE.sort1 -n -@$THREADS
echo "" # print empty row
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.convert.bam # remove temporary files
echo "" # print empty row

# Sort BAM file by chromosome location
samtools sort $SAMPLE_OUTPUT/$NAME_SAMPLE.sort1.bam $SAMPLE_OUTPUT/$NAME_SAMPLE.sorted -@$THREADS
samtools index $SAMPLE_OUTPUT/$NAME_SAMPLE.sorted.bam
echo "" # print empty row
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.sort1.bam # remove temporary files
echo "" # print empty row

samtools view -@$THREADS -h $SAMPLE_OUTPUT/$NAME_SAMPLE.sorted.bam > $SAMPLE_OUTPUT/$NAME_SAMPLE.sam
echo "" # print empty row
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.sorted.bam # remove temporary files
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.sorted.bam.bai # remove temporary files
echo "" # print empty row
echo "Samtools " $NAME_SAMPLE " finished"

# Run filterBamCigar script for selection of intron-spanning reads only
# Subsequently quantify intron-spanning read counts using HTSeq
echo "" # print empty row
$filterBamCigar $SAMPLE_OUTPUT/$NAME_SAMPLE.sam > $SAMPLE_OUTPUT/$NAME_SAMPLE.splice.sam
echo "" # print empty row
$SCRIPT_HTSEQ \
	--mode=union \
	--stranded=no \
	--minaqual=35 \
	--type=exon \
	--idattr=gene_id \
	$SAMPLE_OUTPUT/$NAME_SAMPLE.splice.sam \
	$GTF > $SAMPLE_OUTPUT/$NAME_SAMPLE.htseq.ssv

rm $SAMPLE_OUTPUT/$NAME_SAMPLE.sam # remove temporary files

# convert SAM to BAM files and generate index-file, in order to visualize mapped reads in IGV
samtools view -Sb $SAMPLE_OUTPUT/$NAME_SAMPLE.splice.sam > $SAMPLE_OUTPUT/$NAME_SAMPLE.bam
samtools index $SAMPLE_OUTPUT/$NAME_SAMPLE.bam
rm $SAMPLE_OUTPUT/$NAME_SAMPLE.splice.sam # remove temporary files

echo "" # print empty row
echo "HTSeq " $NAME_SAMPLE " finished"
echo "" # print empty row

done

# summarize all HTSeq files in an output folder
# name: HTSeq_output
if [ ! -d 'HTSeq_output' ]
then
	mkdir 'HTSeq_output'
fi

find $DIR_OUTPUT -name '*htseq.ssv' -exec ln -s {} 'HTSeq_output' \;

echo "All samples processed"

