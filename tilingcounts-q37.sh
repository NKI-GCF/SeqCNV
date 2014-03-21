#!/bin/bash

REF=$1
FAI=$REF.fai
WINDOW=$2
NCPU=2

if [ -z "$3" ]
then
	echo "Defaulting to $NCPU cpu's"
else
	NCPU=$3
	echo "Using $NCPU cpus"
fi

#checks
if [ ! -e $REF ]
then
	echo "Reference file not found"
	exit 1
fi

if [ ! -e $FAI ]
then
	echo "Reference has no fasta index file (.fai)"
	exit 1
fi


#create a temporary file file the windows
POSTMP=windows-$WINDOW.txt

#make the tiling windows
bedtools makewindows -g $FAI -w $WINDOW > $POSTMP

#calculate the gc contents in these windows
bedtools nuc -fi $REF -bed $POSTMP > gccontent-$WINDOW.txt
#count reads
ls -1 *-s.bam |parallel -j $NCPU "samtools view -q37 -bu {} | bedtools coverage -counts -abam stdin -b $POSTMP |bedtools sort > {.}-counts-$WINDOW-q37.txt"


