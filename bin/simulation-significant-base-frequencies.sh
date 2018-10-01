#!/bin/bash

case $# in
    12)
        HOMOGENEOUS_CUTOFF=$1
        MIN_READS=$2
        GENOME_2_MUTATION_RATE=$3
        GENOME_LENGTH=$4
        READ_COUNT=$5
        MEAN_LENGTH=$6
        SD_LENGTH=$7
        READ_MUTATION_RATE=$8
        OUTFILE=$9
        SAMFILE=${10}
        ALGORITHM=${11}
        SORT_ON=${12}
        ;;

    *)
        echo "$(basename $0): expected 10 arguments!" >&2
        exit 1
        ;;
esac

tmp=/tmp/ids.$RANDOM

sam-to-fasta-alignment.py --referenceId genome-1 --samfile $SAMFILE | fasta-ids.py > $tmp

genome1Count=$(grep -c genome-1-read $tmp)
genome2Count=$(grep -c genome-2-read $tmp)
readCount=$(( genome1Count + genome2Count ))

rm $tmp

significant-base-frequencies.py \
    --homogeneousCutoff $HOMOGENEOUS_CUTOFF --minReads $MIN_READS --dropUnmapped \
    --title "Double infection ($GENOME_2_MUTATION_RATE genome separation rate) $ALGORITHM alignment of genome-1<br>$GENOME_LENGTH nt, $READ_COUNT reads from each genome (mean length $MEAN_LENGTH, sd $SD_LENGTH), $READ_MUTATION_RATE read mutation rate.<br>$readCount reads were aligned ($genome1Count from genome 1, $genome2Count from genome 2)" \
    --referenceId genome-1 \
    --titleFontSize 11 --axisFontSize 12 \
    --outfile $OUTFILE --samfile $SAMFILE $SORT_ON
