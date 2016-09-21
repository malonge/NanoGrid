#!/usr/bin/env bash

######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

ASM=`cat asm`
PREFIX=`cat prefix`
SCRIPT_PATH=`cat scripts`
READS=`cat reads`

if [ -e $PREFIX.eventalign.sorted.bam ]; then
   echo "Already done"
   exit
fi

if [ ! -e $ASM.sa ]; then
   bwa index $ASM && touch $PREFIX.index.success
fi
if [ ! -e $PREFIX.map.success ]; then
   (bwa mem -x ont2d -t 8 $ASM $READS && touch $PREFIX.map.success) | samtools view -Sb - | samtools sort  - -o $PREFIX.sorted.bam
   samtools index $PREFIX.sorted.bam
fi

if [ ! -e $PREFIX.eventalign.success ]; then
   # Align the reads in event space
   ($SCRIPT_PATH/nanopolish/nanopolish eventalign -t 8 --sam -r $READS -b $PREFIX.sorted.bam -g $ASM --models $SCRIPT_PATH/nanopolish/etc/r9-models/nanopolish_models.fofn && touch $PREFIX.eventalign.success) | samtools view -Sb - | samtools sort - -o $PREFIX.eventalign.sorted.bam
   samtools index $PREFIX.eventalign.sorted.bam
fi

# run the first mapping job to make sure we have all files needed created, otherwise there is a race condition
sh  $SCRIPT_PATH/nanoParallelSGE.sh 1
