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
READS=`cat readsextracted`

if [ -e $PREFIX.sorted.cram ]; then
   echo "Already done"
else
   if [ ! -e $READS ]; then
      echo "Error no reads found to map, check the output of extract.out"
      exit
   fi
i
   if [ ! -e $ASM.mmi ]; then
      minimap2 -ax map-ont -d $ASM.mmi $ASM && touch $PREFIX.index.success
   fi
   if [ ! -e $PREFIX.map.success ]; then
      minimap2 -ax map-ont -t 16 --secondary=no $ASM.mmi $READS > $PREFIX.sam && touch $PREFIX.map.success)
      if [ -e $PREFIX.map.success ]; then
         samtools sort -@16 -O cram -o $PREFIX.sorted.cram -T $PREFIX.tmp --reference=$ASM $PREFIX.sam
         samtools index $PREFIX.sorted.cram
      fi
   fi
fi

# run the first mapping job to make sure we have all files needed created, otherwise there is a race condition
sh  $SCRIPT_PATH/nanoParallelSGE.sh 1
