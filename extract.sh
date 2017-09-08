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
RAW=`cat raw`
READS=`cat reads`

if [ x$READS != "x" ] && [ -e $READS ]; then
   echo "Already done"
else
   if [ eventless ]; then
      $SCRIPT_PATH/nanopolish nanopolish index -d $RAW reads.fastq && echo "reads" > reads
   else
     $SCRIPT_PATH/nanopolish/nanopolish extract -r $RAW > reads.fa && echo "reads" > reads
     bgzip reads.fa
   fi
fi
