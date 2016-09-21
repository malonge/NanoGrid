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
ASMPREFIX=`cat asmprefix`
SCRIPT_PATH=`cat scripts`
READS=`cat reads`

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
jobid=$1
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

# now figure out which contig we are
NUM_JOBS=`wc -l $ASMPREFIX.fofn |awk '{print $1}'`
if [ $jobid -le 0 ]; then
   echo "Invalid job id, must be 1 or greater"
   exit
fi

if [ $jobid -gt $NUM_JOBS ]; then
   echo "Invalid job id, max is $NUM_JOBS"
   exit
fi

line=`cat $ASMPREFIX.fofn |head -n $jobid |tail -n 1`

if [ -e $ASMPREFIX.$jobid.fa ]; then 
   echo "Already done"
   exit 
fi

$SCRIPT_PATH/nanopolish/nanopolish variants --fix-homopolymers --consensus $ASMPREFIX.$jobid.fa -w $line -r $READS -b $PREFIX.sorted.bam -g $ASM -e $PREFIX.eventalign.sorted.bam -t 4 --min-candidate-frequency 0.1 --models $SCRIPT_PATH/nanopolish/etc/r9-models/nanopolish_models.fofn
