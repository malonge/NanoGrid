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

if [ -e `pwd`/CONFIG ]; then
   CONFIG=`pwd`/CONFIG
else
   CONFIG=${SCRIPT_PATH}/CONFIG
fi
GRID=`cat $CONFIG |grep -v "#" |grep  GRIDENGINE |tail -n 1 |awk '{print $2}'`

if [ $GRID == "SGE" ]; then
   baseid=$SGE_TASK_ID
   offset=$1
elif [ $GRID == "SLURM" ]; then
   baseid=$SLURM_ARRAY_TASK_ID
   offset=$1
fi

if [ x$baseid = x -o x$baseid = xundefined -o x$baseid = x0 ]; then
  baseid=$1
  offset=0
fi

if [ x$offset = x ]; then
  offset=0
fi

jobid=`expr $baseid + $offset`

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

echo Running job $jobid based on command line options.

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

$SCRIPT_PATH/nanopolish/nanopolish variants --faster --consensus=$ASMPREFIX.$jobid.fa -w $line -r $READS -b $PREFIX.sorted.cram -g $ASM -t 16 --min-candidate-frequency 0.01 --fix-homopolymers
