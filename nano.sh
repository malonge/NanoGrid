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
echo "Usage: nano.sh <assembly fasta> <reads.fasta>"
echo "This scrip will run nanopolish in parallel on the grid."
echo "It will index and sort alignments"

MACHINE=`uname`
PROC=`uname -p`
SCRIPT_PATH=$BASH_SOURCE
SCRIPT_PATH=`dirname $SCRIPT_PATH`
JAVA_PATH=$SCRIPT_PATH:.

ASM=$1
READS=$2
ASMPREFIX=`echo $1 |sed s/.fasta//g |sed s/.fna//g |sed s/.fa//g`
PREFIX=`echo $2 |sed s/.fasta//g |sed s/.fna//g |sed s/.fa//g`

syst=`uname -s`
arch=`uname -m`
name=`uname -n`

if [ "$arch" = "x86_64" ] ; then
  arch="amd64"
fi
if [ "$arch" = "Power Macintosh" ] ; then
  arch="ppc"
fi

if [ x$ASM == "x" ]; then
   echo "Error: you must specify an assembly fasta file"
   exit
fi

if [ x$READS == "x" ]; then
   echo "Error: you must specify reads fasta"
   exit
fi

if [ ! -e $ASM ]; then
   echo "Error: couldn't find $ASM, please try again"
   exit
fi
if [ ! -e $READS ]; then
   echo "Error: couldn't find $READS, please try again"
   exit
fi

# nanopolish only understand .fa extension, check for those
if [ "$ASM" != "$ASMPREFIX.fa" ]; then
   echo "Error: $ASM must end in .fa, nanopolish does not support other extensions"
   exit
fi
if [ "$READS" != "$PREFIX.fa" ]; then
   echo "Error: $READS must end in .fa, nanopolish does not support other extensions"
   exit
fi

NUM_CTG=`grep ">" $ASM |wc -l |awk '{print $1}'`
if [ $NUM_CTG -le 0 ]; then
   echo "Error: there are no contigs in the provided fasta file $ASM. Please try again"
   exit
fi

echo "$ASM" > asm
echo "$PREFIX" > prefix
echo "$ASMPREFIX" > asmprefix
echo "$SCRIPT_PATH" > scripts
echo "$READS" > reads

echo "Running with $PREFIX $ASM $READS mappings on $NUM_CTG contigs"

# split the mappings
python $SCRIPT_PATH/nanopolish/scripts/nanopolish_makerange.py  $ASM > $ASMPREFIX.fofn
NUM_JOBS=`wc -l $ASMPREFIX.fofn |awk '{print $1}'`

# now we can submit each range as an individual job and a merge job for the end
qsub -A ${ASMPREFIX}_nanopolish -V -pe thread 8 -l mem_free=2g -cwd -N "${ASMPREFIX}map" -j y -o `pwd`/map.out $SCRIPT_PATH/map.sh
qsub -A ${ASMPREFIX}_nanopolish -V -pe thread 4 -l mem_free=2g  -hold_jid "${ASMPREFIX}map" -t 1-$NUM_JOBS -cwd -N "${ASMPREFIX}nano" -j y  -o `pwd`/\$TASK_ID.polish.out $SCRIPT_PATH/nanoParallelSGE.sh
qsub -A ${ASMPREFIX}_nanopolish -V -pe thread 1 -l mem_free=10g -hold_jid "${ASMPREFIX}nano" -cwd -N "${ASMPREFIX}merge" -j y -o `pwd`/merge.out $SCRIPT_PATH/merge.sh
