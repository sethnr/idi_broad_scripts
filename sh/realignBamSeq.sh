#!/bin/bash


## defaults
QUERY=1
REFERENCE=2
READFRAGSIZE=450
READPAIRS=60000 #~100x for Pf3D7 genome
READLENGTH=250
OUT="realignBam.out"

while getopts "q:r:f:t:o:" opt; do
  case $opt in
    q) QUERY=${OPTARG} ;;
    r) REFERENCE=$OPTARG ;;
    f) REFERENCE=$OPTARG ;;
    t) TMP=$OPTARG ;;
    o) OUT=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
    h) echo "  alignSimSeq -q <query_fasta> -r <ref_fasta>";
       echo "              -s <fragment_size:450> -l <read_length:250>";
       echo "              -n <no_read_pairs:60,000> -t <tmpdir:tmp>";;
  esac
done

#QUERYOUT=`basename $QUERY`
QUERYOUT=${OUT}.realign.bam
QOUT1=${QUERYOUT/%.bam/.1.fq}
QOUT2=${QUERYOUT/%.bam/.2.fq}
QUERYOUT=${QUERYOUT/%.bam/}

echo $QOUT1
echo $QOUT2
echo $REFERENCE


if [[ ! -f ${QUERYOUT}.S.bam ]];
then
    #samtools bam2fq /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bam
    echo "sorting"
    echo samtools sort -n $QUERY ${QUERYOUT}.S
    samtools sort -n $QUERY ${QUERYOUT}.S
    echo "sorted";
fi

if [[ ! -f ${QUERYOUT}.SFS.bam ]];
then
    #samtools bam2fq /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bam
    echo "fixing and sorting"
    samtools fixmate ${QUERYOUT}.S.bam ${QUERYOUT}.SF.bam 
    samtools sort -n ${QUERYOUT}.SF.bam ${QUERYOUT}.SFS
    echo "fixed and sorted";
fi

if [[ ! -f ${QOUT1} ]];
then
    echo "fastqing"
    echo bedtools bamtofastq  -i  ${QUERYOUT}.SFS.bam  -fq $QOUT1 -fq2 $QOUT2
    bedtools bamtofastq  -i  ${QUERYOUT}.SFS.bam  -fq $QOUT1 -fq2 $QOUT2
    echo "fastq'd";
fi

echo bwa mem -M -t 5 $REFERENCE  $QOUT1  $QOUT2 \| samtools view -bS - \> ${OUT}.bam
bwa mem -M -t 5 $REFERENCE  $QOUT1  $QOUT2 | samtools view -bS - > ${OUT}.bam


echo "sort and index bam"
samtools sort -o ${OUT}.s.bam -T temp  ${OUT}.bam
#mv ${OUT}.s.bam ${OUT}.bam
samtools index ${OUT}.s.bam
