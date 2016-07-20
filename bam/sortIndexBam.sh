#!/bin/bash

BAM=$1

echo "sort and index bam"
samtools sort -o ${BAM}.s.bam ${BAM}.bam 
samtools index ${BAM}.s.bam
