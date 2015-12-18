#!/bin/bash

BAM=$1

echo "sort and index bam"
samtools sort ${BAM}.bam ${BAM}.s 
samtools index ${BAM}.s.bam
