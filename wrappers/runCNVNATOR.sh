#!/bin/bash

#fasta files for each chromosome (suffix .fa)
REFS=/seq/plasmodium/sredmond/refs/Pf3D7_v3_chrs/
BINSIZE=200

while getopts "B:b:R:l:G:p:o:m:M:V" opt; do
  case $opt in
    b) BINSIZE=${OPTARG} ;;
    R) REFS=${OPTARG};;
    o) OUTFILE=$OPTARG;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done
shift $(expr $OPTIND - 1 )

#echo $@

echo "dumping aligns"
echo cnvnator -root ${OUTFILE}.tree  -tree $@
cnvnator -root ${OUTFILE}.tree  -tree $@

echo "calculating histogram"
echo cnvnator -root ${OUTFILE}.tree -his $BINSIZE -d $REFS
cnvnator -root ${OUTFILE}.tree -his $BINSIZE -d $REFS

echo "getting stats"
echo cnvnator -root ${OUTFILE}.tree -stat $BINSIZE -d $REFS
cnvnator -root ${OUTFILE}.tree -stat $BINSIZE -d $REFS

echo "partitioning coverages"
echo cnvnator -root ${OUTFILE}.tree -partition $BINSIZE -d $REFS
cnvnator -root ${OUTFILE}.tree -partition $BINSIZE -d $REFS

echo "calling CNVs"
echo cnvnator -root ${OUTFILE}.tree -call $BINSIZE -d $REFS
cnvnator -root ${OUTFILE}.tree -call $BINSIZE -d $REFS > ${OUTFILE}_CNVs
