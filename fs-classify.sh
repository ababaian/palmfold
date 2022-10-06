#!/bin/bash

# fs-classify
#
# ./fs-classify -i <DB> -r <REF_DIR> -o <OUT_DIR> [OPTIONS]
set -eu
VERSION='0.1.0'

# USAGE
function usage {
	echo "FoldSeek palmprint classifier v$VERSION"
	echo
	echo "USAGE: ./fs-classify.sh -i <DB> -r <REF_DIR> -o <OUT_DIR> [OPTIONS]"
	echo
	echo " -h  Print this message"
	echo
	echo " [Required]"
	echo " -i  FoldSeek database to classify"
	echo
	echo " [Optional]"
	echo " -r  Directory to reference palmprint structures [./pol]"
	echo " -o  Output directory [.]"
	echo
	echo " [Miscellenaeous]"
	echo " -c  TM-score cutoff [0.6]"
	echo " -t  CPU threads to use [10]"
	echo " -m  Directory to write temporary files [/tmp]"
	echo
	exit 0
}

# Parse arguments
DIR=$(dirname $(realpath $0))
FDB=""
REF=$DIR/pol
OUT=$(pwd)
CUT="0.6"
CPU="10"
TMP="/tmp"

while getopts i:r:o:c:t:m:h OPT; do
	case $OPT in
		i)
			FDB=$OPTARG
			;;
		r)
			REF=$OPTARG
			;;
		o)
			OUT=$OPTARG
			;;
		c)
			CUT=$OPTARG
			;;
		t)
			CPU=$OPTARG
			;;
		m)
			TMP=$OPTARG
			;;
		h)
			usage
			;;
		\?)
			echo "Unknown option given"
			usage
			;;
	esac
done

# Check input
if [ -z "$FDB" ]
then
	echo "FoldSeek database not given (-i)"
	false
	exit 1
elif [ ! -f "$FDB" ]
then
	echo "FoldSeek database not found"
	false
	exit 1
elif [ ! -f "${FDB}_ca" ]
then
	echo "Invalid FoldSeek database given (C-alpha database not found)"
	false
	exit 1
fi

# Check palmdb
if [ ! -d "$REF/palmprint" ]
then
	echo "Could not find palmprint PDB directory ($REF/palmprint)"
	false
	exit 1
elif [ ! -f "$REF/rdrp.model.list" ] || [ ! -f "$REF/xdxp.model.list" ]
then
	echo "Could not find RdRP/XdXP list ($REF/[rdrp/xdxp].model.list)"
	false
	exit 1
fi

# Check output permission
mkdir -p $OUT
if [ ! -w "$OUT" ]
then
	echo "Unable to write on given output directory: $OUT"
	false
	exit 1
elif [ ! -w "$TMP" ]
then
	echo "Unable to write on given temp directory: $TMP"
	false
	exit 1
fi

# Run foldseek
foldseek easy-search $REF/palmprint $FDB $OUT/result.m8 $TMP --max-seqs 100000 --alignment-type 1 --tmalign-fast 0 --format-output query,target,evalue -e inf --tmscore-threshold $CUT --threads $CPU

# Parse foldseek results
cut -f2 $OUT/result.m8 | sort | uniq | while read ID
do
	RDTM=$(grep -f $REF/rdrp.model.list $OUT/result.m8 | grep $ID | cut -f3 | sort -rn | head -1)
	RDID=$(grep -f $REF/rdrp.model.list $OUT/result.m8 | grep $ID | sort -rnk3,3 | head -1 | cut -d. -f1)
	XDTM=$(grep -f $REF/xdxp.model.list $OUT/result.m8 | grep $ID | cut -f3 | sort -rn | head -1)

	if [ -z $XDTM ]
	then
		XDTM="0.0"
	fi
	
	if [ ! -z $RDTM ]
	then
		if [ $(echo "$RDTM >= $XDTM" | bc -l) ]
		then
			printf "%s\t%s\t%s\n" $ID $RDID $RDTM | tee -a $OUT/rdrp.list
		fi
	fi
done

