#!/bin/bash

# fs-classify
#
# ./fs-classify -i <DB> -r <REF_DIR> -o <OUT_DIR> [OPTIONS]
set -eu
VERSION='0.1.0'

# USAGE
function usage {
	echo "Foldseek palmprint classifier v$VERSION"
	echo
	echo "USAGE: ./fs-classify.sh -i <DB> -r <REF_DIR> -o <OUT_DIR> [OPTIONS]"
	echo
	echo " -h  Print this message"
	echo
	echo " [Required]"
	echo " -i  Foldseek database to classify"
	echo
	echo " [Optional]"
	echo " -r  Directory to reference palmprint structures [./pol]"
	echo " -o  Output directory [.]"
	echo
	echo " [Miscellenaeous]"
	echo " -c  TM-score cutoff [0.6]"
	echo " -x  Maximum results per query sequence allowed to pass the foldseek prefilter (affects sensitivity) [100000]"
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
XSQ="100000"
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
		x)
			XSQ=$OPTARG
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
	echo "Foldseek database not given (-i)"
	echo "USAGE: ./fs-classify.sh -i <DB> -r <REF_DIR> -o <OUT_DIR> [OPTIONS]"
	exit 1
elif [ ! -f "$FDB" ]
then
	echo "Foldseek database not found"
	false
	exit 1
elif [ ! -f "${FDB}_ca" ]
then
	echo "Invalid Foldseek database given (C-alpha database not found)"
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
foldseek easy-search $REF/palmprint $FDB $OUT/result.m8 $TMP --max-seqs $XSQ --alignment-type 1 --tmalign-fast 0 --format-output query,target,evalue -e inf --tmscore-threshold $CUT --threads $CPU

# Parse foldseek results
printf "%s\t%s\t%s\t%s\t%s\n" "ID" "model_XdXP" "score_XdXP" "model_RdRP" "score_RdRP" > $OUT/rdrp.list

cut -f2 $OUT/result.m8 | sort | uniq | while read ID
do
	RDTM=$(grep -f $REF/rdrp.model.list $OUT/result.m8 | grep $ID | cut -f3 | sort -rn | head -1)
	RDID=$(grep -f $REF/rdrp.model.list $OUT/result.m8 | grep $ID | sort -rnk3,3 | head -1 | cut -d. -f1)
	XDTM=$(grep -f $REF/xdxp.model.list $OUT/result.m8 | grep $ID | cut -f3 | sort -rn | head -1)
	XDID=$(grep -f $REF/xdxp.model.list $OUT/result.m8 | grep $ID | sort -rnk3,3 | head -1 | cut -d. -f1)

	if [ -z $XDTM ]
	then
		XDID="N/A"
		XDTM="0.0"
	fi
	
	if [ ! -z $RDTM ]
	then
		if [ $(echo "$RDTM >= $XDTM" | bc -l) -eq 1 ]
		then
			printf "%s\t%s\t%s\t%s\t%s\n" $ID $XDID $XDTM $RDID $RDTM | tee -a $OUT/rdrp.list
		fi
	fi
done

