#!/bin/bash
# palmfold
#
# Structural alignment and scoring against
# a reference set of palmprint structures
#
VERSION='0.2.0'
#
set -eu

# Usage
function usage {
  echo "palmfold v $VERSION"
  echo ""
  echo "Usage: ./palmfold.sh -p <path_to_palmprint> -d <path_to_test_pdb> -o <output_dir> [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Fields"
  echo "    -o    Output files path (created if DNE)"
  echo ""
  echo "    Reference Palmprint and Test Structures"
  echo "    -p    Palmprint Reference Path (pdb) [./pol]"
  echo "            must contain palmprint/*.pdb & rdrp.model.list file"
  echo "    -d    Predicted Fold PDB-file Path [./pdb]"
  echo ""
  echo "    TMalign Parameters"
  echo "    -s    TMalign Cut-off threshold for inclusion [0.5]"
  echo ""
  echo "ex: ./palmfold.sh -p ./pol -d ./pdb -o test_run"
  false
  exit 1
}

# Default parameters
PALMPRINTS="./pol"
PDBS="./pdb"
OUTNAME='' # unset
CUTOFF="0.5"
SPATH=$(dirname $(realpath $0))

while getopts p:d:o:s:h FLAG; do
  case $FLAG in
    p)
      PALMPRINTS=$OPTARG
      ;;
    d)
      PDBS=$OPTARG
      ;;
    o)
      OUTNAME=$OPTARG
      ;;
    s)
      CUTOFF=$OPTARG
      ;;
    h)  #show help ----------
      usage
      ;;
    \?) #unrecognized option - show help
      echo "Input parameter not recognized"
      usage
      ;;
  esac
done
shift $((OPTIND-1))

# Check Required
if [ -z "$OUTNAME" ]; then
    echo "Output directory name not set (-o)"
    false
    exit 1
fi

# Initialize workspace
mkdir -p $OUTNAME
mkdir -p $OUTNAME/pdb_realign
mkdir -p $OUTNAME/png_realign
mkdir -p $OUTNAME/pdb_realign2
mkdir -p $OUTNAME/png_realign2
mkdir -p $OUTNAME/tmfa
mkdir -p $OUTNAME/fa
mkdir -p $OUTNAME/fa/pp
mkdir -p $OUTNAME/fa/rc

# tmp directory
mkdir -p $OUTNAME/tmp

# Initialize TMalign header output
if [ -f $OUTNAME/result.tm ]; then
  echo "Overwriting result files."
fi
echo -e 'PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali' > $OUTNAME/result.tm

# Run TMalign against Reference Palmprints
# =========================================================
for pdbz in $(ls $PDBS/); do

  # If PDB is gz compressed, decompress in place
  GZ='FALSE'
  if [[ $pdbz == *.gz ]]; then
    GZ='TRUE'
    gzip -d $PDBS/$pdbz
    pdb=$(echo $pdbz | sed 's/.gz//g' -)
  else
    pdb=$pdbz
  fi

  # Input PDB File
  #echo $pdb

  # Iterate through reference palmprint
  for pp in $(ls $PALMPRINTS/palmprint/); do
    TMalign -outfmt 2 \
      $PDBS/$pdb \
      $PALMPRINTS/palmprint/$pp \
      >> $OUTNAME/tmp/pdb_raw.tm
  done

  echo DONE SECTION 1

  # Clean-up TM output
  grep '.pdb' $OUTNAME/tmp/pdb_raw.tm \
   | sed 's/\.pdb//g' - \
   | sed "s,$PDBS/,,g" - \
   | sed "s,$PALMPRINTS/palmprint/,,g" - \
   > $OUTNAME/tmp/pdb_clean.tm

  rm $OUTNAME/tmp/pdb_raw.tm

  # Append TMalign CSV to output file
  cat $OUTNAME/tmp/pdb_clean.tm >> $OUTNAME/result.tm

  # Isolate Maximum RdRP TMalign Score
  maxRdRP=$(grep -f $PALMPRINTS/rdrp.model.list $OUTNAME/tmp/pdb_clean.tm \
    | cut -f 2,4 | sort -k 2 -nr -| head -n1)

    maxRdRP_model=$(echo $maxRdRP | cut -d' ' -f 1)
    maxRdRP_score=$(echo $maxRdRP | cut -d' ' -f 2)

  # Isolate Maximum XdXP TMalign Score (NOT RdRP)
  maxXdXP=$(grep -v -f $PALMPRINTS/xdxp.model.list $OUTNAME/tmp/pdb_clean.tm \
    | cut -f 2,4 | sort -k 2 -nr -| head -n1)

    maxXdXP_model=$(echo $maxXdXP | cut -d' ' -f 1)
    maxXdXP_score=$(echo $maxXdXP | cut -d' ' -f 2)

  # Use bc to compare floats
  # Does maxRdRP_score pass CUTOFF value
  if [ 1 -eq "$(echo "$maxRdRP_score >= $CUTOFF" | bc)" ]; then
    
    # A significant RdRP match is present
    # Does maxRdRP surpass maxXdXP

    if [ 1 -eq "$(echo "$maxRdRP_score >= $maxXdXP_score" | bc)" ]; then
      # RdRP Hit is significant and surpasses XdXP

      # Generate Fasta and Re-Align PDB output
      # against TOP HIT only
      
      TMalign -outfmt 1 \
      $PDBS/$pdb \
      $PALMPRINTS/palmprint/$maxRdRP_model.pdb \
      -o $OUTNAME/tmp/realign \
      > $OUTNAME/tmfa/$pdb.fa

      mv $OUTNAME/tmp/realign.pdb $OUTNAME/pdb_realign/$pdb

      # Render PNG output of PDB
      sed "s,\$PDBIN,$OUTNAME/pdb_realign/$pdb,g" palmviz.pml \
        | sed "s,\$PNGOUT,$OUTNAME/png_realign/$pdb.png,g" \
        > tmp.pml

      pymol -c -Q -r tmp.pml
      rm tmp.pml

      # PROCESS TM FASTA FILE TO ISOLATE
      # PALMPRINT AND RDRPCORE
      # python3 palmgrab.py -i <input.tm.fa> -p <palmprint.fa> -r <rdrpcore.fa>
      python3 $SPATH/palmgrab.py $OUTNAME/tmfa/$pdb.fa \
                          $OUTNAME/tmp/$pdb.pp.fa \
                          $OUTNAME/tmp/$pdb.rc.fa

      mv $OUTNAME/tmp/$pdb.pp.fa $OUTNAME/fa/pp/
      mv $OUTNAME/tmp/$pdb.rc.fa $OUTNAME/fa/rc/

      echo -e "$pdb\tPositive: $maxRdRP_model $maxRdRP_score"
    else echo -e "$pdb\tNegative: XdXP score is higher"

      TMalign -outfmt 1 \
        $PDBS/$pdb \
        $PALMPRINTS/palmprint/1hhs_A-cyst.pdb \
        -o $OUTNAME/tmp/realign \
        > $OUTNAME/tmfa/$pdb.fa

      mv $OUTNAME/tmp/realign.pdb $OUTNAME/pdb_realign2/$pdb

      # Render PNG output of PDB
      sed "s,\$PDBIN,$OUTNAME/pdb_realign2/$pdb,g" palmviz.pml \
        | sed "s,\$PNGOUT,$OUTNAME/png_realign2/$pdb.png,g" \
        > tmp.pml

    fi
  else echo -e "$pdb\tNegative: Score below cutoff"

    TMalign -outfmt 1 \
      $PDBS/$pdb \
      $PALMPRINTS/palmprint/1hhs_A-cyst.pdb \
      -o $OUTNAME/tmp/realign \
      > $OUTNAME/tmfa/$pdb.fa

    mv $OUTNAME/tmp/realign.pdb $OUTNAME/pdb_realign2/$pdb

    # Render PNG output of PDB
    sed "s,\$PDBIN,$OUTNAME/pdb_realign2/$pdb,g" palmviz.pml \
      | sed "s,\$PNGOUT,$OUTNAME/png_realign2/$pdb.png,g" \
      > tmp.pml

  fi

  # Recompress, if needed
  if [[ "$GZ" = 'TRUE' ]]; then
    gzip $PDBS/$pdb
  fi

done

# Create merged fasta outputs
cat /dev/null > $OUTNAME/palmprints.fa
cat /dev/null > $OUTNAME/rdrpcores.fa
ls $OUTNAME/fa/pp/ | while read PP; do
  cat $OUTNAME/fa/pp/$PP >> $OUTNAME/palmprints.fa
done
ls $OUTNAME/fa/rc/ | while read RC; do
  cat $OUTNAME/fa/rc/$RC >> $OUTNAME/rdrpcores.fa
done

# Clean-up temporary directory
#rm -rf $OUTNAME/tmp
