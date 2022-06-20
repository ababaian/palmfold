#!/bin/pymol
#
# palmviz.pml 
# palmfold Visualization Script
#

## Input parameters to be read from
## command-line
##
#pdbin  = '/home/ab/Desktop/palmfold/palmfold_v4/palmdb_pilot/pdb_realign/1685_unrelaxed_rank_1_model_3.pdb_realign.pdb'
#pngout = '/home/ab/Desktop/palmfold/palmfold_v4/palmdb_pilot/output.png'
##
## These motif coordinates need to be
## adjusted to the input PDB (fold)
## taken from the structural alignment
#Astart = 1
#Aend   = 12
#Bstart = 51
#Bend   = 64
#Cstart = 79
#Cend   = 86

## Script
##
# Load folded PDB file
## CHANGE INPUT VARIABLE
#load /home/ab/Desktop/palmfold/palmfold_v4/palmdb_pilot/pdb_realign/1685_unrelaxed_rank_1_model_3.pdb_realign.pdb
load $PDBIN

# Set background cartoon
# to white
select /*
hide
show cartoon
color gray40, sele

# Use standard view angle angle
set_view (\
    -0.891759515,    0.449078411,    0.055630598,\
     0.429066062,    0.800084352,    0.419247121,\
     0.143766090,    0.397736073,   -0.906164587,\
     0.000000000,    0.000000000, -173.986724854,\
    19.074066162,  -16.383844376,  -10.496936798,\
   154.656341553,  193.317108154,  -20.000000000 )


# Create motif selections
## CHANGE VARIABLES
select motifA, resi 1-12
select motifB, resi 51-64
select motifC, resi 84-91

color tv_blue, motifA
color tv_green, motifB
color tv_red, motifC

set cartoon_transparency, 0.2
## CHANGE OUTPUT VARIABLE
#png /home/ab/Desktop/palmfold/palmfold_v4/palmdb_pilot/output.png, 10cm, dpi=300, ray=1
png $PNGOUT, 10cm, dpi=300, ray=1