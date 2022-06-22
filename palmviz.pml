#!/bin/pymol
#
# palmviz.pml 
# palmfold Visualization Script
#

## Input parameters to be read from
## command-line
##
#pdbin = '/home/ab/Desktop/palmfold/palmfold_v4/local2/00-1667.tar.gz.pdb'
pdbin = $PDBIN

#pngout = '/home/ab/Desktop/palmfold/palmfold_v4/output.png'
pngout = $PNGOUT

##
## These motif coordinates need to be
## adjusted to the input PDB (fold)
## taken from the structural alignment
Asrt = 1  # plus 11
Bsrt = 51 # plus 13
Csrt = 79 # plus 7

## Script
##
# Load folded PDB file
cmd.load(str(pdbin))

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
cmd.select("motifA", "resi " + str(Asrt) + "-" + str(Asrt+11))
cmd.select("motifB", "resi " + str(Bsrt) + "-" + str(Bsrt+13))
cmd.select("motifC", "resi " + str(Csrt) + "-" + str(Csrt+7))

color tv_blue, motifA
color tv_green, motifB
color tv_red, motifC

set cartoon_transparency, 0.2
## CHANGE OUTPUT VARIABLE
cmd.png(str(pngout), "10cm", "10cm", 300, 1)