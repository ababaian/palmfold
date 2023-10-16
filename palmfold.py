#!/bin/python3
# palmfold
#
# RNA virus RdRp structural classifier.
#
import argparse
import logging as log
import gzip
import shutil
from pathlib import Path
from palmfoldutils import PalmStructs
# Initialization Routine ==================================
def main(inputpath:Path, outpath:Path, palmprints:Path, threshold:float):
    # Verify Input Path exists
    if not inputpath.is_dir():
        log.error("The input directory %s does not exist", inputpath)
        exit(1)

    # Verify Output Path exists
    if not outpath.is_dir():
        log.info("Output directory %s does not exist. Creating", outpath)
        outpath.mkdir()
    else:
        log.info("Output directory %s exists. Overwriting", outpath)
        
    # Extract protein filenames
    names = [
        filename.stem
        for filename in inputpath.iterdir()
        if filename.suffix==".pdb"
    ]

    namesgz = [
        filename.name[:-7]
        for filename in inputpath.iterdir()
        if filename.name.endswith(".pdb.gz")
    ]

    log.info(" %s %s", names, namesgz)
    
    # Verify Input Path contains PDB files
    if not names and not namesgz:
        log.error("Input directory %s doesn't contain any .pdb(.gz) files.", inputpath)
        exit(1)
    
    # Create the Palmprint datastructure
    log.info("Initializing palmfold...")
    ps = PalmStructs(palmprints)

    # Classify the input protein structures
    for protgz in namesgz:
        log.info("")
        log.info("  Analyzing %s...", protgz)
        # Decompress gz
        # TODO should also use `TemporaryDirectory`
        with gzip.open(inputpath/f"{protgz}.pdb.gz", 'rb') as f_in:
            with open(inputpath/f"{protgz}.pdb", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        ps.align(inputpath, outpath, protgz, threshold)
        (inputpath/f"{protgz}.pdb").unlink()

    for prot in names:
        log.info("")
        log.info("  Analyzing %s...", prot)
        ps.align(inputpath, outpath, prot, threshold)

# Help / Argument Parsing =================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='palmfold: RNA virus RdRp structural classifier.'
        )
    parser.add_argument(
        '--inputpath', '-i', type=Path, default='./pdb',
        help='Directory containing "*.pdb(.gz)"" files to be classified. [./pdb]'
        )
    parser.add_argument(
        '--palmprints', '-p', type=Path, default=Path(__file__).parent/'pol',
        help='path to "~/palmfold/pol/" containing palmprint pdb models. [./pol]'
        )
    ## TODO: implement output directory
    parser.add_argument(
        '--outpath', '-o', type=Path, default='./out',
        help='Directory for writing output files. [./out]'
        )
    parser.add_argument(
        '--threshold', '-t', type=float, default=0.4,
        help='Polymerase+ re-alignment threshold for TMAlign. [0.4]'
        )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Print INFO messages, else prints WARN and ERROR only'
        )
    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
        log.info("Verbose output.")
    else:
        log.basicConfig(format="%(levelname)s: %(message)s")

    # Main script routine =================================
    log.info('''
      _____      _            __      _     _ 
     |  __ \    | |          / _|    | |   | |
     | |__) |_ _| |_ __ ___ | |_ ___ | | __| |
     |  ___/ _` | | '_ ` _ \|  _/ _ \| |/ _` |
     | |  | (_| | | | | | | | || (_) | | (_| |
     |_|   \__,_|_|_| |_| |_|_| \___/|_|\__,_|
     =========================================
    ''')
    log.info('''
          #
        @@/%##&
        *########%
        **#######*
        **##   ##*  --=  ####*
        **## P ##*    (****#####*
        **## A ##*      ,***&######@===-#
        **## L ##*      /  #########     @
        **## M ##*          #######       #
        **## F ##*                       @
        .*## O ##*     ___--__       @#
         *## L ##*    @       *%      %
         *## D ##*   %          #       #
         *##   ##*    @@         @       %
         *#######*     *       *@@@@      **
         *#######*    *#*      *@@@@    *&@@*
         *#######*   @###*     *@@@@   **&@@@*
         *#######*  ######*    *@@@@    *&@@@
         *#######*   ####*     *@@@@    *@@@@
        .*#######*   ####*     *@@@@    *@@@@
        @*########   ####*    **####*   *@@@@ 
        @*#######@     (       ****       @@@ 
        @*#######@     @         *         @
        @*#######@    *                    #
         *#######*#   *                   #
          %####& /#    #
              @         #
               @       #
                *###/##
        ''')

    main(args.inputpath, args.outpath, args.palmprints, args.threshold)
    log.info("")
    log.info("palmfold run completed successfully!")
