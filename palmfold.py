#!/bin/python3
# palmfold
#
# RNA virus RdRp structural classifier.
#
from os import path, mkdir, listdir, rename, remove
from sys import stderr

import argparse
import subprocess
import logging as log
import gzip
import shutil

class PalmStructs:

    # Initialization Routine ------------------------------
    # Ensures ./pol dir-tree is populated
    # with the correct files (palmprint pdb, pdb, fasta)
    def __init__(self, polpath):
        # Check if the .pol directory exists
        if not path.exists(polpath):
            log.error("The palmprint directory %s does not exist", polpath)
            exit(1)
        self.polpath = polpath
        
        # Verify the reference RdRp models and fasta files exist
        self.rdrps = []
        rdrp_filelist = path.join(polpath, "rdrp.model.list")

        if not path.exists(rdrp_filelist):
            log.error('The "rdrp.model.list" file not found in %s,', polpath)
            exit(1)

        # Check each input RdRp model fasta, full_pdb, and palmprint_pdb
        with open(rdrp_filelist) as rdrp_fl:
            log.info("  Verifying RdRp Model List")
            # Verify each rdrp file in directory tree
            for line in rdrp_fl:
                line = line.strip()
                if len(line) > 0:
                    name = line
                    log.info("    %s", name)
                    # verify fasta
                    if not path.exists(path.join(polpath, "fa", f"{name}.fa")):
                        log.warning("Fasta file missing for: %s", name)
                        continue
                    # verify full sequence
                    gz = path.join(polpath, "full_length", f"{name}.pdb.gz")
                    if not path.exists(gz):
                        log.warning("Full length structure missing for: %s", name)
                        continue
                    # verify palmprint pdb structucture (required)
                    if not path.exists(path.join(polpath, "palmprint", f"{name}.pdb")):
                        log.error("<palmprint>.pdb missing for: %s", name)
                        exit(1)
                    self.rdrps.append(name)


        # Verify the reference XdXp models and fasta files exist
        self.xdxps = []
        xdxp_filelist = path.join(polpath, "xdxp.model.list")

        if not path.exists(xdxp_filelist):
            log.error('The "xdxp.model.list" file not found in %s,', polpath)
            exit(1)

        # Check each input XdXp model fasta, full_pdb, and palmprint_pdb
        with open(xdxp_filelist) as xdxp_fl:
            log.info("  Verifying XdXp Model List")
            # Verify each XdXp file in directory tree
            for line in xdxp_fl:
                line = line.strip()
                if len(line) > 0:
                    name = line
                    log.info("    %s", name)
                    # verify fasta
                    if not path.exists(path.join(polpath, "fa", f"{name}.fa")):
                        log.warning("Fasta file missing for: %s", name)
                        continue
                    # verify full sequence
                    gz = path.join(polpath, "full_length", f"{name}.pdb.gz")
                    if not path.exists(gz):
                        log.warning("Full length structure missing for: %s", name)
                        continue
                    # verify palmprint pdb structucture (required)
                    if not path.exists(path.join(polpath, "palmprint", f"{name}.pdb")):
                        log.error("<palmprint>.pdb missing for: %s", name)
                        exit(1)
                    self.xdxps.append(name)

        self.all_domains = self.rdrps + self.xdxps

    # Classification Routine ------------------------------
    # Runs TMalign of input against all palmprint pdb
    # Extracts scores and re-aligns if a RdRp palmprint detected
    def align(self, pdbpath, outpath, name, tm_threshold):
        # TODO: this looks like an inefficient way to find files in path
        # it introduces a bug where two files "ABC.pdb" and "ABCD.pdb" will
        # have the same start and endings, and cause multiple returns and skip
        # of ABC.pdb
        pdb_file = [file for file in listdir(pdbpath) if file.startswith(f"{name}") and file.endswith(".pdb")]
        if len(pdb_file) == 0:
            log.warning("No <input>.pdb file found  %s", name)
            return
        elif len(pdb_file) > 1:
            log.warning("Multiple <input>.pdb files found for %s", name)
            log.warning("Abiguous situation, skip molecule")
            return

        # Define input/output files
        out_name = path.join(outpath, path.splitext(pdb_file[0])[0])
        out_tsv = path.join(out_name + ".tm")
        out_ppfa = path.join(out_name + ".pp.fa")
        out_rcfa = path.join(out_name + ".rc.fa")
        out_rpdb = path.join(out_name + "_realigned.pdb")
        pdb_file = path.join(pdbpath, pdb_file[0])

        # Reset maximum scores observed
        max_rdrp_score = max_xdxp_score = 0
        max_rdrp = max_xdxp = None

        # Initialize output TSV file
        log.info("  Writing tsv output to %s", out_tsv)
        with open(out_tsv, "w") as out:
            print(
                "PDBchain1\tPDBchain2\tTM1\tTM2\t" +
                "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali",
                file=out
                )

            # TODO: these should be wrapped in functions to not repeat code
            # Run RdRp-palmprint TMalign and identify max_RdRp_Score
            log.info("  Run TMalign against RdRp palmprints")
            for domain in self.rdrps:
                # TMalign command
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')}"
                log.info("    %s", tmalign_cmd)

                # Runs TMalign and captures output
                ret_val = subprocess.run(
                    tmalign_cmd.split(" "),
                    capture_output=True,
                    text=True,
                    close_fds=False)
                # Parse the captured output for scores
                values = ret_val.stdout.split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
                # Append TMalign output to TSV output
                print("\t".join([str(v) for v in values]), file=out)
                # Save best scores
                score = values[3]
                if score >= max_rdrp_score:
                    max_rdrp_score = score
                    max_rdrp = domain

            # Run XdXp-palmprint TMalign and identify max_XdXp_Score
            log.info("  Run TMalign against XdXp palmprints")
            for domain in self.xdxps:
                # Run TM-Align
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')}"
                log.info("    %s", tmalign_cmd)
                
                # Runs TMalign and captures output
                ret_val = subprocess.run(
                    tmalign_cmd.split(" "),
                    capture_output=True,
                    text=True,
                    close_fds=False)
                # Parse the captured output for scores
                values = ret_val.stdout.split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
                # Append TMalign output to TSV output
                print("\t".join([str(v) for v in values]), file=out)
                # Save best scores
                score = values[3]
                if score >= max_xdxp_score:
                    max_xdxp_score = score
                    max_xdxp = domain

        # Does maxRdRP_score pass CUTOFF value
        log.info("")
        log.info("  Max RdRp palmprint: %s - %s", max_rdrp, max_rdrp_score)
        log.info("  Max XdXp palmprint: %s - %s", max_xdxp, max_xdxp_score)
        log.info("")

        # RdRp classification
        # RdRp TMscore above re-alignment threshold
        if max(max_rdrp_score, max_xdxp_score) >= tm_threshold:
            log.info("  + Polymerase re-alignment threshold score reached: %s", tm_threshold)
            if max_rdrp_score > max_xdxp_score:
                # RdRp+ classification
                log.info("  ++ RdRp-classification for %s", name)
            else:
                log.info("  %s classified as non-RdRp polymerase", name)

            log.info("")
            log.info("  Re-aligning:")

            # Run TMalign to re-align input structure to top-RdRp palmprint
            realign_cmd = f"TMalign -outfmt 1 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')} -o {pdb_file}_tmpalign"
            log.info("  %s", realign_cmd)
            with open(f"{pdb_file}.tmp", "w") as output:
                subprocess.run(realign_cmd.split(" "), stdout=output)

            # Move the realign
            rename(f"{pdb_file}_tmpalign.pdb", f"{out_rpdb}")

            # Run palmgrab.py to sub-select palmprint/core fasta seq
            log.info("")
            log.info("  Run palmgrab.py to extract sub-sequence fasta")
            log.info("  cmd:")
            scriptdir = path.dirname(path.realpath(__file__))
            log.info('  python3 ' + path.join(scriptdir, "palmgrab.py ") + f"{pdb_file}.tmp " + f"{out_ppfa} " + f"{out_rcfa}")

            subprocess.run(['python3', path.join(scriptdir, "palmgrab.py"), f"{pdb_file}.tmp", f"{out_ppfa}", f"{out_rcfa}"])
            log.info("  done")
            log.info("")

            for f in listdir(pdbpath):
                if "_tmpalign" in f:
                    remove(path.join(pdbpath, f))
            remove(f"{pdb_file}.tmp")


# Initialization Routine ==================================
def main(inputpath, outpath, palmprints, threshold):
    # Verify Input Path exists
    if not path.exists(inputpath):
        log.error("The input directory %s does not exist", inputpath)
        exit(1)

    # Verify Output Path exists
    if not path.exists(outpath):
        log.info("Output directory %s does not exist. Creating", outpath)
        mkdir(outpath)
    else:
        log.info("Output directory %s exists. Overwriting", outpath)
        
    # Extract protein filenames
    names = [
        filename[:-4]
        for filename in listdir(inputpath)
        if filename.endswith(".pdb")
    ]

    namesgz = [
        filename[:-7]
        for filename in listdir(inputpath)
        if filename.endswith(".pdb.gz")
    ]

    log.info(" %s %s", names, namesgz)

    # Verify Input Path contains PDB files
    if not names:
        if not namesgz:
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
        with gzip.open(path.join(inputpath, protgz + ".pdb.gz"), 'rb') as f_in:
            with open(path.join(inputpath, protgz + ".pdb"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        ps.align(inputpath, outpath, protgz, threshold)
        remove(path.join(inputpath, protgz + ".pdb"))

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
        '--inputpath', '-i', type=str, default='./pdb',
        help='Directory containing "*.pdb(.gz)"" files to be classified. [./pdb]'
        )
    parser.add_argument(
        '--palmprints', '-p', type=str, default='./pol',
        help='path to "~/palmfold/pol/" containing palmprint pdb models. [./pol]'
        )
    ## TODO: implement output directory
    parser.add_argument(
        '--outpath', '-o', type=str, default='./out',
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

    main(args.inputpath, args.outpath, args.palmprints, args.threshold)
    log.info("")
    log.info("palmfold run completed successfully!")
