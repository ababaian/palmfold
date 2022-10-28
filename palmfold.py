from os import path, mkdir, listdir, rename, remove
from sys import stderr

import argparse
import subprocess
import logging as log

class PalmStructs:

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

        # Check each input RdRp model fasta, full_pdb, and palmprint_pdb
        with open(xdxp_filelist) as xdxp_fl:
            log.info("  Verifying RdRp Model List")
            # Verify each xdxp file in directory tree
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

    def align(self, directory, name, out_tsv, rdrp_threshold):
        pdb_file = [file for file in listdir(directory) if file.startswith(f"{name}") and file.endswith(".pdb")]
        if len(pdb_file) == 0:
            print(f"No pdb file found for molecule {name}", file=stderr)
            return
        elif len(pdb_file) > 1:
            print(f"Multiple pdb file found for molecule {name}.", file=stderr)
            print(f"Abiguous situation, skip molecule", file=stderr)
            return
        pdb_file = path.join(directory, pdb_file[0])

        max_rdrp_score = max_xdxp_score = 0
        max_rdrp = max_xdxp = None
        with open(out_tsv, "w") as out:
            print("PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali", file=out)
            # Serch for the max score with rdrp proteins
            for domain in self.rdrps:
                # Run TM-Align
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')}"
                ret_val = subprocess.run(tmalign_cmd.split(" "), capture_output=True, text=True, close_fds=False)
                # Parse the output
                values = ret_val.stdout.split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
                print("\t".join([str(v) for v in values]), file=out)
                # Save best scores
                score = values[3]
                if score >= max_rdrp_score:
                    max_rdrp_score = score
                    max_rdrp = domain

            # Serch for the max score with xdxp proteins
            for domain in self.xdxps:
                # Run TM-Align
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')}"
                ret_val = subprocess.run(tmalign_cmd.split(" "), capture_output=True, text=True, close_fds=False)
                # Parse the output
                values = ret_val.stdout.split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
                print("\t".join([str(v) for v in values]), file=out)
                # Save best scores
                score = values[3]
                if score >= max_xdxp_score:
                    max_xdxp_score = score
                    max_xdxp = domain

        # Does maxRdRP_score pass CUTOFF value
        if max_rdrp_score >= rdrp_threshold:
            # Does maxRdRP surpass maxXdXP
            if max_rdrp_score > max_xdxp_score:
                # Realign the protein with the ref
                realign_cmd = f"TMalign -outfmt 1 {pdb_file} {path.join(self.polpath, 'palmprint', f'{domain}.pdb')} -o {pdb_file}_tmpalign"
                with open(f"{pdb_file}.tmp", "w") as output:
                    subprocess.run(realign_cmd.split(" "), stdout=output)
                # Get the realign
                rename(f"{pdb_file}_tmpalign.pdb", f"{pdb_file}_realign.pdb")
                # Get the fastas
                scriptdir = path.dirname(path.realpath(__file__))
                subprocess.run(['python3', path.join(scriptdir, "palmgrab.py"), f"{pdb_file}.tmp", f"{pdb_file}.pp.fa", f"{pdb_file}.rc.fa"])
                for f in listdir(directory):
                    if "_tmpalign" in f:
                        remove(path.join(directory, f))
                remove(f"{pdb_file}.tmp")
                remove(pdb_file)


# Main script routine =====================================
def main(inputpath, palmprints, threshold):
    # Verify Input Path exists
    if not path.exists(inputpath):
        log.error("The input directory %s does not exist", inputpath)
        exit(1)

    # Extract protein filenames
    # TODO add support for .pdb and .pdb.gz files (gzip compatibility)
    names = [filename[:-4] for filename in listdir(inputpath) if filename.endswith(".pdb")]
    
    # Verify Input Path contains PDB files
    if not names:
        log.error("Input directory %s doesn't contain any .pdb files.", inputpath)
        exit(1)
    
    # Create the Palmprint datastructure
    log.info("Initializing palmfold...")
    ps = PalmStructs(palmprints)

    # Classify the input protein structures
    for prot in names:
        log.info("  Analyzing %s...", prot)
        ps.align(inputpath, prot, path.join(inputpath, f"{prot}.tm"), threshold)

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
        help='Output directory into which output files are created. [./out]'
        )
    parser.add_argument(
        '--threshold', '-t', type=float, default=0.5,
        help='Polymerase+ classification threshold for TMAlign. [0.5]'
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

    main(args.inputpath, args.palmprints, args.threshold)

    
