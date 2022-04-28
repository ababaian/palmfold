from os import path, mkdir, listdir, rename, remove
from sys import stderr

import argparse
import subprocess


class PalmStructs:

    def __init__(self, folder):
        if not path.exists(folder):
            print(f"No palmprint directory found at {folder}", file=stderr)
            exit(1)
        self.folder = folder
        
        # List rdrps
        self.rdrps = []
        rdrp_filelist = path.join(folder, "rdrp.model.list")
        if not path.exists(rdrp_filelist):
            print(f"No rdrp list found", file=stderr)
            exit(1)
        with open(rdrp_filelist) as rdrp_fl:
            # Verify each rdrp file tree
            for line in rdrp_fl:
                line = line.strip()
                if len(line) > 0:
                    name = line
                    # verify fasta
                    if not path.exists(path.join(folder, "fa", f"{name}.fa")):
                        print(f"Absent fasta for molecule {name}", file=stderr)
                        continue
                    # verify full sequence
                    if not path.exists(path.join(folder, "full_length", f"{name}.pdb.gz")):
                        print(f"Absent full length sequence for molecule {name}", file=stderr)
                    # verify pdb struct
                    if not path.exists(path.join(folder, "palmprint", f"{name}.pdb")):
                        print(f"Absent palmprint for molecule {name}", file=stderr)
                        continue
                    self.rdrps.append(name)
        # list xdxp
        self.xdxps = []
        xdxp_filelist = path.join(args.palmprints, "xdxp.model.list")
        if not path.exists(xdxp_filelist):
            print(f"No rdrp list found", file=stderr)
            exit(1)
        with open(xdxp_filelist) as xdxp_fl:
            # Verify each xdxp file tree
            for line in xdxp_fl:
                line = line.strip()
                if len(line) > 0:
                    name = line
                    # verify fasta
                    if not path.exists(path.join(folder, "fa", f"{name}.fa")):
                        print(f"Absent fasta for molecule {name}", file=stderr)
                        continue
                    # verify full sequence
                    if not path.exists(path.join(folder, "full_length", f"{name}.pdb.gz")):
                        print(f"Absent full length sequence for molecule {name}", file=stderr)
                    # verify pdb struct
                    if not path.exists(path.join(folder, "palmprint", f"{name}.pdb")):
                        print(f"Absent palmprint for molecule {name}", file=stderr)
                        continue
                    self.xdxps.append(name)

        self.all_domains = self.rdrps + self.xdxps


    def align(self, directory, name, out_tsv, rdrp_threshold):
        pdb_file = [file for file in listdir(directory) if file.startswith(name) and file.endswith(".pdb")]
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
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.folder, 'palmprint', f'{domain}.pdb')}"
                ret_val = subprocess.run(tmalign_cmd.split(" "), stdout=subprocess.PIPE)
                # Parse the output
                values = ret_val.stdout.decode('utf-8').split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
                # Save best scores
                score = values[3]
                if score >= max_rdrp_score:
                    max_rdrp_score = score
                    max_rdrp = domain

            # Serch for the max score with xdxp proteins
            for domain in self.xdxps:
                # Run TM-Align
                tmalign_cmd = f"TMalign -outfmt 2 {pdb_file} {path.join(self.folder, 'palmprint', f'{domain}.pdb')}"
                ret_val = subprocess.run(tmalign_cmd.split(" "), stdout=subprocess.PIPE)
                # Parse the output
                values = ret_val.stdout.decode('utf-8').split("\n")[1]
                values = [name, domain] + [float(x) for x in values.split("\t")[2:]]
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
                realign_cmd = f"TMalign -outfmt 1 {pdb_file} {path.join(self.folder, 'palmprint', f'{domain}.pdb')} -o {pdb_file}_tmpalign"
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract statistics from a colab splited directory')
    parser.add_argument('--directory', '-d', type=str, default="./pdb", help='Directory where are stored the pdb files to inspect')
    parser.add_argument('--palmprints', '-p', type=str, default="./pol", help='Directory where are stored the palmprints')
    parser.add_argument('--threshold', '-t', type=float, default=0.5, help='Cutoff threshold for TMAlign')

    args = parser.parse_args()

    # Verify paths
    if not path.exists(args.directory):
        print(f"No pdb directory found at {args.directory}", file=stderr)
        exit(1)
    # Extract protein names
    names = [filename[:-4] for filename in listdir(args.directory) if filename.endswith(".a3m")]

    # Create the Palmprint datastructure
    ps = PalmStructs(args.palmprints)

    # Align the proteins
    for prot in names:
        ps.align(args.directory, prot, path.join(args.directory, f"{prot}.tm"), args.threshold)
