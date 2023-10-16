from .palmgrab import palmgrab
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable
import logging as log
import subprocess

TMALIGN_=Path(__file__).parent/'TMalign'
TMALIGN=TMALIGN_ if TMALIGN_.is_file() else Path('TMalign')

def check_tmalign_version():
    '''
    TODO
    warning: the given link of TMalign binary 
    doesn't have -outfmt parameter,
    
    new-link below works well:
    https://zhanggroup.org/TM-align/TMtools20190822.tar.gz
    
    '''
    return


def tm_align(pdb_file:Path,ref_file:Path):
    '''
    run TMalign between `pdb_file` and `ref_file`
    only output the splited score lines.
    
    TODO: if we only need the score,
    consider this package rather than call a subprecess:
    https://github.com/jvkersch/tmtools  
    '''
    tmalign_cmd = f"{TMALIGN} -outfmt 2 {pdb_file} {ref_file}"
    log.info("    %s", tmalign_cmd)

    # Runs TMalign and captures output
    ret_val = subprocess.run(
        tmalign_cmd.split(" "),
        capture_output=True,
        text=True,
        close_fds=False)
    # Parse the captured output for scores
    values = ret_val.stdout.split("\n")[1]
    values = [pdb_file.stem, ref_file.stem] + [float(x) for x in values.split("\t")[2:]]
    return values


class PalmStructs:
    # Initialization Routine ------------------------------
    # Ensures ./pol dir-tree is populated
    # with the correct files (palmprint pdb, pdb, fasta)
    def __init__(self, polpath:Path):
        '''
        parameters:
        -----------
        polpath:Path(str-OK)
            Ensures the library contains \
            the correct tress:
            
            `rdrp.model.list`
            
            `xdxp.model.list`
            
            `palmprint\`
            
            (optional)`fa\` & `full_length` 
            
        properties:
        ----------
        `polpath`
        
        `all_domains`
        
        `rdrps & xdxps`
        '''
        # Check if the .pol directory exists
        polpath=Path(polpath)
        if not polpath.is_dir():
            log.error("The palmprint directory %s does not exist", polpath)
            exit(1)
        self.polpath = polpath
        
        # Verify the reference RdRp models and fasta files exist
        self.rdrps = []
        rdrp_filelist = polpath / "rdrp.model.list"

        if not rdrp_filelist.is_file():
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
                    if not (polpath/"fa"/f"{name}.fa").is_file():
                        log.warning("Fasta file missing for: %s", name)
                        continue
                    # verify full sequence
                    if not (polpath/"full_length"/f"{name}.pdb.gz").is_file():
                        log.warning("Full length structure missing for: %s", name)
                        continue
                    # verify palmprint pdb structucture (required)
                    if not (polpath/"palmprint"/f"{name}.pdb").is_file():
                        log.error("<palmprint>.pdb missing for: %s", name)
                        exit(1)
                    self.rdrps.append(name)


        # Verify the reference XdXp models and fasta files exist
        self.xdxps = []
        xdxp_filelist = polpath/"xdxp.model.list"

        if not xdxp_filelist.is_file():
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
                    if not (polpath/"fa"/f"{name}.fa").is_file():
                        log.warning("Fasta file missing for: %s", name)
                        continue
                    # verify full sequence
                    if not (polpath/"full_length"/f"{name}.pdb.gz").is_file():
                        log.warning("Full length structure missing for: %s", name)
                        continue
                    # verify palmprint pdb structucture (required)
                    if not (polpath/"palmprint"/f"{name}.pdb").is_file():
                        log.error("<palmprint>.pdb missing for: %s", name)
                        exit(1)
                    self.xdxps.append(name)

        self.all_domains = self.rdrps + self.xdxps

    # Classification Routine ------------------------------
    # Runs TMalign of input against all palmprint pdb
    # Extracts scores and re-aligns if a RdRp palmprint detected
    def align(self, pdbpath:Path, outpath:Path, name:str, tm_threshold:float=0.4):
        '''
        fetch `pdbpath`/`name`.pdb file as target.
        align it to each palm in `self.polpath`/{rdrp,xdxp}.model.list.
        record the align score in `outpath`/{name}.tm
        
        if the highest score surpass given `tm_threshold`
        realign the target pdb to highest matching (`outpath`/{name}_realign.pdb) 
        write the seq of palmprint/core (`outpath`/{name}_{pp,rc}.fa)
        
        parameters:
        -----------
        pdbpath,outpath: pathlib.Path (str compatible).
                         
        name: str
            should not contains suffix.
        
        tm_threshold: float=0.4 
            only when highest score  
            will classification/realignment be executed.
        '''
        pdbpath,outpath=Path(pdbpath),Path(outpath) # for compatibility 

        pdb_file = pdbpath / f'{name}.pdb'
        if not pdb_file.is_file():
            log.warning("No <input>.pdb file found  %s", pdb_file)
            return
        
        # Define input/output files
        outf:Callable[[str], Path] = lambda x: outpath / (pdb_file.stem+x)
        
        out_tsv,out_ppfa,out_rcfa,out_rpdb = (
            outf('.tm'),outf(".pp.fa"),outf(".rc.fa"),outf("_realigned.pdb"))

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

            log.info("  Run TMalign against RdRp palmprints")
            for domain in self.rdrps:
                values = tm_align(pdb_file,self.polpath/'palmprint'/f'{domain}.pdb')
                print("\t".join([str(v) for v in values]), file=out)
                score= values[3]
                if score >= max_rdrp_score:
                    max_rdrp_score = score
                    max_rdrp = domain

            # Run XdXp-palmprint TMalign and identify max_XdXp_Score
            log.info("  Run TMalign against XdXp palmprints")
            for domain in self.xdxps:
                # Run TM-Align
                values = tm_align(pdb_file,self.polpath/'palmprint'/f'{domain}.pdb')
                print("\t".join([str(v) for v in values]), file=out)
                score= values[3]
                if score >= max_xdxp_score:
                    max_xdxp_score = score
                    max_xdxp = domain

        # Does maxRdRP_score pass CUTOFF value
        #TODO write top matches in a log file. (now only in stdout in verbose mode)
        log.info("")
        log.info("  Max RdRp palmprint: %s - %s", max_rdrp, max_rdrp_score)
        log.info("  Max XdXp palmprint: %s - %s", max_xdxp, max_xdxp_score)
        log.info("")

        # RdRp classification
        # RdRp TMscore above re-alignment threshold
        if max(max_rdrp_score, max_xdxp_score) >= tm_threshold:
            # MAIN TODO output a pse file with:
            # 1.original pdb & best match palm pdb, aligned.
            # 2.highlighted by colors: motif 1-5
            # 3.highlighted by licorice: ion-clamps, aromatic/polar gates.
            
            log.info("  + Polymerase re-alignment threshold score reached: %s", tm_threshold)
            if max_rdrp_score > max_xdxp_score:
                # RdRp+ classification
                log.info("  ++ RdRp-classification for %s", name)
                ## fix a bug here
                domain=max_rdrp
            else:
                log.info("  %s classified as non-RdRp polymerase", name)
                domain=max_xdxp
            log.info("")
            log.info("  Re-aligning:")

            # Run TMalign to re-align input structure to top-RdRp palmprint

            realign_cmd = (f"{TMALIGN} -outfmt 1 {pdb_file} "
                           f"{self.polpath/'palmprint'/f'{domain}.pdb'} -o {out_rpdb.with_suffix('')}")
            log.info("  %s", realign_cmd)
                
            with TemporaryDirectory() as tdir:
                temp_file = Path(tdir)/f'{name}.pdb'
                log.info('execute `palmgrab` to sub-select palmprint/core fasta seq')
                with open(temp_file,'w+') as f:
                    subprocess.run(realign_cmd.split(" "), stdout=f)
                palmgrab(temp_file,f"{out_ppfa}",f"{out_rcfa}")
                
            for item in outpath.iterdir():
                if '_realigned' in item.stem and item.suffix!='.pdb':
                    item.unlink()
            log.info("  done")
            log.info("")

