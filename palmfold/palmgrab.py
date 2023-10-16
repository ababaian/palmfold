import os
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
def palmgrab(inputfa:Path, palmout:Path, rdrpout:Path) :
    inputfa=Path(inputfa)
    palmout=Path(palmout)
    rdrpout=Path(rdrpout)
    PfinList = [] # palmprint fasta
    RfinList = [] # rdrpcore  fasta
    # Resulttrunc.append(SeqRecord(Seq(splitseq), id=name))
    if inputfa.is_file():
        #for input in glob.glob(os.path.join(os.path.join(inputfa,'tmfa'),'*.fa')) :
        # Import 2 TMalign fasta output
        seqList = []
        fasta_sequences = SeqIO.parse(open(inputfa), 'fasta')
        pdb_id = inputfa.stem #os.path.basename(inputfa.split('.fa')[0])

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            # Remove comments
            sequence = str(fasta.seq).split('#')[0]
            seqList.append(sequence)

        # Determine palmprint (sequence[1]) coordinates in alignment
        find=0
        lind=0
        # start coordinate
        for i, ef in enumerate(seqList[1]) :
            if ef != '-':
                find = i
                break
        # end coordinate
        #for i, el in enumerate(seqList[1][::1].split('\n')[-1]) :
        for i, el in enumerate(reversed(seqList[1])) :
            if el != '-' and el != ' ':
                lind = i
                break

        PfinList.append(SeqRecord(Seq(seqList[0][find:len(seqList[1])-lind]),
            id=pdb_id,
            description= 'pp:' + str(find) + '-' + str(len(seqList[1])-lind)))

        #print(PfinList[0])

        # Determine rdrpcore extensions
        #truncList = seqList[0].split(seqList[0][find:len(seqList[1])-lind])
        if find >= 100 :
            if lind >= 50 :
                # Both upstream/downstream extension are available
                RfinList.append(
                    SeqRecord(Seq(seqList[0][find-100:len(seqList[1])-lind+50]),
                        id=pdb_id,
                        description= 'rc:' + str(find-100) + '-' + str(len(seqList[1])-lind+50)))
            else :
                # Upstream extension available, Downstream not available
                # use end of sequence
                RfinList.append(
                    SeqRecord(Seq(seqList[0][find-100:len(seqList[1])]),
                        id=pdb_id,
                        description= 'rc:' + str(find-100) + '-' + str(len(seqList[1])) ))
        else :
            if lind >= 50:
                # Upstream extension not available, Downstream is available
                # use start of sequence
                RfinList.append(
                    SeqRecord(Seq(seqList[0][:len(seqList[1])-lind+50]),
                        id=pdb_id,
                        description= 'rc:' + str(find) + '-' + str(len(seqList[1])-lind+50)))
            else :
                # Upstream/downstream extension not available
                # use start/end of sequence
                RfinList.append(
                    SeqRecord(Seq(seqList[0]),
                        id=pdb_id,
                        description= 'rc:' + str(find) + '-' + str(len(seqList[1]))))

        #print(RfinList[0])

    # Write output files
    SeqIO.write(PfinList, str(palmout), "fasta")
    SeqIO.write(RfinList, str(rdrpout), "fasta")

if __name__ == '__main__':
    palmgrab(sys.argv[1],sys.argv[2],sys.argv[3])
