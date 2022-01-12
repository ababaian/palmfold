import os
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

def palmgrab(inputdir, palmout, rdrpout) :
    PfinList = []
    RfinList = []
    # Resulttrunc.append(SeqRecord(Seq(splitseq), id=name))
    for input in glob.glob(os.path.join(os.path.join(inputdir,'tmfa'),'*.fa')) :
        seqList = []
        fasta_sequences = SeqIO.parse(open(input), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqList.append(sequence)
        find=0
        lind=0
        for i, el in enumerate(seqList[1]) :
            if el != '-':
                find = i
                break
        for i, el in enumerate(seqList[1][::1].split('\n')[-1]) :
            if el != '-' and el != ' ':
                lind = i
                break
        PfinList.append(SeqRecord(Seq(seqList[0][find:len(seqList[1])-lind]),id=input.split('.fa')[0]))
        truncList = seqList[0].split(seqList[0][find:len(seqList[1])-lind])
        if len(truncList[0]) >= 100 :
            if len(truncList[1]) >= 50 :
                RfinList.append(
                  SeqRecord(Seq(truncList[0][100:] + seqList[0][find:len(seqList[1])-lind] + truncList[1][:50]),
                            id=input.split('.fa')[0]))
            else :
                RfinList.append(
                  SeqRecord(Seq(truncList[0][100:] + seqList[0][find:len(seqList[1])-lind] + truncList[1]),
                            id=input.split('.fa')[0]))
        else :
            if len(truncList[1]) >= 50:
                RfinList.append(
                    SeqRecord(Seq(truncList[0] + seqList[0][find:len(seqList[1]) - lind] + truncList[1][:50]),
                              id=input.split('.fa')[0]))
            else:
                RfinList.append(
                    SeqRecord(Seq(truncList[0] + seqList[0][find:len(seqList[1] - lind)] + truncList[1]),
                              id=input.split('.fa')[0]))
    
    writeDir = os.path.join(inputdir, 'resfa')
    os.mkdir(writeDir)
        
    palm_out = FastaIO.FastaWriter(os.path.join(writeDir, palmout), wrap=None)
    palm_out.write_file(PfinList)
    rdrp_out =  FastaIO.FastaWriter(os.path.join(writeDir, rdrpout), wrap=None)
    rdrp_out.write_file(RfinList)

if __name__ == '__main__':
    palmgrab(sys.argv[1],sys.argv[2],sys.argv[3])
