# Search for a given DNA sequence in human genome
# Arguments: target sequence 
# Output: frequency of DNA

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse

parser = argparse.ArgumentParser(description='DNA sequence finder')
parser.add_argument('-s', '--seq', help='DNA nucleotide sequence', required=True)
args = parser.parse_args()

# Search sub in string, and output the index of all occurances in listindex. offset is used to find overlapped occurances.
def allindices(string, sub, offset=0):
	listindex=[]
	i = string.find(sub, offset)
	while i >= 0:
        	listindex.append(i)
		i = string.find(sub, i + 1)
	return listindex

###########################################################

tar_seq = args.seq

# Load genome FASTA file
genome_fa_path = "/nfs/cancer_ref01/human/37/genome.fa"

# Open output file
outputname = "_".join([tar_seq, 'output.txt'])
output = open(outputname, "w")
headerlist = ['chrom', 'num_seq']
header = "\t".join(headerlist)
output.write(header+"\n")
listseqlen=[]
# Parsing FASTA file
for seq_record in SeqIO.parse(genome_fa_path, "fasta"):
	print(seq_record.id)
	print(len(seq_record))
	all_mid_seqindex=allindices(seq_record.seq, tar_seq, 0)
	print(len(all_mid_seqindex))
	listseqlen.append(len(all_mid_seqindex))
	vals = [seq_record.id,("%d" % len(all_mid_seqindex))]
	row_line = "\t".join(vals)
	output.write(row_line + "\n")
print(sum(listseqlen))
vals = ['Total',("%d" % sum(listseqlen))]
row_line = "\t".join(vals)
output.write(row_line + "\n")
output.close()

