from Bio.Align.Applications import TCoffeeCommandline
from Bio import AlignIO

cline = TCoffeeCommandline(infile="unaligned.fasta", output="fasta", outfile="aligned_dengue_virus_dna.fasta")
stdout, stderr = cline()
align = AlignIO.read("aligned_dengue_virus_dna.fasta", "fasta")
