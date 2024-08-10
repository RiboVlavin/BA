from Bio.Seq import Seq

dna = Seq('TATAGCGC')
print(dna)
print(dna.complement())
print(dna.reverse_complement())

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

for sequence in SeqIO.parse("dengue_virus_dna.fasta", "fasta"):
    print(sequence.id)
    print(repr(sequence.seq))
    print(repr(sequence.reverse_complement().seq))
    print(repr(sequence.translate().seq))
    print("length:", len(sequence))
    # count is not overlapping (does not matter with just one letter)
    print("A:", sequence.count("A"))
    print("T:", sequence.count("T"))
    print("G:", sequence.count("G"))
    print("C:", sequence.count("C"))
    print("GC%:", gc_fraction(sequence))

for sequence in SeqIO.parse("aligned_dengue_virus_dna.fasta", "fasta"):
    print(sequence.id)
    print(repr(sequence.seq))
    print("length:", len(sequence))
    print("gap:", sequence.count("-"))
