from Bio import SeqIO

extreme_value_hydropathy = 0
hydropathy_scores = {
    "A":  1.8,  # A: alanine
    "B": -3.5,  # B: aspartate (aspartic acid)/asparagine
    "C":  2.5,  # C: cystine
    "D": -3.5,  # D: aspartate (aspartic acid)
    "E": -3.5,  # E: glutamate (glutamic acid)
    "F":  2.8,  # F: phenylalanine
    "G": -0.4,  # G: glycine
    "H": -3.2,  # H: histidine
    "I":  4.5,  # I: isoleucine
    "K": -3.9,  # K: lysine
    "L":  3.8,  # L: leucine
    "M":  1.9,  # M: methionine
    "N": -3.5,  # N: asparagine
    "P": -1.6,  # P: proline
    "Q": -3.5,  # Q: glutamine
    "R": -4.5,  # R: arginine
    "S": -0.8,  # S: serine
    "T": -0.7,  # T: threonine
    "U":  2.5,  # U: selenocysteine -> different than cysteine?
    "V":  4.2,  # V: valine
    "W": -0.9,  # W: tryptophan
    "Y": -1.3,  # Y: tyrosine
    "Z": -3.5,  # Z: glutamate (glutamic acid)/glutamine
    "X":    extreme_value_hydropathy,  # X: any
    "*":    extreme_value_hydropathy,  # *: translation stop
    "-":    extreme_value_hydropathy,  # -: gap of indeterminate length
}
min_value_hydropathy = min(hydropathy_scores.values())
max_value_hydropathy = max(hydropathy_scores.values())


def get_hydropathy(aminoacid):
    return hydropathy_scores[aminoacid]


# verändert den Wertebereich der hydropathy_scores von -4,5 bis 4,5 auf 0 bis 1
def get_normalized_hydropathy(aminoacid):
    return (hydropathy_scores[aminoacid] - min_value_hydropathy) / (max_value_hydropathy - min_value_hydropathy)


extreme_value_polar_requirements = 8.95 # Mittelwert der restlichen Scores
polar_requirements_scores = {
    "A":  7.0,  # A: alanine
    "B": 11.5,  # B: aspartate (aspartic acid)/asparagine -> Achtung: hier wurde der Mittelwert zwischen den beiden Säuren verwendet
    "C": 11.5,  # C: cystine
    "D": 13.0,  # D: aspartate (aspartic acid)
    "E": 12.5,  # E: glutamate (glutamic acid)
    "F":  5.0,  # F: phenylalanine
    "G":  7.9,  # G: glycine
    "H":  8.4,  # H: histidine
    "I":  4.9,  # I: isoleucine
    "K": 10.1,  # K: lysine
    "L":  4.9,  # L: leucine
    "M":  5.3,  # M: methionine
    "N": 10.0,  # N: asparagine
    "P":  6.6,  # P: proline
    "Q":  8.6,  # Q: glutamine
    "R":  9.1,  # R: arginine
    "S":  7.5,  # S: serine
    "T":  6.6,  # T: threonine
    "U": 11.5,  # U: selenocysteine -> different than cysteine?
    "V":  5.6,  # V: valine
    "W":  5.3,  # W: tryptophan
    "Y":  5.7,  # Y: tyrosine
    "Z": 10.55,  # Z: glutamate (glutamic acid)/glutamine-> Achtung: hier wurde der Mittelwert zwischen den beiden Säuren verwendet
    "X":    extreme_value_polar_requirements,  # X: any
    "*":    extreme_value_polar_requirements,  # *: translation stop
    "-":    extreme_value_polar_requirements,  # -: gap of indeterminate length
}
min_value_polar_requirements = min(polar_requirements_scores.values())
max_value_polar_requirements = max(polar_requirements_scores.values())


def get_polar_requirements(aminoacid):
    return polar_requirements_scores[aminoacid]


# verändert den Wertebereich der polar requirements von 0 bis 13 auf 0 bis 1
def get_normalized_polar_requirements(aminoacid):
    return (polar_requirements_scores[aminoacid] - min_value_polar_requirements) / (max_value_polar_requirements - min_value_polar_requirements)


def calculate_hydropathy_score(protein_sequence):
    score = 0
    for aminoacid in protein_sequence:
        score += get_normalized_hydropathy(aminoacid)
    return score/len(protein_sequence)


def calculate_polar_requirements_score(protein_sequence):
    score = 0
    for aminoacid in protein_sequence:
        score += get_normalized_polar_requirements(aminoacid)
    return score/len(protein_sequence)


for sequence in SeqIO.parse("dengue_virus.fasta", "fasta"):
    print(sequence.id)
    print(repr(sequence.seq))
    print(repr(sequence.translate().seq))
    print(calculate_hydropathy_score(sequence.translate()))
    print(calculate_polar_requirements_score(sequence.translate()))
