from Bio import SeqIO
import math
import matplotlib.pyplot as plt
from pathlib import Path

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
    "B": 11.5,  # B: easpartate (aspartic acid)/asparagine -> Achtung: hier wurde der Mittelwert zwischen den biden Säuren verwendet
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


# Funktioniert nur, wenn alle Sequenzen in sequences die gleiche Länge haben
def calculate_mutations_per_length(sequences):
    mutations = 0
    sequence_length = len(sequences[0])
    for i in range(sequence_length):
        same = True
        for sequence in sequences[1:]:
            if sequences[0][i] != sequence[i]:
                same = False
                break
        if not same:
            mutations += 1
    return mutations/sequence_length


# returns a value between 0 and 1.
# 0 bedeutet, dass es keine Änderungen zwischen den Seuquenzen gab
# 1 bedeutet, dass es an jeder Stelle der Sequenzen zu den größten Abständen an mindestens 2 der eingegeben Sequenzen
# gab. Größter Abstand bedeutet hier, dass sowohl bei den polar_requirements als auch bei der hydropathy die jeweils
# größten und kleinsten Werte der Skala vertreten waren.
def calculate_score(sequences):
    sequence_length = len(sequences[0])
    cumulative_polar_requirements_difference = 0
    cumulative_hydropathy_difference = 0
    for i in range(sequence_length):
        min_polar_requirements_score = 1
        max_polar_requirements_score = 0
        min_hydropathy_score = 1
        max_hydropathy_score = 0
        for sequence in sequences:
            if get_normalized_polar_requirements(sequence[i]) < min_polar_requirements_score:
                min_polar_requirements_score = get_normalized_polar_requirements(sequence[i])
            if get_normalized_polar_requirements(sequence[i]) > max_polar_requirements_score:
                max_polar_requirements_score = get_normalized_polar_requirements(sequence[i])
            if get_normalized_hydropathy(sequence[i]) < min_hydropathy_score:
                min_hydropathy_score = get_normalized_hydropathy(sequence[i])
            if get_normalized_hydropathy(sequence[i]) > max_hydropathy_score:
                max_hydropathy_score = get_normalized_hydropathy(sequence[i])
        cumulative_polar_requirements_difference += max_polar_requirements_score - min_polar_requirements_score
        cumulative_hydropathy_difference += max_hydropathy_score - min_hydropathy_score
    return (cumulative_polar_requirements_difference/sequence_length + cumulative_hydropathy_difference/sequence_length)/2


# erster Abschnitt beginnt bei 0 und der letzte Abschnitt endet bei der sequence_length + 1
# die Abschnitt sind also immer inklusive des ersten Wertes zu sehen und exklusive des Letzten
def automatic_sequence_sections(sequence_length, amount_of_sequences):
    sequence_sections = []
    step_width = sequence_length / amount_of_sequences
    for i in range(1, amount_of_sequences + 1):
        sequence_sections.append((math.ceil((i-1) * step_width), math.ceil(i * step_width)))
    sequence_sections[-1] = (
        math.ceil((amount_of_sequences - 1) * step_width),
        math.ceil(amount_of_sequences * step_width) + 1
    )
    return sequence_sections


def calculate_scores_and_mutations(dataset_name, dataset_type, sequence_sections):
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        sequences.append(sequence.seq)
    score_per_section = []
    mutations_per_section = []
    for section in sequence_sections:
        partial_sequences = []
        for sequence in sequences:
            partial_sequences.append(sequence[section[0]:section[1]])
        score_per_section.append(calculate_score(partial_sequences))
        mutations_per_section.append(calculate_mutations_per_length(partial_sequences))
    return score_per_section, mutations_per_section


# eine Hilfsfunktion, die die Länge der ersten Sequenz einer fasta Datei berechnet
def calculate_sequence_length(file_path):
    with open(file_path, 'r') as file:
        # Erste Zeile lesen (Sequenz-ID)
        file.readline()
        # Sequenzdaten lesen
        sequence = ''
        for line in file:
            if line.startswith('>'):
                break  # Nächste Sequenz-ID erreicht
            sequence += line.strip()
        # Länge der Sequenz berechnen
        seq_length = len(sequence)
        return seq_length


# diese Funktion nutzt einen Ordner als Eingabe und berechnet die Scores/ Mutations und gibt die Namen zurück für die Eingabe in ein Diagramm
def generate_output_with_a_folder(output_folder):
    file_names = []
    score_per_file = []
    mutations_per_file = []
    output_dir = Path(output_folder)
    # Stelle sicher, dass der Ordner existiert
    if not output_dir.exists():
        print(f"Der Ordner {output_dir} existiert nicht.")
    else:
        # Iteriere über alle Dateien im Ordner
        for file in output_dir.iterdir():
            # Überprüfe, ob die Datei eine FASTA-Datei ist
            if file.suffix == ".fasta":
                file_names.append(file.stem)
                sequences = []
                # Lese die FASTA-Datei und speichere die Sequenzen
                with open(file, 'r') as handle:
                    for sequence in SeqIO.parse(handle, "fasta"):
                        sequences.append(sequence.seq)
                score_per_file.append(calculate_score(sequences))
                mutations_per_file.append(calculate_mutations_per_length(sequences))
    return file_names, score_per_file, mutations_per_file


# ab hier verarbeiten wir einmal den dengue_virus Datensatz mit dem naiven automatisierten Sequenzteilen und einer
# Anzahl von 10 Seuqenzteilen
# wir printen zwei Diagramme (Balkendiagramme) und ein kombiniertes Diagram mit beiden y-Werten
outer_scope_amount_of_sequences = 10
outer_scope_sequence_length = calculate_sequence_length("aligned_polyprotein_dengue_virus.fasta")
result_automatic_sequence_sections = automatic_sequence_sections(outer_scope_sequence_length, outer_scope_amount_of_sequences)
scores_and_mutations = calculate_scores_and_mutations("aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)

# plot erstellen mit den Scores der Sequenzteile
categories = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
plt.bar(categories, scores_and_mutations[0])
plt.xlabel('Abschnitte als Indizes')
plt.ylabel('Score')
plt.title('Balkendiagramm: Scores')
# x-Achsenbeschriftungen drehen
plt.xticks(rotation=45, ha='right')  # Schrift um 45 Grad drehen, horizontal ausrichten
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.2)  # Mehr Platz unten hinzufügen
plt.show()

# plot erstellen mit den Anzahl an Muationen
categories = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
plt.bar(categories, scores_and_mutations[1])
plt.xlabel('Abschnitte als Indizes')
plt.ylabel('Anteil an Mutationen')
plt.title('Balkendiagramm: Mutationen')
# x-Achsenbeschriftungen drehen
plt.xticks(rotation=45, ha='right')  # Schrift um 45 Grad drehen, horizontal ausrichten
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.2)  # Mehr Platz unten hinzufügen
plt.show()

# plot erstellen mit Scores und Mutationen
# Beispieldaten erstellen
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y1 = scores_and_mutations[0]
y2 = scores_and_mutations[1]
# Diagramm erstellen
fig, ax1 = plt.subplots()
# Erste y-Achse plotten
color = 'tab:red'
ax1.set_xlabel('Abschnitte als Indizes')
ax1.set_ylabel('Score', color=color)
ax1.plot(x, y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)
# Zweite y-Achse erstellen
ax2 = ax1.twinx()  # Zweite y-Achse teilt sich die x-Achse mit der ersten
# Zweite y-Achse plotten
color = 'tab:blue'
ax2.set_ylabel('Anteil an Mutationen', color=color)
ax2.plot(x, y2, color=color)
ax2.tick_params(axis='y', labelcolor=color)
# x-Achsenbeschriftungen drehen
ax1.set_xticks(ax1.get_xticks())  # Setzen der x-Achsen-Ticks explizit
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')  # Drehen der x-Achsenbeschriftungen
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.2)  # Mehr Platz unten hinzufügen
# Titel und Gitternetz hinzufügen
fig.suptitle('Diagramm mit Scores und Mutationen')
plt.show()


# ab hier wird der Ordner aligned_dengue_virus ausgewertet
# wir printen zwei Diagramme (Balkendiagramme) und ein kombiniertes Diagram mit beiden y-Werten
outer_scope_file_names, outer_scope_score_per_file, outer_scope_mutations_per_file = generate_output_with_a_folder("aligned_dengue_virus")

# plot erstellen mit den Scores der Sequenzteile
plt.bar(outer_scope_file_names, outer_scope_score_per_file)
plt.xlabel('Namen der Regionen/ Proteine')
plt.ylabel('Score')
plt.title('Balkendiagramm: Scores')
# x-Achsenbeschriftungen drehen
plt.xticks(rotation=45, ha='right')  # Schrift um 45 Grad drehen, horizontal ausrichten
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.5)  # Mehr Platz unten hinzufügen
plt.show()

# plot erstellen mit den Anzahl an Muationen
plt.bar(outer_scope_file_names, outer_scope_mutations_per_file)
plt.xlabel('Namen der Regionen/ Proteine')
plt.ylabel('Anteil an Mutationen')
plt.title('Balkendiagramm: Mutationen')
# x-Achsenbeschriftungen drehen
plt.xticks(rotation=45, ha='right')  # Schrift um 45 Grad drehen, horizontal ausrichten
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.5)  # Mehr Platz unten hinzufügen
plt.show()

# plot erstellen mit Scores und Mutationen
# Beispieldaten erstellen
x = outer_scope_file_names
y1 = outer_scope_score_per_file
y2 = outer_scope_mutations_per_file
# Diagramm erstellen
fig, ax1 = plt.subplots()
# Erste y-Achse plotten
color = 'tab:red'
ax1.set_xlabel('Namen der Regionen/ Proteine')
ax1.set_ylabel('Score', color=color)
ax1.plot(x, y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)
# Zweite y-Achse erstellen
ax2 = ax1.twinx()  # Zweite y-Achse teilt sich die x-Achse mit der ersten
# Zweite y-Achse plotten
color = 'tab:blue'
ax2.set_ylabel('Anteil an Mutationen', color=color)
ax2.plot(x, y2, color=color)
ax2.tick_params(axis='y', labelcolor=color)
# x-Achsenbeschriftungen drehen
ax1.set_xticks(ax1.get_xticks())  # Setzen der x-Achsen-Ticks explizit
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')  # Drehen der x-Achsenbeschriftungen
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.5)  # Mehr Platz unten hinzufügen
# Titel und Gitternetz hinzufügen
fig.suptitle('Diagramm mit Scores und Mutationen')
plt.show()
