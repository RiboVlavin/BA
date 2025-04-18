from Bio import SeqIO
import math
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import random

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
    "B": 11.5,  # B: aspartate (aspartic acid)/asparagine -> Achtung: hier wurde der Mittelwert zwischen den biden Säuren verwendet
    "C":  5.5,  # C: cystine
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
    "Z":10.55,  # Z: glutamate (glutamic acid)/glutamine-> Achtung: hier wurde der Mittelwert zwischen den beiden Säuren verwendet
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
# gab.
def calculate_cumulative_polar_requirements_difference(sequences):
    sequence_length = len(sequences[0])
    cumulative_polar_requirements_difference = 0
    for i in range(sequence_length):
        min_polar_requirements_score = 1
        max_polar_requirements_score = 0
        for sequence in sequences:
            if sequence[i] != '-':
                if get_normalized_polar_requirements(sequence[i]) < min_polar_requirements_score:
                    min_polar_requirements_score = get_normalized_polar_requirements(sequence[i])
                if get_normalized_polar_requirements(sequence[i]) > max_polar_requirements_score:
                    max_polar_requirements_score = get_normalized_polar_requirements(sequence[i])
        cumulative_polar_requirements_difference += max_polar_requirements_score - min_polar_requirements_score
    return cumulative_polar_requirements_difference/sequence_length


# returns a value between 0 and 1.
# 0 bedeutet, dass es keine Änderungen zwischen den Seuquenzen gab
# 1 bedeutet, dass es an jeder Stelle der Sequenzen zu den größten Abständen an mindestens 2 der eingegeben Sequenzen
# gab.
def calculate_cumulative_hydropathy_difference(sequences):
    sequence_length = len(sequences[0])
    cumulative_hydropathy_difference = 0
    for i in range(sequence_length):
        min_hydropathy_score = 1
        max_hydropathy_score = 0
        for sequence in sequences:
            if sequence[i] != '-':
                if get_normalized_hydropathy(sequence[i]) < min_hydropathy_score:
                    min_hydropathy_score = get_normalized_hydropathy(sequence[i])
                if get_normalized_hydropathy(sequence[i]) > max_hydropathy_score:
                    max_hydropathy_score = get_normalized_hydropathy(sequence[i])
        cumulative_hydropathy_difference += max_hydropathy_score - min_hydropathy_score
    return cumulative_hydropathy_difference/sequence_length


# returns a value between 0 and 1.
# 0 bedeutet, dass es keine Änderungen zwischen den Seuquenzen gab
# 1 bedeutet, dass es an jeder Stelle der Sequenzen zu den größten Abständen an mindestens 2 der eingegeben Sequenzen
# gab. Größter Abstand bedeutet hier, dass sowohl bei den polar_requirements als auch bei der hydropathy die jeweils
# größten und kleinsten Werte der Skala vertreten waren.
def calculate_score(sequences):
    cumulative_polar_requirements_difference = calculate_cumulative_polar_requirements_difference(sequences)
    cumulative_hydropathy_difference = calculate_cumulative_hydropathy_difference(sequences)
    return (cumulative_polar_requirements_difference + cumulative_hydropathy_difference)/2


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


def calculate_polar_requirements_scores(dataset_name, dataset_type, sequence_sections):
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        sequences.append(sequence.seq)
    score_per_section = []
    for section in sequence_sections:
        partial_sequences = []
        for sequence in sequences:
            partial_sequences.append(sequence[section[0]:section[1]])
        score_per_section.append(calculate_cumulative_polar_requirements_difference(partial_sequences))
    return score_per_section


def calculate_hydropathy_scores(dataset_name, dataset_type, sequence_sections):
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        sequences.append(sequence.seq)
    score_per_section = []
    for section in sequence_sections:
        partial_sequences = []
        for sequence in sequences:
            partial_sequences.append(sequence[section[0]:section[1]])
        score_per_section.append(calculate_cumulative_hydropathy_difference(partial_sequences))
    return score_per_section


def calculate_combined_scores(dataset_name, dataset_type, sequence_sections):
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        sequences.append(sequence.seq)
    score_per_section = []
    for section in sequence_sections:
        partial_sequences = []
        for sequence in sequences:
            partial_sequences.append(sequence[section[0]:section[1]])
        score_per_section.append(calculate_score(partial_sequences))
    return score_per_section


def calculate_mutations(dataset_name, dataset_type, sequence_sections):
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        sequences.append(sequence.seq)
    mutations_per_section = []
    for section in sequence_sections:
        partial_sequences = []
        for sequence in sequences:
            partial_sequences.append(sequence[section[0]:section[1]])
        mutations_per_section.append(calculate_mutations_per_length(partial_sequences))
    return mutations_per_section


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


def random_sequence_permutation(dataset_name, dataset_type, number_of_swaps, output_file):
    headers = []
    sequences = []
    for sequence in SeqIO.parse(dataset_name, dataset_type):
        headers.append(sequence.id)
        sequences.append(list(sequence.seq))
    for _ in range(number_of_swaps):
        index1, index2 = random.randint(0, len(sequences[0]) - 1), random.randint(0, len(sequences[0]) - 1)
        for sequence in sequences:
            seq1 = sequence[index1]
            seq2 = sequence[index2]
            sequence[index1], sequence[index2] = seq2, seq1
    modified_sequences = []
    for seq in sequences:
        modified_sequences.append("".join(seq))  # Konvertiere die Liste von Zeichen zurück in eine Sequenz
    with open(output_file, "w") as file:
        for header, seq in zip(headers, modified_sequences):
            file.write(f">{header}\n{seq}\n")


# hier sind die plots, die wir verwenden

# Berechnungen für naiven Ansatz
outer_scope_amount_of_sequences = 20
outer_scope_sequence_length = calculate_sequence_length("aligned_polyprotein_dengue_virus.fasta")
result_automatic_sequence_sections = automatic_sequence_sections(outer_scope_sequence_length, outer_scope_amount_of_sequences)
outer_scope_polar_requirements_scores = calculate_polar_requirements_scores("aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
outer_scope_hydropathy_scores = calculate_hydropathy_scores("aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
outer_scope_combined_scores = calculate_combined_scores("aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
outer_scope_mutations = calculate_mutations("aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
print("aligned_polyprotein_dengue_virus:")

# Create subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Mutationen über Sequenzen (naiver Ansatz)
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y = outer_scope_mutations
plt.plot(x, y, color='blue')
# Durchschnittswerte berechnen
mean_y = np.mean(y)
axs[0, 0].plot(x, y)
axs[0, 0].axhline(y=mean_y, color='tab:blue', linestyle='--', label='Durchschnitt')
axs[0, 0].set_xlabel('Abschnitte als Indizes')
axs[0, 0].set_ylabel('Anteil an Mutationen')
axs[0, 0].set_title('Mutationen im Dengue-Virus')
axs[0, 0].legend()
axs[0, 0].tick_params(axis='x', rotation=45)
axs[0, 0].text(-0.1, 1.05, 'A', transform=axs[0, 0].transAxes, fontsize=14, fontweight='bold')

# Plot 2: polar requirements über Sequenzen (naiver Ansatz)
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y = outer_scope_polar_requirements_scores
plt.plot(x, y, color='blue')
# Durchschnittswerte berechnen
mean_y = np.mean(y)
axs[0, 1].plot(x, y)
axs[0, 1].axhline(y=mean_y, color='tab:blue', linestyle='--', label='Durchschnitt')
axs[0, 1].set_xlabel('Abschnitte als Indizes')
axs[0, 1].set_ylabel('Durchschnittliche Polar-requirements \n Differenzen (normiert)')
axs[0, 1].set_title('Polar requirements im Dengue-Virus')
axs[0, 1].legend()
axs[0, 1].tick_params(axis='x', rotation=45)
axs[0, 1].text(-0.1, 1.05, 'B', transform=axs[0, 1].transAxes, fontsize=14, fontweight='bold')

# Plot 3: hydropathie über Sequenzen (naiver Ansatz)
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y = outer_scope_hydropathy_scores
plt.plot(x, y, color='blue')
# Durchschnittswerte berechnen
mean_y = np.mean(y)
axs[1, 0].plot(x, y)
axs[1, 0].axhline(y=mean_y, color='tab:blue', linestyle='--', label='Durchschnitt')
axs[1, 0].set_xlabel('Abschnitte als Indizes')
axs[1, 0].set_ylabel('Durchschnittliche Hydropathie \n Differenzen (normiert)')
axs[1, 0].set_title('Hydropathie im Dengue-Virus')
axs[1, 0].legend()
axs[1, 0].tick_params(axis='x', rotation=45)
axs[1, 0].text(-0.1, 1.05, 'C', transform=axs[1, 0].transAxes, fontsize=14, fontweight='bold')

# Plot 4: kombinierter Score von polar requirements und hydropathie über Sequenzen (naiver Ansatz)
axs[1, 1].clear()
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y = outer_scope_combined_scores
plt.plot(x, y, color='blue')
# Durchschnittswerte berechnen
mean_y = np.mean(y)
axs[1, 1].plot(x, y)
axs[1, 1].axhline(y=mean_y, color='tab:blue', linestyle='--', label='Durchschnitt')
axs[1, 1].set_xlabel('Abschnitte als Indizes')
axs[1, 1].set_ylabel('Mutationsscores')
axs[1, 1].set_title('Mutationsscores im Dengue-Virus')
axs[1, 1].legend()
axs[1, 1].tick_params(axis='x', rotation=45)
axs[1, 1].text(-0.1, 1.05, 'D', transform=axs[1, 1].transAxes, fontsize=14, fontweight='bold')

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()

# Print statistics for each plot
print("mutations:")
max_value = max(outer_scope_mutations)
max_index = outer_scope_mutations.index(max_value)
min_value = min(outer_scope_mutations)
min_index = outer_scope_mutations.index(min_value)
print(f"Maximum value mutations: {max_value} (index: {max_index})")
print(f"Minimum value mutations: {min_value} (index: {min_index})")
print(f"Mean mutations: {np.mean(outer_scope_mutations)}")

print("polar_requirements_scores:")
max_value = max(outer_scope_polar_requirements_scores)
max_index = outer_scope_polar_requirements_scores.index(max_value)
min_value = min(outer_scope_polar_requirements_scores)
min_index = outer_scope_polar_requirements_scores.index(min_value)
print(f"Maximum value polar_requirements_scores: {max_value} (index: {max_index})")
print(f"Minimum value polar_requirements_scores: {min_value} (index: {min_index})")
print(f"Mean polar requirements: {np.mean(outer_scope_polar_requirements_scores)}")

print("hydropathy_scores:")
max_value = max(outer_scope_hydropathy_scores)
max_index = outer_scope_hydropathy_scores.index(max_value)
min_value = min(outer_scope_hydropathy_scores)
min_index = outer_scope_hydropathy_scores.index(min_value)
print(f"Maximum value hydropathy: {max_value} (index: {max_index})")
print(f"Minimum value hydropathy: {min_value} (index: {min_index})")
print(f"Mean hydropathy: {np.mean(outer_scope_hydropathy_scores)}")

print("mutationscore:")
max_value = max(outer_scope_combined_scores)
max_index = outer_scope_combined_scores.index(max_value)
min_value = min(outer_scope_combined_scores)
min_index = outer_scope_combined_scores.index(min_value)
print(f"Maximum value mutationscore: {max_value} (index: {max_index})")
print(f"Minimum value mutationscore: {min_value} (index: {min_index})")
print(f"Mean mutationscore: {np.mean(outer_scope_combined_scores)}")

# plots mit normaler Sequenzverteilung und 20 Abschnitten (kombinierte Graphik mit Scores und Mutations; naiver Ansatz)
# Beispieldaten erstellen
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y1 = outer_scope_combined_scores
y2 = outer_scope_mutations
# Durchschnittswerte berechnen
mean_y1 = np.mean(y1)
mean_y2 = np.mean(y2)
# Maxima und Minima berechnen
max_y1, min_y1 = np.max(y1), np.min(y1)
max_y2, min_y2 = np.max(y2), np.min(y2)
# Indizes der Maxima und Minima finden
max_index_y1 = np.argmax(y1)
min_index_y1 = np.argmin(y1)
max_index_y2 = np.argmax(y2)
min_index_y2 = np.argmin(y2)
# Diagramm erstellen
fig, ax1 = plt.subplots()
# Erste y-Achse plotten
color = 'tab:red'
ax1.set_xlabel('Abschnitte als Indizes')
ax1.set_ylabel('Mutationsscore', color=color)
ax1.plot(x, y1, color=color, zorder=1)
ax1.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die erste y-Achse
ax1.axhline(mean_y1, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationsscore', zorder=1)
# Maxima und Minima für die erste y-Achse hinzufügen
ax1.scatter(x[max_index_y1], max_y1, color='red', zorder=1)  # Punkt für Maximum
ax1.scatter(x[min_index_y1], min_y1, color='red', zorder=1)  # Punkt für Minimum
# Zweite y-Achse erstellen
ax2 = ax1.twinx()  # Zweite y-Achse teilt sich die x-Achse mit der ersten
# Zweite y-Achse plotten
color = 'tab:blue'
ax2.set_ylabel('Anteil an Mutationen', color=color)
ax2.plot(x, y2, color=color, zorder=1)
ax2.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die zweite y-Achse
ax2.axhline(mean_y2, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationen', zorder=1)
# Maxima und Minima für die zweite y-Achse hinzufügen
ax2.scatter(x[max_index_y2], max_y2, color='blue', zorder=1)  # Punkt für Maximum
ax2.scatter(x[min_index_y2], min_y2, color='blue', zorder=1)  # Punkt für Minimum
# x-Achsenbeschriftungen drehen
ax1.set_xticks(ax1.get_xticks())  # Setzen der x-Achsen-Ticks explizit
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')  # Drehen der x-Achsenbeschriftungen
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.2)  # Mehr Platz unten hinzufügen
# Titel und Gitternetz hinzufügen
fig.suptitle('Mutationsscores und Mutationen im Dengue-Virus')
# Legende hinzufügen (über das gesamte Diagramm)
fig.legend(loc='upper center', bbox_to_anchor=(0.515, 0.9), ncol=2)
plt.tight_layout()
plt.show()


# plots mit random Sequenzverteilung und 20 Abschnitten (kombinierte Graphik mit Scores und Mutations; naiver Ansatz)
outer_scope_number_of_swaps = 10000
#random_sequence_permutation("aligned_polyprotein_dengue_virus.fasta", "fasta", outer_scope_number_of_swaps, "random_aligned_polyprotein_dengue_virus.fasta")
outer_scope_amount_of_sequences = 20
outer_scope_sequence_length = calculate_sequence_length("random_aligned_polyprotein_dengue_virus.fasta")
result_automatic_sequence_sections = automatic_sequence_sections(outer_scope_sequence_length, outer_scope_amount_of_sequences)
outer_scope_combined_scores = calculate_combined_scores("random_aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
outer_scope_mutations = calculate_mutations("random_aligned_polyprotein_dengue_virus.fasta", "fasta", result_automatic_sequence_sections)
print("random_aligned_polyprotein_dengue_virus:")
print(result_automatic_sequence_sections)
print(outer_scope_combined_scores)
print(outer_scope_mutations)
print("combined_scores:")
max_value = max(outer_scope_combined_scores)
max_index = outer_scope_combined_scores.index(max_value)
min_value = min(outer_scope_combined_scores)
min_index = outer_scope_combined_scores.index(min_value)
print(f"Maximum value: {max_value} (index: {max_index})")
print(f"Minimum value: {min_value} (index: {min_index})")
print("mutations:")
max_value = max(outer_scope_mutations)
max_index = outer_scope_mutations.index(max_value)
min_value = min(outer_scope_mutations)
min_index = outer_scope_mutations.index(min_value)
print(f"Maximum value: {max_value} (index: {max_index})")
print(f"Minimum value: {min_value} (index: {min_index})")
# plot erstellen mit Scores und Mutationen
# Beispieldaten erstellen
x = [f"{start}-{end}" for start, end in result_automatic_sequence_sections]
y1 = outer_scope_combined_scores
y2 = outer_scope_mutations
# Durchschnittswerte berechnen
mean_y1 = np.mean(y1)
mean_y2 = np.mean(y2)
# Maxima und Minima berechnen
max_y1, min_y1 = np.max(y1), np.min(y1)
max_y2, min_y2 = np.max(y2), np.min(y2)
# Indizes der Maxima und Minima finden
max_index_y1 = np.argmax(y1)
min_index_y1 = np.argmin(y1)
max_index_y2 = np.argmax(y2)
min_index_y2 = np.argmin(y2)
# Diagramm erstellen
fig, ax1 = plt.subplots()
# Erste y-Achse plotten
color = 'tab:red'
ax1.set_xlabel('Abschnitte als Indizes')
ax1.set_ylabel('Mutationsscore', color=color)
ax1.plot(x, y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die erste y-Achse
ax1.axhline(mean_y1, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationsscore')
# Maxima und Minima für die erste y-Achse hinzufügen
ax1.scatter(x[max_index_y1], max_y1, color='red', zorder=5)  # Punkt für Maximum
ax1.scatter(x[min_index_y1], min_y1, color='red', zorder=5)  # Punkt für Minimum
# Zweite y-Achse erstellen
ax2 = ax1.twinx()  # Zweite y-Achse teilt sich die x-Achse mit der ersten
# Zweite y-Achse plotten
color = 'tab:blue'
ax2.set_ylabel('Anteil an Mutationen', color=color)
ax2.plot(x, y2, color=color)
ax2.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die zweite y-Achse
ax2.axhline(mean_y2, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationen')
# Maxima und Minima für die zweite y-Achse hinzufügen
ax2.scatter(x[max_index_y2], max_y2, color='blue', zorder=5)  # Punkt für Maximum
ax2.scatter(x[min_index_y2], min_y2, color='blue', zorder=5)  # Punkt für Minimum
# x-Achsenbeschriftungen drehen
ax1.set_xticks(ax1.get_xticks())  # Setzen der x-Achsen-Ticks explizit
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')  # Drehen der x-Achsenbeschriftungen
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.2)  # Mehr Platz unten hinzufügen
# Titel und Gitternetz hinzufügen
fig.suptitle('Mutationsscores und Mutationen mit zufälliger Permutation')
# Legende hinzufügen (über das gesamte Diagramm)
fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.9), ncol=2)
plt.tight_layout()
plt.show()
print(f"Mean mutationscore: {mean_y1}")
print(f"Mean mutations: {mean_y2}")


# plots mit Bereichen/ Proteinen (kombinierte Graphik mit Scores und Mutations)
outer_scope_file_names, outer_scope_score_per_file, outer_scope_mutations_per_file = generate_output_with_a_folder("aligned_dengue_virus")
print("proteins_aligned_polyprotein_dengue_virus:")
print("mutationscores:")
max_value = max(outer_scope_score_per_file)
max_index = outer_scope_score_per_file.index(max_value)
min_value = min(outer_scope_score_per_file)
min_index = outer_scope_score_per_file.index(min_value)
print(f"Maximum value mutationscores: {max_value} (index: {max_index})")
print(f"Minimum value mutationscores: {min_value} (index: {min_index})")
print("mutations:")
max_value = max(outer_scope_mutations_per_file)
max_index = outer_scope_mutations_per_file.index(max_value)
min_value = min(outer_scope_mutations_per_file)
min_index = outer_scope_mutations_per_file.index(min_value)
print(f"Maximum value mutations: {max_value} (index: {max_index})")
print(f"Minimum value mutations: {min_value} (index: {min_index})")
# plot erstellen mit Scores und Mutationen
# Beispieldaten erstellen
x = outer_scope_file_names
y1 = outer_scope_score_per_file
y2 = outer_scope_mutations_per_file
# Durchschnittswerte berechnen
mean_y1 = np.mean(y1)
mean_y2 = np.mean(y2)
# Maxima und Minima berechnen
max_y1, min_y1 = np.max(y1), np.min(y1)
max_y2, min_y2 = np.max(y2), np.min(y2)
# Indizes der Maxima und Minima finden
max_index_y1 = np.argmax(y1)
min_index_y1 = np.argmin(y1)
max_index_y2 = np.argmax(y2)
min_index_y2 = np.argmin(y2)
# Diagramm erstellen
fig, ax1 = plt.subplots()
# Erste y-Achse plotten
color = 'tab:red'
ax1.set_xlabel('Namen der Regionen/ Proteine')
ax1.set_ylabel('Mutationsscore', color=color)
ax1.plot(x, y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die erste y-Achse
ax1.axhline(mean_y1, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationsscore')
# Maxima und Minima für die erste y-Achse hinzufügen
ax1.scatter(x[max_index_y1], max_y1, color='red', zorder=5)  # Punkt für Maximum
ax1.scatter(x[min_index_y1], min_y1, color='red', zorder=5)  # Punkt für Minimum
# Zweite y-Achse erstellen
ax2 = ax1.twinx()  # Zweite y-Achse teilt sich die x-Achse mit der ersten
# Zweite y-Achse plotten
color = 'tab:blue'
ax2.set_ylabel('Anteil an Mutationen', color=color)
ax2.plot(x, y2, color=color)
ax2.tick_params(axis='y', labelcolor=color)
# Durchschnittslinie für die zweite y-Achse
ax2.axhline(mean_y2, color=color, linestyle='--', linewidth=1, label='Durchschn. Mutationen')
# Maxima und Minima für die zweite y-Achse hinzufügen
ax2.scatter(x[max_index_y2], max_y2, color='blue', zorder=5)  # Punkt für Maximum
ax2.scatter(x[min_index_y2], min_y2, color='blue', zorder=5)  # Punkt für Minimum
# x-Achsenbeschriftungen drehen
ax1.set_xticks(ax1.get_xticks())  # Setzen der x-Achsen-Ticks explizit
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')  # Drehen der x-Achsenbeschriftungen
# Mehr Platz für die x-Achsenbeschriftungen schaffen
plt.subplots_adjust(bottom=0.48)  # Mehr Platz unten hinzufügen
# Titel und Gitternetz hinzufügen
fig.suptitle('Mutationsscores und Mutationen im Dengue-Virus')
# Legende hinzufügen (über das gesamte Diagramm)
fig.legend(loc='upper center', bbox_to_anchor=(0.52, 0.9), ncol=2)
plt.tight_layout()
plt.show()
print(f"Mean mutationscore: {mean_y1}")
print(f"Mean mutations: {mean_y2}")
