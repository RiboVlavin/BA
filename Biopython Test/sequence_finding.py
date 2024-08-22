import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Get the directory of the current script
script_dir = Path(__file__).parent

# Define the output directory
output_dir = script_dir / "output"

# Create the output directory if it doesn't exist
output_dir.mkdir(exist_ok=True)

# Path to your local GenPept file
genpept_file = "sequence.gp"  # Replace with your actual file path


# Function to create a safe filename
def safe_filename(name):
    return "".join([c for c in name if c.isalpha() or c.isdigit() or c == ' ']).rstrip()


# Read the GenPept record from the local file
record = SeqIO.read(genpept_file, "genbank")

# Process features and save individual proteins
for feature in record.features:
    if feature.type == "Region" or feature.type == "mat_peptide":
        if "product" in feature.qualifiers:
            product = feature.qualifiers["product"][0]
            start = feature.location.start
            end = feature.location.end

            # Extract the protein sequence
            protein_seq = feature.extract(record.seq)

            # Create a SeqRecord for the protein
            protein_record = SeqRecord(protein_seq,
                                       id=f"{record.id}_{start + 1}_{end}",
                                       description=product)

            # Save the protein sequence to a FASTA file in the output directory
            filename = f"{safe_filename(product)}.fasta"
            file_path = output_dir / filename
            with open(file_path, "w") as output_handle:
                SeqIO.write(protein_record, output_handle, "fasta")

            print(f"Saved {product} to {file_path}")
            print(f"Sequence: {protein_seq}\n")

print("All protein sequences have been saved to individual FASTA files in the output directory.")