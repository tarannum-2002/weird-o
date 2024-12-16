import base64

DNA_REVERTED_SEQUENCE_FILE = "scripts/dna.txt"

DNA_SEQUENCE_FILE = "scripts/dna_sequence.txt"

PROTEIN_FILE = "scripts/protein.txt"

table = {
    'ATA': 'I1', 'ATC': 'I2', 'ATT': 'I3', 'ATG': 'M1',
    'ACA': 'T1', 'ACC': 'T2', 'ACG': 'T3', 'ACT': 'T4',
    'AAC': 'N1', 'AAT': 'N2', 'AAA': 'K1', 'AAG': 'K2',
    'AGC': 'S1', 'AGT': 'S2', 'AGA': 'R1', 'AGG': 'R2',
    'CTA': 'L1', 'CTC': 'L2', 'CTG': 'L3', 'CTT': 'L4',
    'CCA': 'P1', 'CCC': 'P2', 'CCG': 'P3', 'CCT': 'P4',
    'CAC': 'H1', 'CAT': 'H2', 'CAA': 'Q1', 'CAG': 'Q2',
    'CGA': 'R3', 'CGC': 'R4', 'CGG': 'R5', 'CGT': 'R6',
    'GTA': 'V1', 'GTC': 'V2', 'GTG': 'V3', 'GTT': 'V4',
    'GCA': 'A1', 'GCC': 'A2', 'GCG': 'A3', 'GCT': 'A4',
    'GAC': 'D1', 'GAT': 'D2', 'GAA': 'E1', 'GAG': 'E2',
    'GGA': 'G1', 'GGC': 'G2', 'GGG': 'G3', 'GGT': 'G4',
    'TCA': 'S3', 'TCC': 'S4', 'TCG': 'S5', 'TCT': 'S6',
    'TTC': 'F1', 'TTT': 'F2', 'TTA': 'L5', 'TTG': 'L6',
    'TAC': 'Y1', 'TAT': 'Y2', 'TAA': '_1', 'TAG': '_2',
    'TGC': 'C1', 'TGT': 'C2', 'TGA': '_3', 'TGG': 'W1',
}

reversed_table = {'I1': 'ATA', 'I2': 'ATC', 'I3': 'ATT', 'M1': 'ATG', 'T1': 'ACA', 'T2': 'ACC', 'T3': 'ACG',
                  'T4': 'ACT',
                  'N1': 'AAC', 'N2': 'AAT', 'K1': 'AAA', 'K2': 'AAG', 'S1': 'AGC', 'S2': 'AGT', 'R1': 'AGA',
                  'R2': 'AGG',
                  'L1': 'CTA', 'L2': 'CTC', 'L3': 'CTG', 'L4': 'CTT', 'P1': 'CCA', 'P2': 'CCC', 'P3': 'CCG',
                  'P4': 'CCT',
                  'H1': 'CAC', 'H2': 'CAT', 'Q1': 'CAA', 'Q2': 'CAG', 'R3': 'CGA', 'R4': 'CGC', 'R5': 'CGG',
                  'R6': 'CGT',
                  'V1': 'GTA', 'V2': 'GTC', 'V3': 'GTG', 'V4': 'GTT', 'A1': 'GCA', 'A2': 'GCC', 'A3': 'GCG',
                  'A4': 'GCT',
                  'D1': 'GAC', 'D2': 'GAT', 'E1': 'GAA', 'E2': 'GAG', 'G1': 'GGA', 'G2': 'GGC', 'G3': 'GGG',
                  'G4': 'GGT',
                  'S3': 'TCA', 'S4': 'TCC', 'S5': 'TCG', 'S6': 'TCT', 'F1': 'TTC', 'F2': 'TTT', 'L5': 'TTA',
                  'L6': 'TTG',
                  'Y1': 'TAC', 'Y2': 'TAT', '_1': 'TAA', '_2': 'TAG', 'C1': 'TGC', 'C2': 'TGT', '_3': 'TGA',
                  'W1': 'TGG'}

amino_acid_acronyms = {
    'I1': 'Isoleucine', 'I2': 'Isoleucine', 'I3': 'Isoleucine',
    'M1': 'Methionine',
    'T1': 'Threonine', 'T2': 'Threonine', 'T3': 'Threonine', 'T4': 'Threonine',
    'N1': 'Asparagine', 'N2': 'Asparagine',
    'K1': 'Lysine', 'K2': 'Lysine',
    'S1': 'Serine', 'S2': 'Serine',
    'R1': 'Arginine', 'R2': 'Arginine', 'R3': 'Arginine', 'R4': 'Arginine', 'R5': 'Arginine', 'R6': 'Arginine',
    'L1': 'Leucine', 'L2': 'Leucine', 'L3': 'Leucine', 'L4': 'Leucine', 'L5': 'Leucine', 'L6': 'Leucine',
    'P1': 'Proline', 'P2': 'Proline', 'P3': 'Proline', 'P4': 'Proline',
    'H1': 'Histidine', 'H2': 'Histidine',
    'Q1': 'Glutamine', 'Q2': 'Glutamine',
    'V1': 'Valine', 'V2': 'Valine', 'V3': 'Valine', 'V4': 'Valine',
    'A1': 'Alanine', 'A2': 'Alanine', 'A3': 'Alanine', 'A4': 'Alanine',
    'D1': 'Aspartic Acid', 'D2': 'Aspartic Acid',
    'E1': 'Glutamic Acid', 'E2': 'Glutamic Acid',
    'G1': 'Glycine', 'G2': 'Glycine', 'G3': 'Glycine', 'G4': 'Glycine',
    'S3': 'Serine', 'S4': 'Serine', 'S5': 'Serine', 'S6': 'Serine',
    'F1': 'Phenylalanine', 'F2': 'Phenylalanine',
    'Y1': 'Tyrosine', 'Y2': 'Tyrosine',
    'C1': 'Cysteine', 'C2': 'Cysteine',
    '_1': 'Stop Codon', '_2': 'Stop Codon', '_3': 'Stop Codon',
    'W1': 'Tryptophan'
}



# Step 1: Convert Image to Base64
def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode("utf-8")
    return encoded_string


# Step 2: Convert Base64 to Binary
def base64_to_binary(base64_string):
    # Convert the base64 string to bytes and then to a binary string
    binary_string = ''.join(format(byte, '08b') for byte in base64.b64decode(base64_string))
    return binary_string


# Step 3: Map Binary to DNA
def binary_to_dna(binary_string):
    binary_to_dna_map = {
        "10": "A", "00": "T", "01": "C", "11": "G"
    }
    dna_sequence = ''.join([binary_to_dna_map[binary_string[i:i + 2]] for i in range(0, len(binary_string), 2)])
    return dna_sequence


# Step 4: Save DNA Sequence to a File
def save_dna_to_file(dna_sequence, filename=DNA_SEQUENCE_FILE):
    with open(filename, "w") as file:
        file.write(dna_sequence)


# Step6: Save protein from DNA
def save_protein_to_a_file(filename=DNA_SEQUENCE_FILE):
    with open(filename, "r") as file:
        dna_sequence = file.read()

    protein = ""
    remainder = len(dna_sequence) % 3  # Check if the sequence is divisible by 3

    # Pad the sequence if necessary
    if remainder != 0:
        print(f"DNA sequence length ({len(dna_sequence)}) is not divisible by 3. Padding with A.")
        dna_sequence += "A" * (3 - remainder)  # Add padding to make the length divisible by 3

    # Translate DNA to protein
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        protein += table.get(codon, "?")  # Use "?" if codon is invalid

    with open(PROTEIN_FILE, "w") as file:
        file.write(protein)
    return protein


def protein_to_dna(protein_sequence, reverse_table):
    dna_sequence = ""
    amino_acids = []
    # Iterate over the protein sequence in steps of 2
    for i in range(0, len(protein_sequence), 2):
        # Read two characters (or less if at the end of the sequence)
        protein_chunk = protein_sequence[i:i + 2]

        if protein_chunk in reverse_table:
            dna_sequence += reverse_table[protein_chunk]

        if protein_chunk in amino_acid_acronyms:
            amino_acids.append(amino_acid_acronyms[protein_chunk])

        else:
            raise ValueError(f"Unknown protein {protein_chunk} encountered in sequence.")

        # Write the result to the file if output_file is specified
    output_file = DNA_REVERTED_SEQUENCE_FILE
    with open(output_file, "w") as f:
        f.write(dna_sequence)

    with open("scripts/amino_acids.txt", "a") as f:
        for amino_acid in amino_acids:
            f.write(amino_acid + "\n")


    return dna_sequence


# Step 5: Load DNA Sequence from File
def load_dna_from_file(filename=DNA_REVERTED_SEQUENCE_FILE):
    with open(filename, "r") as file:
        dna_sequence = file.read()
    return dna_sequence


# Step 6: Map DNA back to Binary
def dna_to_binary(dna_sequence):
    dna_to_binary_map = {
        "A": "10", "T": "00", "C": "01", "G": "11"
    }
    binary_string = ''.join([dna_to_binary_map[base] for base in dna_sequence])
    return binary_string


# Step 7: Convert Binary to Base64
def binary_to_base64(binary_string):
    # Split the binary string into 8-bit chunks
    byte_array = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
    byte_data = bytearray(int(byte, 2) for byte in byte_array)
    base64_string = base64.b64encode(byte_data).decode("utf-8")
    return base64_string


# Step 8: Convert Base64 back to Image
def base64_to_image(base64_string, output_image_path):
    image_data = base64.b64decode(base64_string)
    with open(output_image_path, "wb") as image_file:
        image_file.write(image_data)


# Main Workflow
if __name__ == "__main__":
    # 1. Convert Image to Base64
    image_path = "/Users/tarannums/projects/weird-o/scripts/Screenshot.png"  # Update with your image file path
    base64_string = image_to_base64(image_path)

    # 2. Convert Base64 to Binary
    binary_string = base64_to_binary(base64_string)

    # 3. Convert Binary to DNA
    dna_sequence = binary_to_dna(binary_string)

    # 4. Save DNA sequence to a file
    save_dna_to_file(dna_sequence)

    protein = save_protein_to_a_file()

    dna = protein_to_dna(protein_sequence=protein, reverse_table=reversed_table)

    print("DNA sequence saved to dna_sequence.txt")

    # 5. Load DNA sequence from file
    dna_sequence_from_file = load_dna_from_file()

    # 6. Convert DNA back to Binary
    binary_string_reconstructed = dna_to_binary(dna_sequence_from_file)

    # 7. Convert Binary back to Base64
    base64_string_reconstructed = binary_to_base64(binary_string_reconstructed)

    # 8. Convert Base64 back to Image
    output_image_path = "scripts/reconstructed_image.jpg"  # Specify the output path for the reconstructed image
    base64_to_image(base64_string_reconstructed, output_image_path)

    print(f"Image reconstructed and saved as {output_image_path}")
