import base64
import cv2
import numpy as np
import matplotlib.pyplot as plt


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
    dna_sequence = ''.join([binary_to_dna_map[binary_string[i:i+2]] for i in range(0, len(binary_string), 2)])
    return dna_sequence

# Step 4: Save DNA Sequence to a File
def save_dna_to_file(dna_sequence, filename="dna_sequence.txt"):
    with open(filename, "w") as file:
        file.write(dna_sequence)

# Step 5: Load DNA Sequence from File
def load_dna_from_file(filename="dna_sequence.txt"):
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
    byte_array = [binary_string[i:i+8] for i in range(0, len(binary_string), 8)]
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
    image_path = "pngtree-blue-bird-vector-or-color-illustration-png-image_2013004_11zon_11zon.jpg"  # Update with your image file path
    base64_string = image_to_base64(image_path)

    # 2. Convert Base64 to Binary
    binary_string = base64_to_binary(base64_string)

    # 3. Convert Binary to DNA
    dna_sequence = binary_to_dna(binary_string)

    # 4. Save DNA sequence to a file
    save_dna_to_file(dna_sequence)

    print("DNA sequence saved to dna_sequence.txt")

    # 5. Load DNA sequence from file
    dna_sequence_from_file = load_dna_from_file()

    # 6. Convert DNA back to Binary
    binary_string_reconstructed = dna_to_binary(dna_sequence_from_file)

    # 7. Convert Binary back to Base64
    base64_string_reconstructed = binary_to_base64(binary_string_reconstructed)

    # 8. Convert Base64 back to Image
    output_image_path = "reconstructed_image.jpg"  # Specify the output path for the reconstructed image
    base64_to_image(base64_string_reconstructed, output_image_path)

    print(f"Image reconstructed and saved as {output_image_path}")

