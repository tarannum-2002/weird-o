import numpy as np
import cv2

# Step 1: Load Image
def load_image(path, grayscale=True):
    if grayscale:
        image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    else:
        image = cv2.imread(path)
    return image

# Step 2: Apply Fourier Transform
def apply_fourier_transform(image):
    f_transform = np.fft.fft2(image)
    f_shift = np.fft.fftshift(f_transform)  # Shift zero frequency to the center
    magnitude_spectrum = np.abs(f_shift)  # Compute magnitude spectrum
    return f_shift, magnitude_spectrum

# Step 3: Normalize and Quantize Fourier Data
def quantize_data(data, levels=256):
    # Normalize data to fit into the range [0, levels-1]
    normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
    quantized = (normalized * (levels - 1)).astype(np.uint8)  # Convert to 0-255 range
    return quantized

# Step 4: Convert Quantized Data to Binary
def convert_to_binary(data):
    # Convert each integer to an 8-bit binary string
    binary_data = [f"{val:08b}" for val in data.flatten()]
    return binary_data

# Step 5: Map Binary Data to DNA Sequence
def binary_to_dna(binary_data):
    # Define binary to DNA mapping
    binary_to_dna_map = {"00": "A", "01": "T", "10": "C", "11": "G"}
    dna_sequence = []
    for binary in binary_data:
        dna = "".join([binary_to_dna_map[binary[i:i+2]] for i in range(0, len(binary), 2)])
        dna_sequence.append(dna)
    return "".join(dna_sequence)

# Step 6: Save DNA Sequence to a File
def save_dna_to_file(dna_sequence, filename="dna_sequence.txt"):
    with open(filename, "w") as file:
        file.write(dna_sequence)

# Main Workflow
if __name__ == "__main__":
    img_path = "pngtree-blue-bird-vector-or-color-illustration-png-image_2013004_11zon_11zon.jpg"  # Replace with your image path
    image = load_image(img_path)

    # Apply Fourier Transform
    f_shift, magnitude_spectrum = apply_fourier_transform(image)

    # Quantize the Fourier data
    quantized_data = quantize_data(magnitude_spectrum)

    # Convert quantized data to binary
    binary_data = convert_to_binary(quantized_data)

    # Map binary data to DNA
    dna_sequence = binary_to_dna(binary_data)

    # Save the DNA sequence to a file
    save_dna_to_file(dna_sequence)

    print("DNA sequence saved to dna_sequence.txt")
