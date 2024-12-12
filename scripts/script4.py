import numpy as np
import cv2
import matplotlib.pyplot as plt

# Step 1: Load DNA Sequence from File
def load_dna_from_file(filename="dna_sequence.txt"):
    with open(filename, "r") as file:
        dna_sequence = file.read()
    return dna_sequence

# Step 2: Map DNA to Binary
def dna_to_binary(dna_sequence):
    # Define DNA to binary mapping
    dna_to_binary_map = {"A": "00", "T": "01", "C": "10", "G": "11"}
    binary_data = []

    # Split the DNA sequence into chunks of 2 bases (representing 2 binary bits)
    for base in dna_sequence:
        if base in dna_to_binary_map:
            binary_data.append(dna_to_binary_map[base])

    # Join binary data into a single string and split it into 8-bit segments (bytes)
    binary_string = "".join(binary_data)
    binary_values = [binary_string[i:i+8] for i in range(0, len(binary_string), 8)]

    return binary_values

# Step 3: Convert Binary to Quantized Data
def binary_to_quantized(binary_values, levels=256):
    # Convert binary strings to integers
    int_values = [int(b, 2) for b in binary_values]
    # Convert to original quantized data (0-255 range)
    quantized_data = np.array(int_values, dtype=np.uint8)
    return quantized_data

# Step 4: Reconstruct the Fourier Transform
def reconstruct_fourier(quantized_data, shape, levels=256):
    # Ensure the quantized_data size matches the target shape
    required_size = np.prod(shape)

    if quantized_data.size != required_size:
        # Resize the data to the required shape (using truncation or interpolation)
        quantized_data = np.resize(quantized_data, required_size)

    # Reshape quantized data into the original magnitude spectrum shape
    magnitude_spectrum = quantized_data.reshape(shape)

    # Reverse the normalization and quantization
    magnitude_spectrum = (magnitude_spectrum / (levels - 1)) * (np.max(magnitude_spectrum) - np.min(magnitude_spectrum)) + np.min(magnitude_spectrum)

    # Reconstruct the complex Fourier transform (assuming no phase information, which is an approximation)
    f_shift = magnitude_spectrum * np.exp(1j * np.zeros_like(magnitude_spectrum))  # Assuming zero phase

    return f_shift

# Step 5: Inverse Fourier Transform
def inverse_fourier_transform(f_shift):
    # Reverse the shift
    f_ishift = np.fft.ifftshift(f_shift)
    img_reconstructed = np.fft.ifft2(f_ishift)
    img_reconstructed = np.abs(img_reconstructed)
    return img_reconstructed

# Step 6: Display the Reconstructed Image
def display_reconstructed_image(image):
    plt.imshow(image, cmap='gray')
    plt.title("Reconstructed Image from DNA Sequence")
    plt.axis("off")
    plt.show()

# Main Workflow
if __name__ == "__main__":
    # Load DNA sequence
    dna_sequence = load_dna_from_file("dna_sequence.txt")

    # Step 2: Map DNA to Binary
    binary_values = dna_to_binary(dna_sequence)

    # Step 3: Convert Binary to Quantized Data
    quantized_data = binary_to_quantized(binary_values)

    # Step 4: Reconstruct the Fourier Transform (assuming the original image shape was known)
    # Example shape, replace with actual image shape used during encoding
    original_shape = (256, 256)  # Update with the actual shape of the image's magnitude spectrum
    f_shift = reconstruct_fourier(quantized_data, original_shape)

    # Step 5: Inverse Fourier Transform
    reconstructed_image = inverse_fourier_transform(f_shift)

    # Step 6: Display the Reconstructed Image
    display_reconstructed_image(reconstructed_image)
