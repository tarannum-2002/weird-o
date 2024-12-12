import numpy as np
import cv2
import matplotlib.pyplot as plt

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

# Step 7: Reverse the DNA Sequence and Convert Back to Binary
def dna_to_binary(dna_sequence):
    # Define DNA to binary mapping
    dna_to_binary_map = {"A": "00", "T": "01", "C": "10", "G": "11"}
    binary_data = []
    for base in dna_sequence:
        binary_data.append(dna_to_binary_map[base])

    # Join binary data into a single string and split it into 8-bit segments (bytes)
    binary_string = "".join(binary_data)
    binary_values = [binary_string[i:i+8] for i in range(0, len(binary_string), 8)]
    return binary_values

# Step 8: Convert Binary to Quantized Data
def binary_to_quantized(binary_values, levels=256):
    # Convert binary strings to integers
    int_values = [int(b, 2) for b in binary_values]
    # Convert to original quantized data (0-255 range)
    quantized_data = np.array(int_values, dtype=np.uint8)
    return quantized_data

# Step 9: Reconstruct the Fourier Transform
def reconstruct_fourier(quantized_data, shape, levels=256):
    # Ensure the quantized_data size matches the target shape
    required_size = np.prod(shape)

    # if quantized_data.size != required_size:
    #     # Resize the data to the required shape (using truncation or interpolation)
    #     quantized_data = np.resize(quantized_data, required_size)

    # Reshape quantized data into the original magnitude spectrum shape
    magnitude_spectrum = quantized_data.reshape(shape)

    # Reverse the normalization and quantization
    magnitude_spectrum = (magnitude_spectrum / (levels - 1)) * (np.max(magnitude_spectrum) - np.min(magnitude_spectrum)) + np.min(magnitude_spectrum)

    # Reconstruct the complex Fourier transform (assuming no phase information, which is an approximation)
    f_shift = magnitude_spectrum * np.exp(1j * np.zeros_like(magnitude_spectrum))  # Assuming zero phase

    return f_shift

# Step 10: Inverse Fourier Transform
def inverse_fourier_transform(f_shift):
    # Reverse the shift
    f_ishift = np.fft.ifftshift(f_shift)
    img_reconstructed = np.fft.ifft2(f_ishift)
    img_reconstructed = np.abs(img_reconstructed)
    return img_reconstructed

# Step 11: Display the Reconstructed Image
def display_reconstructed_image(image):
    plt.imshow(image, cmap='gray')
    plt.title("Reconstructed Image from DNA Sequence")
    plt.axis("off")
    plt.show()

# Main Workflow
if __name__ == "__main__":
    # Load the image
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

    # Reverse the process

    # Step 7: Reverse the DNA Sequence and Convert Back to Binary
    dna_sequence = open("dna_sequence.txt").read()
    binary_values = dna_to_binary(dna_sequence)

    # Step 8: Convert Binary to Quantized Data
    quantized_data_reversed = binary_to_quantized(binary_values)

    # Step 9: Reconstruct the Fourier Transform (assuming the original image shape was known)
    original_shape = magnitude_spectrum.shape  # Use the shape of the magnitude spectrum
    f_shift_reconstructed = reconstruct_fourier(quantized_data_reversed, original_shape)

    # Step 10: Inverse Fourier Transform
    reconstructed_image = inverse_fourier_transform(f_shift_reconstructed)

    # Step 11: Display the Reconstructed Image
    display_reconstructed_image(reconstructed_image)
