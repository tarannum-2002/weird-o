import numpy as np
import cv2
import matplotlib.pyplot as plt

# Step 1: Load the image
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
    magnitude_spectrum = 20 * np.log(np.abs(f_shift))  # Log scale for better visibility
    return f_shift, magnitude_spectrum

# Step 3: Display obfuscated image
def display_fourier(magnitude_spectrum):
    plt.imshow(magnitude_spectrum, cmap='gray')
    plt.title("Fourier Transform (Obfuscated)")
    plt.axis("off")
    plt.show()

# Step 4: Inverse Fourier Transform
def inverse_fourier_transform(f_shift):
    f_ishift = np.fft.ifftshift(f_shift)  # Reverse the shift
    img_reconstructed = np.fft.ifft2(f_ishift)
    img_reconstructed = np.abs(img_reconstructed)
    return img_reconstructed

# Example Usage
if __name__ == "__main__":
    # Load an image
    img_path = "pngtree-blue-bird-vector-or-color-illustration-png-image_2013004_11zon_11zon.jpg"  # Update with your image path
    image = load_image(img_path)

    # Apply Fourier Transform
    f_shift, magnitude_spectrum = apply_fourier_transform(image)

    # Display obfuscated image
    display_fourier(magnitude_spectrum)

    # Reconstruct original image
    reconstructed_image = inverse_fourier_transform(f_shift)

    # Display reconstructed image
    plt.imshow(reconstructed_image, cmap='gray')
    plt.title("Reconstructed Image")
    plt.axis("off")
    plt.show()
