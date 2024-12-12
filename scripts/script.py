import qrcode
import cv2
from PIL import Image
from qrcode.console_scripts import error_correction


def decode_qr(qr_path, output_image_path):
    """
    Decodes a QR code and writes the original image back to a file.
    """
    # Read the QR code image
    qr_image = cv2.imread(qr_path)
    detector = cv2.QRCodeDetector()

    # Decode the QR code
    data, _, _ = detector.detectAndDecode(qr_image)

    if data:
        # Write the decoded data back to an image file
        with open(output_image_path, "wb") as f:
            f.write(data.encode('latin1'))  # Encoding as latin1 ensures byte accuracy
        print(f"Image decoded and saved to {output_image_path}")
    else:
        print("Failed to decode QR code.")


import base64
import qrcode

import base64
import qrcode

from PIL import Image


def decode_qr(qr_path, output_image_path):
    qr_image = cv2.imread(qr_path)
    detector = cv2.QRCodeDetector()

    data, _, _ = detector.detectAndDecode(qr_image)

    if data:
        with open(output_image_path, "wb") as f:
            f.write(data.encode('latin1'))  # Encoding as latin1 ensures byte accuracy
        print(f"Image decoded and saved to {output_image_path}")
    else:
        print("Failed to decode QR code.")


def compress_image(input_path, output_path, max_size=(30, 30)):
    img = Image.open(input_path)
    img.thumbnail(max_size)
    img.save(output_path, format="PNG")
    print(f"Image compressed and saved to {output_path}")


compress_image("pngtree-blue-bird-vector-or-color-illustration-png-image_2013004_11zon_11zon.jpg", "compressed_image.png")
with open("compressed_image.png", "rb") as image_file:
    encoded_string = base64.b64encode(image_file.read())
qr = qrcode.QRCode(version=40, error_correction=qrcode.ERROR_CORRECT_L)
qr.add_data(encoded_string)
qr.make(fit=True)
qr_image = qr.make_image(fill_color="black", back_color="white")
qr_image.save("encoded.png")

decode_qr("encoded.png", "decoded.png")
