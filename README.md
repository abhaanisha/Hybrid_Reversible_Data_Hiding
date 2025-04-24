# Hybrid Encrypted Reversible Data Hiding (HERDH)

## Overview
This repository contains the implementation and results of the Hybrid Encrypted Reversible Data Hiding (HERDH) method, a novel approach for secure and reversible data hiding in images.

## Method Summary
HERDH combines Digital Wavelet Transform (DWT), Singular Value Decomposition (SVD), and encryption techniques to hide confidential data within cover images while ensuring:

1. **Security**: The confidential data is encrypted before embedding
2. **Reversibility**: Both the original cover image and confidential data can be perfectly recovered
3. **Quality**: The embedded image maintains high visual similarity to the original cover image

## Technical Approach

### Hiding Phase
1. Apply DWT to the cover image (CI) to obtain four components (LL, LH, HL, HH)
2. Extract the HH component (high-frequency details)
3. Apply SVD to the HH component to get its singular values
4. Encrypt the confidential data (CDI) using a secret key
5. Apply SVD to the encrypted confidential data (ECDI)
6. Combine the singular values from steps 3 and 5
7. Reconstruct a new HH component using inverse SVD
8. Apply inverse DWT with the modified HH component to create the embedded cover image (ECI)

### Extraction Phase
1. Apply DWT to the embedded image (ECI)
2. Apply SVD to the HH component of the embedded image
3. Subtract the original cover image's HH singular values
4. Reconstruct the encrypted confidential data using inverse SVD
5. Decrypt to recover the original confidential data
6. Separately reconstruct the original cover image

## Experimental Results

The method was tested on several standard test images with different types of confidential data. Results show excellent performance in terms of:

| Image | Embedding PSNR (dB) | Recovery PSNR (dB) | Confidential Data Type |
|-------|---------------------|--------------------|-----------------------|
| Barbara | 50.69 / 51.38 | 100.00 / 100.00 | Diagonal pattern / Color gradient |
| Boat | 50.64 / 51.31 | 86.52 / 90.28 | Diagonal pattern / Color gradient |
| Lena | 50.66 / 51.31 | 83.40 / 84.46 | Diagonal pattern / Color gradient |
| Pentagon | 50.63 / 51.32 | 100.00 / 100.00 | Diagonal pattern / Color gradient |

### Visual Results

#### Barbara Image Results
![Barbara Results](barbara_results.png)

#### Boat Image Results
![Boat Results](boat_results.png)

#### Lena Image Results
![Lena Results](lena_results.png)

#### Pentagon Image Results
![Pentagon Results](pentagon_results.png)

As shown in the images above, each test case includes:
- Original cover image
- Embedded image with PSNR values (showing minimal visual difference)
- Recovered confidential data with high PSNR values
- The confidential data (diagonal pattern in top row, color gradient in bottom row)

### Key Findings:
- The embedded images maintain high quality with PSNR values consistently above 50 dB
- Perfect or near-perfect recovery of confidential data (PSNR values of 83-100 dB)
- The method works effectively with different types of cover images and confidential data
- Two different confidential data types were tested: diagonal patterns and color gradients

## Advantages
- Higher embedding capacity compared to traditional methods
- Superior visual quality of embedded images
- Perfect reversibility for both cover image and confidential data
- Enhanced security through encryption

## Applications
- Medical image security
- Military communications
- Copyright protection
- Secure transmission of sensitive information
