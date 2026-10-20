# Revisiting AES Transciphering in BGV with Hybrid Byte-Bit Packing

This repository contains the implementation accompanying our research on efficient AES transciphering over the BGV fully homomorphic encryption scheme.

The main idea is to use a hybrid byte-bit representation during homomorphic AES evaluation. Multiplicative inversion in the AES S-box is evaluated in a byte-oriented representation, while the affine transformation and subsequent linear layers are evaluated in a bit-oriented representation over GF(2).

The implementation is based on HElib and is intended to support the experimental evaluation presented in our paper.

## Overview

Existing BGV-based AES transciphering implementations typically retain a byte-oriented ciphertext representation throughout the AES evaluation.

Our implementation introduces an intra-round representation transition:

```text
Byte-oriented representation
        |
        | Multiplicative inversion
        v
Byte-to-bit unpacking
        |
        | Affine transformation
        | ShiftRows
        | MixColumns
        | AddRoundKey
        v
Bit-oriented representation
        |
        | Bit-to-byte packing
        v
Byte-oriented representation


## Repository Structure

```text
.
├── CMakeLists.txt
├── homAES.cpp
├── homAES.h
├── combined_matrix.h
├── decrypt_matrix_v2.h
├── README.md
└── .gitignore
