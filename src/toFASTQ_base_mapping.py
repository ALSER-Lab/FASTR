import numpy as np


def create_base_map(gray_N=0, gray_A=3, gray_C=66, gray_G=129, gray_T=192):
    """
    Create lookup table for subtracting base grayscale values during quality reconstruction.
    Returns: np.ndarray of shape (256,) with int32 values
    """
    base_table = np.zeros(256, dtype=np.int32)
    base_table[gray_N:gray_A] = gray_N
    base_table[gray_A:gray_C] = gray_A
    base_table[gray_C:gray_G] = gray_C
    base_table[gray_G:gray_T] = gray_G
    base_table[gray_T:255] = gray_T
    return base_table


def reverse_base_map(gray_N=0, gray_A=3, gray_C=66, gray_G=129, gray_T=192):
    """
    Create lookup table mapping grayscale values to ASCII base characters (N/A/G/C/T).
    Returns: np.ndarray of shape (256,) with uint8 ASCII values
    """
    reverse_map = np.full(256, ord("N"), dtype=np.uint8)
    reverse_map[gray_N:gray_A] = ord("N")
    reverse_map[gray_A:gray_C] = ord("A")
    reverse_map[gray_C:gray_G] = ord("C")
    reverse_map[gray_G:gray_T] = ord("G")
    reverse_map[gray_T:255] = ord("T")
    return reverse_map
