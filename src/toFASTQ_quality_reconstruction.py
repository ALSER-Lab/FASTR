import numpy as np
from numba import njit


def build_inverse_lut(forward_lut: np.ndarray) -> np.ndarray:
    """
    Invert a forward LUT (phred -> scaled_value) into a reverse LUT
    forward_lut is indexed by phred score, values are scaled integers.
    """
    max_scaled = int(forward_lut.max())
    inverse = np.zeros(max_scaled + 1, dtype=np.uint8)
    for phred in range(len(forward_lut)):
        scaled = int(forward_lut[phred])
        if scaled > 0:
            inverse[scaled] = phred
    return inverse


@njit
def reverse_scaling_to_quality(binary_values, subtract_table, inverse_table):
    """
    Convert binary encoded values back to PHRED quality scores.
    WARNING: This function is JIT-compiled with @njit. Do not use Python objects,
    lists, dicts, or advanced numpy operations. Only basic numpy arrays and operations are supported
    Returns: np.ndarray of reconstructed PHRED quality scores
    """
    remapped = np.empty_like(binary_values)
    for i in range(binary_values.shape[0]):
        remapped[i] = 10 if binary_values[i] == 255 else binary_values[i]
    y = remapped - subtract_table[remapped]
    # Clip between 0 and 63
    y_clipped = np.empty_like(y)
    for i in range(y.shape[0]):
        val = y[i]
        if val < 0:
            val = 0
        elif val > 63:
            val = 63
        y_clipped[i] = val
    # Map through inverse table
    result = np.empty_like(y_clipped)
    for i in range(y_clipped.shape[0]):
        val = inverse_table[y_clipped[i]]
        if val < 0:
            val = 0
        result[i] = val
    return result


def quality_to_ascii(quality_scores: np.ndarray, phred_offset: int = 33) -> bytes:
    """Convert numeric quality scores to ASCII string"""
    return bytes((quality_scores + phred_offset).astype(np.uint8))
