import re

import numpy as np
from numba import njit


def build_formula_func(formula: str):
    """
    Create callable function from quality scaling formula string.
    Returns: Function accepting nparray, returning scaled values
    """
    cleaned = re.sub(r"^\s*f\s*\(\s*x\s*\)\s*=\s*", "", formula.strip()).replace(
        "^", "**"
    )

    safe_dict = {
        "ln": np.log,
        "log": np.log,
        "log10": np.log10,
        "exp": np.exp,
        "sqrt": np.sqrt,
        "abs": np.abs,
        "min": np.minimum,
        "max": np.maximum,
        "np": np,
        "__builtins__": {},
    }

    def formula_func(x):
        local = dict(safe_dict)
        local["x"] = (
            x.astype(np.float32) if isinstance(x, np.ndarray) else np.float32(x)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            result = eval(cleaned, local)
            if isinstance(result, np.ndarray):
                result = np.nan_to_num(result, nan=1.0, posinf=63.0, neginf=1.0)
            return result

    return formula_func


@njit
def reverse_scaling_to_quality(binary_values, subtract_table, inverse_table, max_phred):
    """
    Convert binary encoded values back to PHRED quality scores.

    WARNING: This function is JIT-compiled with @njit. Do not use Python objects,
    lists, dicts, or advanced numpy operations. Only basic numpy arrays and operations are supported

    Returns: np.ndarray of reconstructed PHRED quality scores
    """
    y = binary_values - subtract_table[binary_values]

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
        elif val > max_phred:
            val = max_phred
        result[i] = val

    return result


@njit
def build_inverse_quality_table(scaled_int_for_q, q_possible, max_range):
    """
    Build inverse lookup table for quality score reconstruction with gap-filling.
    We use linear interpolation to fill gaps, making sure the CLOSEST value is preferred...

    WARNING: This function is JIT-compiled with @njit. Do not use Python objects,
    lists, dicts, or advanced numpy operations. Only basic numpy arrays and operations are supported.

    Returns: np.ndarray inverse lookup table mapping scaled values to quality scores
    """
    inverse_table = np.full(max_range + 1, -1, dtype=np.int32)

    for i in range(q_possible.shape[0]):
        s = scaled_int_for_q[i]
        if 0 <= s <= max_range:
            inverse_table[s] = int(q_possible[i])

    # Fill gaps
    valid_indices = []
    for i in range(max_range + 1):
        if inverse_table[i] != -1:
            valid_indices.append(i)

    if len(valid_indices) > 0:
        # Fill beginning
        for i in range(valid_indices[0]):
            inverse_table[i] = inverse_table[valid_indices[0]]

        # Fill gaps
        for i in range(len(valid_indices) - 1):
            start_idx = valid_indices[i]
            end_idx = valid_indices[i + 1]
            start_val = inverse_table[start_idx]
            end_val = inverse_table[end_idx]
            gap = end_idx - start_idx
            for j in range(1, gap):
                inverse_table[start_idx + j] = int(
                    round(start_val + (end_val - start_val) * j / gap)
                )

        # Fill end
        for i in range(valid_indices[-1] + 1, max_range + 1):
            inverse_table[i] = inverse_table[valid_indices[-1]]
    else:
        for i in range(max_range + 1):
            inverse_table[i] = 0

    return inverse_table


def quality_to_ascii(quality_scores: np.ndarray, phred_offset: int = 33) -> bytes:
    """Convert numeric quality scores to ASCII string"""
    return bytes((quality_scores + phred_offset).astype(np.uint8))
