import logging
import re

import numpy as np

logger = logging.getLogger(__name__)


def load_quality_lookup_table(filepath: str) -> np.ndarray:
    """
    Load LUT from attached LUT file.
    Any unmapped quality scores default to 0.
    """
    arr = np.zeros(256, dtype=np.uint8)

    with open(filepath, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Line {lineno}: expected 'quality scaled_value', got: '{line}'"
                )

            q, v = int(parts[0]), int(parts[1])

            if not (1 <= q <= 94):
                raise ValueError(
                    f"Line {lineno}: quality score {q} out of range (1-94)"
                )
            if not (1 <= v <= 63):
                raise ValueError(f"Line {lineno}: scaled value {v} out of range (1-63)")

            arr[q] = v

    mapped = np.count_nonzero(arr)
    logger.info(
        f"Loaded quality lookup table from '{filepath}': "
        f"{mapped} entries mapped, range {arr[arr > 0].min() if mapped else 0}-{arr.max()}"
    )
    return arr


def create_phred_quality_map(phred_offset=33, phred_alphabet_max=41):
    """Create mapping from ASCII quality characters to numeric quality scores"""
    phred_map = np.zeros(256, dtype=np.uint8)

    clamped_low_count = 0
    clamped_high_count = 0

    for ascii_val in range(256):
        quality_score = ascii_val - phred_offset

        if quality_score < 0:
            quality_score = 0
            clamped_low_count += 1
        elif quality_score > phred_alphabet_max:
            quality_score = phred_alphabet_max
            clamped_high_count += 1

        phred_map[ascii_val] = quality_score

    valid_range = (phred_offset, phred_offset + phred_alphabet_max)
    logger.info(
        f"Created Phred quality map: offset={phred_offset}, max_quality={phred_alphabet_max}"
    )
    logger.info(
        f"Valid ASCII range for quality: {valid_range[0]}-{valid_range[1]} ('{chr(valid_range[0])}' to '{chr(valid_range[1])}')"
    )

    if clamped_low_count > 0:
        logger.info(
            f"ASCII values 0-{phred_offset-1} will be clamped to quality 0 if encountered"
        )
    if clamped_high_count > 0:
        logger.info(
            f"ASCII values {phred_offset + phred_alphabet_max + 1}-255 will be clamped to quality {phred_alphabet_max} if encountered"
        )

    return phred_map


def apply_quality_to_bases(
    base_values, quality_scores, base_map, quality_lookup_table=None
):
    if quality_lookup_table is not None:
        scale_factors = quality_lookup_table[quality_scores]
    else:
        scale_factors = np.clip(quality_scores, 0, 62).astype(np.uint8)

    n_val = base_map[ord("N")]
    n_mask = base_values == n_val

    result = base_values + scale_factors

    if np.any(n_mask):
        n_scale = (scale_factors[n_mask] * 2 // 62).astype(np.uint8)
        result[n_mask] = n_val + n_scale

    return result.astype(np.uint8)
