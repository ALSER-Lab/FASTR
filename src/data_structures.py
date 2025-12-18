from dataclasses import dataclass
import numpy as np


@dataclass
class FASTQRecord:
    index: int
    header: bytes
    sequence: np.ndarray
    quality: np.ndarray
    pair_number: int = 0

def detect_pair_number(header: str) -> int:
    """
    Detect if a read is pair 1 or pair 2 from its header (for paired-end compatability).
    Common patterns in paired-end fastq files are: /1, /2, or _1, _2
    """
    # Check for /1 or /2 pattern
    if '/1' in header:
        return 1
    elif '/2' in header:
        return 2
    elif '_1' in header:
        return 1
    elif '_2' in header:
        return 2
    
    return 0
