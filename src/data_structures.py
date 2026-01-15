from dataclasses import dataclass
import numpy as np
import numpy.typing as npt


@dataclass
class FASTQRecord:
    index: int
    header: bytes
    sequence: npt.NDArray[np.uint8]
    quality: npt.NDArray[np.uint8]
    quality_string: bytes = b''
    original_header: bytes = b''

