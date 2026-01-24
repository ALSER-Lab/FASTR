from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass
class FASTQRecord:
    """
    Stores FASTQ record information.
    """

    index: int
    header: bytes
    sequence: npt.NDArray[np.uint8]
    quality: npt.NDArray[np.uint8]
    quality_string: bytes = b""
    original_header: bytes = b""


@dataclass
class MetadataBlock:
    """Stores metadata for a flowcell/run at top of file"""

    structure_template: str
    sequencer_type: str
    scaling_equation: str
    start_index: int
    end_index: int
