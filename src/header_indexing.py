import logging
import os
import tempfile

import numpy as np

logger = logging.getLogger(__name__)


def build_header_index(headers_file_path: str):
    """
    Build mmap index of header file offsets for fast random access.
    Returns: (mmap_array, mmap_temp_path)
    """
    logger.info(f"Building header index for {headers_file_path}...")

    line_count = 0
    with open(headers_file_path, "rb") as hf:
        for _ in hf:
            line_count += 1

    logger.info(f"Found {line_count:,} headers, creating index...")

    temp_dir = tempfile.gettempdir()
    mmap_path = os.path.join(temp_dir, f"header_index_{os.getpid()}.mmap")

    offsets_mmap = np.memmap(mmap_path, dtype=np.int64, mode="w+", shape=(line_count,))

    with open(headers_file_path, "rb") as hf:
        offset = 0
        idx = 0
        while True:
            line = hf.readline()
            if not line:
                break
            offsets_mmap[idx] = offset
            offset = hf.tell()
            idx += 1

    offsets_mmap.flush()
    logger.info(f"Header index built: {line_count:,} headers")

    return offsets_mmap, mmap_path
