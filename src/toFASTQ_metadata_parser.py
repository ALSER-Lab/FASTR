import logging
import re
from typing import List, Optional, Tuple

from data_structures import MetadataBlock

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def find_structure_prefix(structure_template: str) -> str:
    """
    Find constant prefix from structure template before first {REPEATING_X} placeholder.
    Returns: str prefix for header
    """
    if not structure_template:
        return ""

    # Find first occurrence of {REPEATING_}
    match = re.search(r"\{REPEATING_\d+\}", structure_template)

    if match:
        prefix = structure_template[
            : match.start()
        ]  # Return everything before the first placeholder
        if prefix.endswith(":") or prefix.endswith(
            "/"
        ):  # Remove trailing delimiter if present
            prefix = prefix[:-1]
        return prefix

    return structure_template


def parse_metadata_header(
    data: bytes, mode: int
) -> Tuple[List[MetadataBlock], int, Optional[str], Optional[int]]:
    """
    Parse FASTR metadata headers to find reconstruction information.
    Returns: (metadata_blocks, data_start_byte, sra_accession, phred_alphabet_max,
              detected_mode, length_flag, second_head_flag, safe_mode_flag)
    """
    metadata_blocks = []
    sra_accession = None
    length_flag = False
    second_head_flag = False
    phred_alphabet_from_metadata = None
    detected_mode = mode
    safe_mode_flag = False

    if mode == 1:
        # Mode 1: Only bases compressed (no header compression)
        # Parse metadata for quality reconstruction, then find where actual sequences start

        # All modes look at the metadata header
        first_header = data.find(b"\n@")
        if first_header == -1:
            first_header = data.find(b"@")
            if first_header == -1:
                first_header = 0

        # Parse metadata section before first @ header
        if first_header > 0:
            header_section = data[:first_header].decode("utf-8", errors="ignore")
            lines = [
                line.strip() for line in header_section.split("\n") if line.strip()
            ]

            logger.info(f"Parsed {len(lines)} header lines")
            for i, line in enumerate(lines[:5]):
                logger.info(f"  Line {i}: {line[:80]}")

            line_idx = 0

            # Check for SRA accession
            if lines and lines[0].startswith("#") and "=" not in lines[0]:
                sra_accession = lines[0][1:]
                logger.info(f"Found SRA accession: {sra_accession}")
                line_idx = 1

            # Process metadata lines
            current_seq_type = None
            current_structure = None
            current_qual_scale = None

            while line_idx < len(lines):
                line = lines[line_idx]

                if line.startswith("@"):
                    logger.info(f"Reached sequence headers at line {line_idx}")
                    break

                if line.startswith("#MODE="):
                    detected_mode = int(line.split("=", 1)[1].strip())
                    line_idx += 1
                    continue

                if line.startswith("#SEQ-TYPE="):
                    current_seq_type = line.split("=", 1)[1].strip()
                    line_idx += 1
                    continue

                if line.startswith("#PHRED-ALPHABET="):
                    phred_str = line.split("=", 1)[1].strip()
                    if phred_str.startswith("PHRED_"):
                        try:
                            phred_alphabet_from_metadata = (
                                int(phred_str.split("_")[1]) - 1
                            )
                            logger.info(
                                f"Found PHRED alphabet: {phred_alphabet_from_metadata}"
                            )
                        except:
                            pass
                    line_idx += 1
                    continue

                if line.startswith("#LENGTH="):
                    length_str = line.split("=", 1)[1].strip()
                    if length_str.lower() == "y":
                        length_flag = True
                        logger.info(f"Found LENGTH flag: {length_flag}")
                    line_idx += 1
                    continue

                if line.startswith("#SECOND_HEAD="):
                    second_head_str = line.split("=", 1)[1].strip()
                    if second_head_str.lower() == "y":
                        second_head_flag = True
                        logger.info(f"Found SECOND_HEAD flag: {second_head_flag}")
                    line_idx += 1
                    continue

                if line.startswith("#SAFE_MODE="):
                    safe_mode_str = line.split("=", 1)[1].strip()
                    if safe_mode_str.lower() == "y":
                        safe_mode_flag = True
                        logger.info(f"Found SAFE_MODE flag: {second_head_flag}")
                    line_idx += 1
                    continue

                if line.startswith("#STRUCTURE:") or line.startswith("#STRUCTURE="):
                    # Always split on '=' first since that's our actual delimiter
                    if line.startswith("#STRUCTURE="):
                        current_structure = line.split("=", 1)[1].strip()
                    else:
                        current_structure = line.split(":", 1)[1].strip()
                    logger.info(f"Found STRUCTURE metadata: {current_structure}")
                    line_idx += 1
                    continue

                if line.startswith("#QUAL_SCALE="):
                    current_qual_scale = line.split("=", 1)[1].strip()
                    logger.info(f"Found equation: {current_qual_scale}")

                    if current_seq_type and current_qual_scale:
                        metadata_blocks.append(
                            MetadataBlock(
                                structure_template=current_structure or "",
                                sequencer_type=current_seq_type,
                                scaling_equation=current_qual_scale,
                                start_index=0,
                                end_index=-1,
                            )
                        )

                    line_idx += 1
                    continue

                line_idx += 1

        # Find actual data start which is the first @ header
        actual_data_start = first_header if first_header > 0 else 0
        return (
            metadata_blocks,
            actual_data_start,
            sra_accession,
            phred_alphabet_from_metadata,
            detected_mode,
            length_flag,
            second_head_flag,
            safe_mode_flag,
        )

    # Mode 0, 2, and 3: Header compression enabled
    # We reserve \xff (255 in hex) for start of sequence indicator, only for mode 3 (given it doesn't have an '@' indicator)
    first_seq_marker = data.find(b"\xff") if mode == 3 else data.find(b"\n@")
    if first_seq_marker == -1:
        return (
            metadata_blocks,
            0,
            sra_accession,
            phred_alphabet_from_metadata,
            detected_mode,
            length_flag,
            second_head_flag,
            safe_mode_flag,
        )

    search_start = max(0, first_seq_marker - 1000)
    header_before_seq = data[search_start:first_seq_marker]
    last_at = header_before_seq.rfind(b"\n@")

    if last_at != -1:
        actual_data_start = search_start + last_at + 1
    else:
        if data[0:1] == b"@":
            actual_data_start = 0
        else:
            actual_data_start = first_seq_marker

    header_section = data[:actual_data_start].decode("utf-8", errors="ignore")
    lines = [line.strip() for line in header_section.split("\n") if line.strip()]

    logger.info(f"Parsed {len(lines)} header lines")
    for i, line in enumerate(lines[:5]):  # Show first 5 lines
        logger.info(f"  Line {i}: {line[:80]}")

    line_idx = 0

    # Check for SRA accession
    if lines and lines[0].startswith("#") and "=" not in lines[0]:
        sra_accession = lines[0][1:]
        logger.info(f"Found SRA accession: {sra_accession}")
        line_idx = 1

    # Process metadata lines
    current_seq_type = None
    current_structure = None
    current_qual_scale = None

    while line_idx < len(lines):
        line = lines[line_idx]

        if line.startswith("<"):
            break

        if line.startswith("#MODE="):
            detected_mode = int(line.split("=", 1)[1].strip())
            line_idx += 1
            continue

        if line.startswith("#SEQ-TYPE="):
            current_seq_type = line.split("=", 1)[1].strip()
            line_idx += 1
            continue

        if line.startswith("#PHRED-ALPHABET="):
            phred_str = line.split("=", 1)[1].strip()
            if phred_str.startswith("PHRED_"):
                try:
                    phred_alphabet_from_metadata = int(phred_str.split("_")[1]) - 1
                    logger.info(f"Found PHRED alphabet: {phred_alphabet_from_metadata}")
                except:
                    pass
            line_idx += 1
            continue

        if line.startswith("#LENGTH="):
            length_str = line.split("=", 1)[1].strip()
            if length_str.lower() == "y":
                length_flag = True
                logger.info(f"Found LENGTH flag: {length_flag}")
            line_idx += 1
            continue

        if line.startswith("#SECOND_HEAD="):
            second_head_str = line.split("=", 1)[1].strip()
            if second_head_str.lower() == "y":
                second_head_flag = True
                logger.info(f"Found SECOND_HEAD flag: {second_head_flag}")
            line_idx += 1
            continue

        if line.startswith("#SAFE_MODE="):
            safe_mode_str = line.split("=", 1)[1].strip()
            if safe_mode_str.lower() == "y":
                safe_mode_flag = True
                logger.info(f"Found SAFE_MODE flag: {second_head_flag}")
            line_idx += 1
            continue

        if line.startswith("#STRUCTURE:") or line.startswith("#STRUCTURE="):
            if line.startswith("#STRUCTURE="):
                current_structure = line.split("=", 1)[1].strip()
            else:
                current_structure = line.split(":", 1)[1].strip()
            logger.info(f"Found STRUCTURE metadata: {current_structure}")
            line_idx += 1
            continue

        if line.startswith("#QUAL_SCALE="):
            current_qual_scale = line.split("=", 1)[1].strip()
            logger.info(f"Found equation: {current_qual_scale}")
            line_idx += 1
            continue

        line_idx += 1

    if current_seq_type and current_qual_scale:
        metadata_blocks.append(
            MetadataBlock(
                structure_template=current_structure or "",
                sequencer_type=current_seq_type,
                scaling_equation=current_qual_scale,
                start_index=0,
                end_index=-1,
            )
        )

    return (
        metadata_blocks,
        actual_data_start,
        sra_accession,
        phred_alphabet_from_metadata,
        detected_mode,
        length_flag,
        second_head_flag,
        safe_mode_flag,
    )
