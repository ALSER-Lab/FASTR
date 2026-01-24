from typing import Dict, Tuple


def get_delimiter_for_sequencer(sequencer_type: str) -> str:
    """
    Get delimiter character for sequencer type.
    Returns: ':' for Illumina/ONT, '/' for PacBio, ' ' for SRA
    """
    if sequencer_type == "illumina" or sequencer_type == "old_illumina":
        return ":"
    elif sequencer_type.startswith("pacbio") and not sequencer_type.endswith("_sra"):
        return "/"
    elif sequencer_type == "ont":
        return ":"
    elif sequencer_type == "sra" or sequencer_type.endswith("_sra"):
        return " "
    return ":"


def parse_ont_unique_id(unique_id: str):
    """
    Parse ONT unique_id containing key/value pairs separated by colons.
    Returns: (prefix, key_value_dict)
    """
    parts = unique_id.strip().split(":")

    if (
        parts and "=" not in parts[0]
    ):  # First part might not be a kvp but rather a prefix
        prefix = parts[0]
        kv_parts = parts[1:]
    else:
        prefix = ""
        kv_parts = parts

    kvs = {}
    for part in kv_parts:
        if "=" in part:
            k, v = part.split("=", 1)
            kvs[k] = v

    return prefix, kvs


def reconstruct_header_from_structure(
    structure: str, unique_id: str, sequencer_type: str, pair_number: int = 0
) -> str:
    """
    Reconstruct full FASTQ header from structure template and shortened ID.
    Formats used are: (Illumina, PacBio, ONT, SRA, adaptive).
    Returns: Full header string starting with '@'
    """
    if sequencer_type == "adaptive":
        # Detect primary delimiter from structure
        if " " in structure:
            primary_delimiter = " "
        elif "/" in structure:
            primary_delimiter = "/"
        elif ":" in structure:
            primary_delimiter = ":"
        else:
            primary_delimiter = " "

        unique_parts = unique_id.split(primary_delimiter)

        if primary_delimiter == " " and "." in structure.split(" ")[0]:
            structure_temp = structure.replace(".", " ")
            structure_parts = structure_temp.split(" ")

            result_parts = []
            unique_idx = 0
            for s_part in structure_parts:
                if s_part.startswith("{REPEATING_"):
                    if unique_idx < len(unique_parts):
                        result_parts.append(unique_parts[unique_idx])
                        unique_idx += 1
                else:
                    result_parts.append(s_part)

            if len(result_parts) >= 2:
                result = result_parts[0] + "." + " ".join(result_parts[1:])
            else:
                result = " ".join(result_parts)
        else:
            result = structure
            for i, part in enumerate(unique_parts, 1):
                placeholder = f"{{REPEATING_{i}}}"
                result = result.replace(placeholder, part, 1)

    elif sequencer_type == "ont":
        result = unique_id

        if pair_number > 0:
            result = f"{result}/{pair_number}"

        return f"@{result}"

    elif sequencer_type == "ont_sra":
        result = structure.replace("{REPEATING_1}", unique_id)

        if pair_number > 0:
            result = f"{result}/{pair_number}"

        return f"@{result}"

    elif sequencer_type == "illumina_sra":
        unique_parts_by_space = unique_id.split(" ")

        if len(unique_parts_by_space) < 2:
            result = structure.replace("{REPEATING_1}", unique_id)
        else:
            spot = unique_parts_by_space[0]
            illumina_coords = unique_parts_by_space[1]
            coord_parts = illumina_coords.split(":")

            result = structure.replace("{REPEATING_1}", spot)

            for i, coord_part in enumerate(coord_parts):
                placeholder = f"{{REPEATING_{i+2}}}"
                result = result.replace(placeholder, coord_part, 1)

            if len(unique_parts_by_space) > 2:
                extra_fields = " ".join(unique_parts_by_space[2:])
                result = result.replace("{REPEATING_5}", extra_fields, 1)

    elif sequencer_type in ["pacbio_hifi_sra", "pacbio_clr_sra"]:
        parts_by_space = unique_id.split(" ")

        if len(parts_by_space) < 2:
            # Fallback
            result = structure
            for i, part in enumerate(parts_by_space, 1):
                result = result.replace(f"{{REPEATING_{i}}}", part, 1)
        else:
            spot = parts_by_space[0]
            pacbio_path = parts_by_space[1]
            result = structure.replace("{REPEATING_1}", spot, 1)
            path_parts = pacbio_path.split("/")

            if sequencer_type == "pacbio_clr_sra" and len(path_parts) == 2:
                hole = path_parts[0]
                coords = path_parts[1]

                result = result.replace("{REPEATING_2}", hole, 1)

                if "_" in coords:
                    coord_parts = coords.split("_")
                    result = result.replace("{REPEATING_3}", coord_parts[0], 1)
                    if len(coord_parts) > 1:
                        result = result.replace("{REPEATING_4}", coord_parts[1], 1)
                else:
                    result = result.replace("{REPEATING_3}", coords, 1)

            elif sequencer_type == "pacbio_hifi_sra":
                if len(path_parts) >= 1:
                    result = result.replace("{REPEATING_2}", path_parts[0], 1)

            if len(parts_by_space) > 2:
                extra = " ".join(parts_by_space[2:])
                result = result.replace("{REPEATING_5}", extra, 1)
                result = result.replace("{REPEATING_3}", extra, 1)

    elif sequencer_type == "sra":
        unique_parts = unique_id.split(" ")
        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part, 1)

    else:  # Non-adaptive, non-SRA sequencer types
        delimiter = get_delimiter_for_sequencer(sequencer_type)

        # Special handling for PacBio formats that use underscore sub-delimiter
        if sequencer_type in ["pacbio_clr", "pacbio_subread"]:
            # Split by primary delimiter first
            parts = unique_id.split(delimiter)
            unique_parts = []
            for part in parts:
                # If part contains underscore, split it too
                if "_" in part:
                    unique_parts.extend(part.split("_"))
                else:
                    unique_parts.append(part)
        else:
            unique_parts = unique_id.split(delimiter)

        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part)

    # Add pair number if present (applies to all formats)
    if pair_number > 0:
        result = f"{result}/{pair_number}"

    return f"@{result}"
