import re
from typing import Dict, Tuple, List, Optional
from collections import Counter


# Compiled regex patterns for different sequencer types
ILLUMINA_PATTERN = re.compile(r'@([^:]+):([^:]+):([^:]+):(\d+):(\d+):(\d+):(\d+)\s+(\d+):([YN]):(\d+):(.*)')
PACBIO_HIFI_PATTERN = re.compile(r'@+([^/]+)/(\d+)/ccs(?:/(\w+))?')
PACBIO_SUBREAD_PATTERN = re.compile(r'@+([^/]+)/(\d+)/(\d+)_(\d+)')
PACBIO_CLR_PATTERN = re.compile(r'@+([^/]+)/(\d+)/(\d+)_(\d+)(?:\s+RQ=[\d.]+)?')
# Handles: @HWUSI-EAS100R:6:73:941:1973#0/1
OLD_ILLUMINA_PATTERN = re.compile(r'@([^:]+):(\d+):(\d+):(\d+):(\d+)#([A-Za-z0-9]+)/(\d+)')
SRR_PATTERN = re.compile(r'@([A-Z]+)(\d+)\.(\d+)\s+(\d+)') 


def get_delimiter_for_sequencer(sequencer_type: str) -> str:
    """Get the delimiter character used by each sequencer type"""
    if sequencer_type == 'illumina' or sequencer_type == 'old_illumina':
        return ':'
    elif sequencer_type.startswith('pacbio'):
        return '/'
    elif sequencer_type == 'ont':
        return ':'
    elif sequencer_type == 'srr':
        return '.'
    return ':'


def parse_illumina_header(header: str) -> Tuple[Dict, str, str]:
    """Parse Illumina sequencer headers"""
    match = ILLUMINA_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    groups = match.groups()
    common = {
        'instrument': groups[0],
        'run_id': groups[1],
        'flowcell': groups[2],
        'lane': groups[3],
    }
    
    # Structure template with placeholders for repeating parts
    structure = f"{groups[0]}:{groups[1]}:{groups[2]}:{groups[3]}:{{REPEATING_1}}:{{REPEATING_2}}:{{REPEATING_3}} {{REPEATING_4}}:{{REPEATING_5}}:{{REPEATING_6}}:{{REPEATING_7}}"
    
    # The unique/repeating parts that change per read
    unique_id = f"{groups[4]}:{groups[5]}:{groups[6]}:{groups[7]}:{groups[8]}:{groups[9]}:{groups[10]}"
    
    return common, unique_id, structure



def parse_pacbio_hifi_header(header: str) -> Tuple[Dict, str, str]:
    """Parse PacBio HiFi headers (CCS with optional strand info)"""
    match = PACBIO_HIFI_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    common = {'movie': match.group(1), 'read_type': 'hifi'}
    
    if match.group(3):
        structure = f"{match.group(1)}/{{REPEATING_1}}/ccs/{{REPEATING_2}}"
        unique_id = f"{match.group(2)}/{match.group(3)}"
    else:
        structure = f"{match.group(1)}/{{REPEATING_1}}/ccs"
        unique_id = match.group(2)
    
    return common, unique_id, structure


def parse_pacbio_subread_header(header: str) -> Tuple[Dict, str, str]:
    """Parse PacBio subread headers (with start/end coordinates)"""
    match = PACBIO_SUBREAD_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    common = {'movie': match.group(1), 'read_type': 'subread'}
    structure = f"{match.group(1)}/{{REPEATING_1}}/{{REPEATING_2}}_{{REPEATING_3}}"
    unique_id = f"{match.group(2)}/{match.group(3)}_{match.group(4)}"
    
    return common, unique_id, structure


def parse_pacbio_clr_header(header: str) -> Tuple[Dict, str, str]:
    """Parse PacBio CLR (Continuous Long Read) headers"""
    match = PACBIO_CLR_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    common = {'movie': match.group(1), 'read_type': 'clr'}
    structure = f"{match.group(1)}/{{REPEATING_1}}/{{REPEATING_2}}_{{REPEATING_3}}"
    unique_id = f"{match.group(2)}/{match.group(3)}_{match.group(4)}"
    
    return common, unique_id, structure


def parse_ont_header(header: str) -> Tuple[Dict, str, str]:
    """Parse Oxford Nanopore headers"""
    parts = header.strip().split()
    if not parts:
        return {}, header, ""
    
    read_id = parts[0][1:] if parts[0].startswith('@') else parts[0]
    
    # Parse key-value pairs
    kvs = {}
    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            kvs[key] = value
    
    # Common metadata
    common_keys = ['runid', 'sampleid', 'model_version_id', 'basecall_model_version_id']
    common = {k: kvs[k] for k in common_keys if k in kvs}
    
    # Build structure template
    unique_keys = ['read', 'ch', 'start_time']
    structure_parts = [read_id]
    unique_parts = []
    
    for key in unique_keys:
        if key in kvs:
            structure_parts.append(f"{key}={{REPEATING_{len(unique_parts)+1}}}")
            unique_parts.append(kvs[key])
    
    structure = ':'.join(structure_parts) if len(structure_parts) > 1 else read_id
    unique_id = ':'.join(unique_parts) if unique_parts else read_id
    
    return common, unique_id, structure


def parse_srr_header(header: str) -> Tuple[Dict, str, str]:
    """Parse SRA/SRR format headers"""
    match = SRR_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    prefix, accession, read_index, spot = match.groups()
    common = {'prefix': prefix, 'accession': accession}
    structure = f"{prefix}{accession}.{{REPEATING_1}} {{REPEATING_2}}"
    unique_id = f"{read_index}:{spot}"
    
    return common, unique_id, structure


def parse_old_illumina_header(header: str) -> Tuple[Dict, str, str]:
    """Parse pre-Casava 1.8 Illumina headers"""
    match = OLD_ILLUMINA_PATTERN.match(header)
    if not match:
        return {}, header, ""
    
    groups = match.groups()
    common = {
        'instrument': groups[0],
        'lane': groups[1],
        # Old format often lacks explicit run_ids in the header, 
        # so we group by instrument+lane
    }
    
    structure = f"{groups[0]}:{groups[1]}:{{REPEATING_1}}:{{REPEATING_2}}:{{REPEATING_3}}#{{REPEATING_4}}/{{REPEATING_5}}"
    # Unique: tile:x:y#index/read_num
    unique_id = f"{groups[2]}:{groups[3]}:{groups[4]}#{groups[5]}/{groups[6]}"
    
    return common, unique_id, structure

def detect_delimiter(headers: List[str]) -> str:
    """
    Automatically detect the most common delimiter in headers.
    Prioritizes space if key=value pairs are detected (ONT format).
    """
    # Check if headers contain key/val pairs, this means it is likely ont
    if any('=' in h for h in headers[:5]):
        return ' ' # ONT headers use space as the primary delimiter btween their x=y's
    
    # For other formats, count delimiter frequency
    delimiter_candidates = [':', '/', '.', ' ']
    delimiter_counts = Counter()
    
    for header in headers:
        for delim in delimiter_candidates:
            delimiter_counts[delim] += header.count(delim)
    
    if delimiter_counts:
        return delimiter_counts.most_common(1)[0][0]
    return ':'

def analyze_headers_for_pattern(headers: List[str], sample_size: int = 100) -> Tuple[Optional[str], Optional[str], Optional[Dict]]:
    """
    Analyze a sample of headers to detect common patterns and build structure template.
    """
    if not headers:
        return None, None, None
    
    sample = headers[:min(sample_size, len(headers))]
    delimiter = detect_delimiter(sample)
    
    split_headers = [h.split(delimiter) for h in sample] # Split all headers by delimiter
    
    secondary_delimiter = None
    if delimiter == ' ' and len(split_headers[0]) > 0 and all('.' in h[0] for h in split_headers[:min(5, len(split_headers))]):
        secondary_delimiter = '.'
        new_split_headers = []
        for header_parts in split_headers:
            first_field_parts = header_parts[0].split('.')
            new_split_headers.append(first_field_parts + header_parts[1:])
        split_headers = new_split_headers
    
    num_fields = len(split_headers[0])
    
    if not all(len(h) == num_fields for h in split_headers):
        return None, None, None
    
    field_variance = []
    for pos in range(num_fields):
        values = [h[pos] for h in split_headers]
        unique_values = set(values)
        variance = 0 if len(unique_values) == 1 else len(unique_values) / len(values)
        field_variance.append({
            'position': pos,
            'variance': variance,
            'unique_count': len(unique_values),
            'most_common': Counter(values).most_common(1)[0][0] if values else None
        })
    
    structure_parts = []
    common_metadata = {}
    repeating_counter = 1
    
    for i, field_info in enumerate(field_variance):
        if field_info['variance'] == 0:
            common_metadata[f'field_{i}'] = field_info['most_common']
            structure_parts.append(field_info['most_common'])
        else:
            
            structure_parts.append(f"{{REPEATING_{repeating_counter}}}") # REPEATING placeholder creatin
            repeating_counter += 1
    
    if secondary_delimiter:
        structure_template = structure_parts[0] + secondary_delimiter + delimiter.join(structure_parts[1:])
    else:
        structure_template = delimiter.join(structure_parts)
    
    return structure_template, delimiter, common_metadata

def extract_unique_id_from_header(header: str, structure: str, delimiter: str) -> str:
    """
    Extract unique ID from header using structure template.
    """
    # Handle SRA format w/ secondary delimiter
    if delimiter == ' ' and '.' in structure.split(' ')[0]:
        header_parts = header.split(' ')
        if len(header_parts) > 0 and '.' in header_parts[0]:
            first_parts = header_parts[0].split('.')
            all_parts = first_parts + header_parts[1:]
        else:
            all_parts = header_parts
        
        structure_temp = structure.replace('.', ' ')
        structure_parts = structure_temp.split(' ')
    else:
        # Normal single delimiter, being split both by the same delimiter
        all_parts = header.split(delimiter)
        structure_parts = structure.split(delimiter)
    
    if len(all_parts) != len(structure_parts):
        return header
    
    unique_parts = []
    for h_part, s_part in zip(all_parts, structure_parts):
        if s_part.startswith('{REPEATING_'):
            unique_parts.append(h_part)
    
    return delimiter.join(unique_parts)

def reconstruct_header_from_adaptive(structure: str, unique_id: str, delimiter: str) -> str:
    """
    Reconstruct full header from adaptive structure and unique ID.
    """
    unique_parts = unique_id.split(delimiter)
    
    if delimiter == ' ' and '.' in structure.split(' ')[0]:
        structure_temp = structure.replace('.', ' ')
        structure_parts = structure_temp.split(' ')
        
        result_parts = []
        unique_idx = 0
        for s_part in structure_parts:
            if s_part.startswith('{REPEATING_'):
                if unique_idx < len(unique_parts):
                    result_parts.append(unique_parts[unique_idx])
                    unique_idx += 1
            else:
                result_parts.append(s_part)
        
        if len(result_parts) >= 2:
            result = result_parts[0] + '.' + ' '.join(result_parts[1:])
        else:
            result = ' '.join(result_parts)
    else:
        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part, 1)
    
    return f"@{result}"

def adaptive_compress_header(header: str, structure: str, delimiter: str, 
                            common_metadata: Dict) -> Tuple[Dict, str]:
    """
    Compress header using adaptive structure.
    """
    unique_id = extract_unique_id_from_header(header, structure, delimiter)
    return common_metadata, unique_id

def compress_header(header: str, sequencer_type: str) -> Tuple[Dict, str, str]:
    """
    Compress a single header based on sequencer type.
    Returns (common_metadata_dict, unique_id, structure_template)
    """
    if sequencer_type == 'illumina':
        return parse_illumina_header(header)
    elif sequencer_type == 'old_illumina':
        return parse_old_illumina_header(header)
    elif sequencer_type == 'pacbio_hifi':
        return parse_pacbio_hifi_header(header)
    elif sequencer_type == 'pacbio_clr':
        return parse_pacbio_clr_header(header)
    elif sequencer_type == 'pacbio_subread':
        return parse_pacbio_subread_header(header)
    elif sequencer_type == 'ont':
        return parse_ont_header(header)
    elif sequencer_type == 'srr':
        return parse_srr_header(header)
    else:
        return {}, header, ""


def reconstruct_header_from_structure(structure: str, unique_id: str, sequencer_type: str, pair_number: int = 0) -> str:
    """
    Reconstruct full header from structure template and unique ID.
    """
    delimiter = get_delimiter_for_sequencer(sequencer_type)
    unique_parts = unique_id.split(delimiter)
    
    result = structure
    for i, part in enumerate(unique_parts, 1):
        placeholder = f"{{REPEATING_{i}}}"
        result = result.replace(placeholder, part)
    
    if pair_number > 0:
        result = f"{result}/{pair_number}"
    
    return f"@{result}"


def format_metadata_header(common_metadata: Dict, structure: str, sequencer_type: str) -> str:
    """Format common metadata header based on sequencer type"""
    return structure


def metadata_dict_equals(dict1: Dict, dict2: Dict) -> bool:
    """Compare two metadata dictionaries for equality"""
    if set(dict1.keys()) != set(dict2.keys()):
        return False
    for key in dict1:
        if dict1[key] != dict2[key]:
            return False
    return True