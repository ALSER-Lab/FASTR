import re
from typing import Dict, Tuple


# Compiled regex patterns for different sequencer types
ILLUMINA_PATTERN = re.compile(r'@([^:]+):(\d+):([^:]+):(\d+):(\d+):(\d+):(\d+)\s+(\d+):([YN]):(\d+):(.*)')
PACBIO_CCS_PATTERN = re.compile(r'@+([^/]+)/(\d+)/ccs')
PACBIO_HIFI_PATTERN = re.compile(r'@+([^/]+)/(\d+)/ccs(?:/(\w+))?')
PACBIO_SUBREAD_PATTERN = re.compile(r'@+([^/]+)/(\d+)/(\d+)_(\d+)')
PACBIO_CLR_PATTERN = re.compile(r'@+([^/]+)/(\d+)/(\d+)_(\d+)(?:\s+RQ=[\d.]+)?')
SRR_PATTERN = re.compile(r'@([A-Z]+)(\d+)\.(\d+)\s+(\d+)')


def parse_illumina_header(header: str) -> Tuple[Dict, str]:
    """Parse Illumina sequencer headers"""
    match = ILLUMINA_PATTERN.match(header)
    if not match:
        return {}, header
    
    groups = match.groups()
    common = {
        'instrument': groups[0],
        'run_id': groups[1],
        'flowcell': groups[2],
        'lane': groups[3],
    }
    unique_id = f"{groups[4]}:{groups[5]}:{groups[6]}:{groups[7]}:{groups[8]}:{groups[9]}:{groups[10]}"
    return common, unique_id


def parse_pacbio_ccs_header(header: str) -> Tuple[Dict, str]:
    """Parse PacBio CCS (Circular Consensus Sequence) headers"""
    match = PACBIO_CCS_PATTERN.match(header)
    if not match:
        return {}, header
    common = {'movie': match.group(1), 'read_type': 'ccs'}
    unique_id = match.group(2)
    return common, unique_id


def parse_pacbio_hifi_header(header: str) -> Tuple[Dict, str]:
    """Parse PacBio HiFi headers (CCS with optional strand info)"""
    match = PACBIO_HIFI_PATTERN.match(header)
    if not match:
        return {}, header
    common = {'movie': match.group(1), 'read_type': 'hifi'}
    unique_id = match.group(2)
    if match.group(3):
        unique_id = f"{unique_id}/{match.group(3)}"
    return common, unique_id


def parse_pacbio_subread_header(header: str) -> Tuple[Dict, str]:
    """Parse PacBio subread headers (with start/end coordinates)"""
    match = PACBIO_SUBREAD_PATTERN.match(header)
    if not match:
        return {}, header
    common = {'movie': match.group(1), 'read_type': 'subread'}
    unique_id = f"{match.group(2)}/{match.group(3)}_{match.group(4)}"
    return common, unique_id


def parse_pacbio_clr_header(header: str) -> Tuple[Dict, str]:
    """Parse PacBio CLR (Continuous Long Read) headers"""
    match = PACBIO_CLR_PATTERN.match(header)
    if not match:
        return {}, header
    common = {'movie': match.group(1), 'read_type': 'clr'}
    unique_id = f"{match.group(2)}/{match.group(3)}_{match.group(4)}"
    return common, unique_id


def parse_ont_header(header: str) -> Tuple[Dict, str]:
    """Parse Oxford Nanopore headers"""
    parts = header.strip().split()
    if not parts:
        return {}, header
    
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
    
    # Extract unique fields
    unique_keys = ['read', 'ch', 'start_time']
    unique_parts = [f"{k}={kvs[k]}" for k in unique_keys if k in kvs]
    
    if unique_parts:
        unique_id = f"{read_id}:{':'.join(unique_parts)}"
    else:
        unique_id = read_id
    
    return common, unique_id


def parse_srr_header(header: str) -> Tuple[Dict, str]:
    """Parse SRA/SRR format headers"""
    match = SRR_PATTERN.match(header)
    if not match:
        return {}, header
    
    prefix, accession, read_index, spot = match.groups()
    common = {'prefix': prefix, 'accession': accession}
    unique_id = f"{read_index}:{spot}"
    return common, unique_id


def compress_header(header: str, sequencer_type: str) -> Tuple[Dict, str]:
    """
    Compress a single header based on sequencer type.
    Returns (common_metadata_dict, unique_id)
    """
    if sequencer_type == 'illumina':
        return parse_illumina_header(header)
    elif sequencer_type == 'pacbio_hifi':
        return parse_pacbio_hifi_header(header)
    elif sequencer_type == 'pacbio_clr':
        return parse_pacbio_clr_header(header)
    elif sequencer_type == 'pacbio_ccs':
        return parse_pacbio_ccs_header(header)
    elif sequencer_type == 'pacbio_subread':
        return parse_pacbio_subread_header(header)
    elif sequencer_type == 'ont':
        return parse_ont_header(header)
    elif sequencer_type == 'srr':
        return parse_srr_header(header)
    else:
        return {}, header


def format_metadata_header(common_metadata: Dict, sequencer_type: str) -> str:
    """Format common metadata header based on sequencer type"""
    if sequencer_type == 'illumina':
        return f"{common_metadata.get('instrument', '')}:{common_metadata.get('run_id', '')}:{common_metadata.get('flowcell', '')}:{common_metadata.get('lane', '')}"
    elif sequencer_type in ['pacbio', 'pacbio_ccs', 'pacbio_hifi', 'pacbio_subread', 'pacbio_clr']:
        movie = common_metadata.get('movie', '')
        read_type = common_metadata.get('read_type', '')
        if read_type:
            return f"{movie}:{read_type}"
        return movie
    elif sequencer_type == 'ont':
        parts = []
        for key in ['runid', 'sampleid', 'model_version_id', 'basecall_model_version_id']:
            if key in common_metadata:
                parts.append(f"{key}={common_metadata[key]}")
        return ':'.join(parts)
    elif sequencer_type == 'srr':
        return f"{common_metadata.get('prefix', '')}{common_metadata.get('accession', '')}"
    else:
        return ''


def metadata_dict_equals(dict1: Dict, dict2: Dict) -> bool:
    """Compare two metadata dictionaries for equality"""
    if set(dict1.keys()) != set(dict2.keys()):
        return False
    for key in dict1:
        if dict1[key] != dict2[key]:
            return False
    return True