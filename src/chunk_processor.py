import io
import traceback
import numpy as np
from typing import Tuple, Optional
from fastq_parser import parse_fastq_records_from_buffer
from fastq_writer import process_and_write_records


def process_chunk_worker(chunk_data, base_map, phred_map, compress_headers,
                        sequencer_type, paired_end, keep_bases, 
                        keep_quality, quality_scaling, custom_formula,
                        phred_alphabet_max, min_quality, BYTE_LOOKUP, binary,
                        remove_repeating_header, adaptive_structure, 
                        adaptive_delimiter, adaptive_sample_size,
                        extract_headers=False, mode=None, safe_mode=False):
    """
    Worker function that processes a single chunk in parallel.
    Parses FASTQ records and applies transformations.
    Returns: (chunk_id, processed_bytes, metadata, structure_template, record_count)
    """
    try:
        chunk_id, buffer, start_index = chunk_data
        # print(f"Worker processing chunk {chunk_id} with {len(buffer)} bytes")
        
        # Parse the records from buffer
        records, _, metadata, structure, delimiter, count = parse_fastq_records_from_buffer(
            buffer, start_index, base_map, phred_map,
            compress_headers, sequencer_type, paired_end,
            keep_bases, keep_quality,
            adaptive_structure=adaptive_structure,
            adaptive_delimiter=adaptive_delimiter,
            adaptive_sample_size=adaptive_sample_size
        )
        
        # print(f"Worker parsed {count} records from chunk {chunk_id}")
        
        # Process records and write to in-memory buffer
        output_buffer = io.BytesIO()
        headers_buffer = io.BytesIO() if extract_headers else None

        if records:
            process_and_write_records(
                records, output_buffer, base_map, quality_scaling,
                custom_formula, phred_alphabet_max, min_quality, keep_bases,
                binary, keep_quality,
                remove_repeating_header, compress_headers, BYTE_LOOKUP,
                headers_buffer=headers_buffer,
                mode=mode,
                safe_mode=safe_mode
            )
        
       #  print(f"Worker completed chunk {chunk_id}")
        headers_data = headers_buffer.getvalue() if headers_buffer else None 
        return (chunk_id, output_buffer.getvalue(), metadata, structure, delimiter, count, headers_data)
        
    except Exception as e:
        print(f"ERROR in worker processing chunk {chunk_id}: {e}")
        traceback.print_exc()
        raise


def chunk_generator(fastq_path: str, chunk_size_bytes: int):
    """
    Generator function to yield file chunks as they're read.
    Makes sure chunks end on complete record boundaries. 
    Yields: (chunk_id, buffer_data, start_index)
    """
    buffer = b''
    chunk_id = 0
    start_index = 0
    
    with open(fastq_path, 'rb', buffering=chunk_size_bytes) as infile:
        while True:
            chunk = infile.read(chunk_size_bytes)
            if not chunk and not buffer:
                break
            
            buffer += chunk
            
            # Find last complete record boundary
            last_at = buffer.rfind(b'\n@')
            if last_at == -1 or not chunk:
                process_buffer = buffer
                buffer = b''
            else:
                process_buffer = buffer[:last_at+1]
                buffer = buffer[last_at+1:]
            
            if process_buffer:
                yield (chunk_id, process_buffer, start_index)
                
                # Estimate number of records for next chunk's start index
                estimated_records = process_buffer.count(b'\n@')
                start_index += estimated_records
                chunk_id += 1
                
                if chunk_id % 10 == 0:
                    print(f"Read {chunk_id} chunks...")