import pandas as pd
import re
import sys
import os

def merge_all_results(size_csv, perf_csv, output_csv, sample_id, tag):
    # Load data
    size_df = pd.read_csv(size_csv)
    perf_df = pd.read_csv(perf_csv)
    
    # Filter out the summary files if they appear inside the size list
    size_df = size_df[~size_df['Filename'].str.contains('summary.csv', na=False)]

    # 1. Normalize Performance keys
    perf_df['key'] = perf_df['Filename'].str.replace(
        '-stats-time-memory.txt', '', regex=False
    )

    # Precompile regex for performance & safety
    sample_prefix = re.compile(rf'^{re.escape(sample_id)}[.-]')

    # 2. Map Size filenames to Performance keys
    def get_matching_key(fname):
        # Strip the sample ID prefix
        base = sample_prefix.sub('', fname)

        core_name = base
        ext = ""
        for e in ['.fastq', '.fastr.xz', '.fastr', '.xz', '.gz',
                  '.SAM', '.BAM', '.CRAM']:
            if base.endswith(e):
                ext = e
                core_name = base[:-len(e)]
                break

        # Specific mappings for minimap2 results
        if core_name == 'BAM_decom':
            return f'FASTR-minimap2-original-{tag}-BAM_decom'
        if core_name == 'SAM_decom':
            return f'FASTR-minimap2-original-{tag}-SAM_decom'
        if core_name == 'CRAM_decom':
            return f'FASTR-minimap2-original-{tag}-CRAM_decom'

        # Handle different output formats for the original-{tag} process
        if core_name == f'FASTR-minimap2-original-{tag}':
            if ext == '.BAM':
                return f'FASTR-minimap2-original-{tag}-BAM'
            if ext == '.CRAM':
                return f'FASTR-minimap2-original-{tag}-CRAM'
            if ext == '.SAM':
                return f'FASTR-minimap2-original-{tag}'

        return core_name

    size_df['key'] = size_df['Filename'].apply(get_matching_key)

    # 3. Outer join
    merged = pd.merge(
        perf_df,
        size_df,
        on='key',
        how='outer',
        suffixes=('_perf', '_size')
    )

    # Final table
    final_df = merged[['key', 'SizeBytes', 'Total_CPU_Sec', 'Peak_Mem_KB']]
    final_df.columns = [
        'Method_Process',
        'File_Size_Bytes',
        'CPU_Time_Sec',
        'Peak_Memory_KB'
    ]

    final_df = final_df.sort_values(by='Method_Process')
    final_df.to_csv(output_csv, index=False)

    return final_df


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage:")
        print("  python merger.py <size_csv> <perf_csv> <output_csv> <sample_id> <tag>")
        print("Example:")
        print("  python merger.py HiFi-ERR13491966-size_summary.csv "
              "HiFi-ERR13491966-performance_summary.csv "
              "HiFi-ERR13491966-results.csv "
              "ERR13491966 hifi")
    else:
        merged_table = merge_all_results(
            sys.argv[1],
            sys.argv[2],
            sys.argv[3],
            sys.argv[4],  # sample_id
            sys.argv[5],  # tag
        )
        print(merged_table)
