'''
Convert sql export from CGD to g. variant list
    What uses this?   
Created on Jan 25, 2026

@author: pleyte
'''
import argparse
import pandas as pd
from rinc.util import chromosome_map

def _parse_args():
    parser = argparse.ArgumentParser(description='Convert csv from sql dump to variant list')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--input", help="SQL dump (csv)", required=True)    
    parser.add_argument("--output", help="Variant list (txt)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    args = _parse_args()
    
    # Read dataframe 
    df = pd.read_csv(args.input).drop_duplicates(subset=['chromosome', 'position_start', 'reference_base', 'variant_base'])

    # Create a list of the variants formatted as "chr-pos-ref-alt"
    formatted_variants = (
        df['chromosome'].astype(str).map(chromosome_map.get_refseq) + 
        ':' +
        'g.' +
        df['position_start'].astype(str) +  
        df['reference_base'] + '>' + 
        df['variant_base']
    )

    # Write out to a file (one variant per line, no header)
    formatted_variants.to_csv(args.output, index=False, header=False)
    print(f"Wrote {df.shape[0]} rows to {args.output}")    

if __name__ == '__main__':
    main()