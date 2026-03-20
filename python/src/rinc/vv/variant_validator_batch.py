'''
Find the variants and transcripts where CGD and Tfx have the same c. and write 
out variant and transcript lists for submission to variantvalidator.org

Why two separate lists? When you submit the variant list to vv you can optionally
provide a list of transcripts you want to limit the results to. We limit the 
transcripts to just those that are known to CGD and Transcript Effects. 

Created on Feb 5, 2026

@author: pleyte
'''

import argparse
import logging.config
import pandas as pd
from rinc.util.log_config import LogConfig
from rinc.util import chromosome_map

class VariantValidatorBatch(object):
    def __init__(self):
        self._logger = logging.getLogger(__name__)
    
    def get_joined_dataframes(self, cgd_nomenclature_file, tfx_nomenclature_file):
        cgd_df = pd.read_csv(cgd_nomenclature_file, dtype=str)
        self._logger.info(f"Read {cgd_df.shape[0]} rows from {cgd_nomenclature_file}")
        
        tfx_df = pd.read_csv(tfx_nomenclature_file, dtype=str)
        self._logger.info(f"Read {tfx_df.shape[0]} rows from {tfx_nomenclature_file}")
        
        join_columns = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        merged_df = pd.merge(cgd_df, tfx_df, on=join_columns, how='inner')
        self._logger.info(f"Found {merged_df.shape[0]} common variant transcripts")
        
        return merged_df
    
    def add_simple_g_dot(self, df):
        """
        Add a simple g_dot field to the dataframe 
        """
        df['vv_g_dot'] = (
            df['chromosome'].astype(str).apply(chromosome_map.get_refseq) + 
            ':g.' + 
            df['position'].astype(str) + 
            df['reference'] + 
            '>' + 
            df['alt'])
        
        return df

    def write_variants(self, variant_output_file, df):
        """
        Write the simple g dots to file
        """
        v_df = df['vv_g_dot'].drop_duplicates()
        v_df.to_csv(variant_output_file, index=False, header=False)
        self._logger.info(f"Wrote {v_df.shape[0]} variants in simple g. format to {variant_output_file}")
    
    def write_transcripts(self, transcript_output_file, df):
        """
        Write the transcripts out to file 
        """
        t_df = df['cdna_transcript'].drop_duplicates()
        t_df.to_csv(transcript_output_file, index=False, header=False)
        self._logger.info(f"Wrote {t_df.shape[0]} cDNA transcripts in simple g. format to {transcript_output_file}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Create a list of variants and transcripts for submission to variantvalidator.org')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--cgd_nomenclature", help="CGD nomenclature (csv)", required=True)
    parser.add_argument("--tfx_nomenclature", help="Tfx nomenclature (csv)", required=True)
    parser.add_argument("--variant_output", help="List of variants (csv)", required=True)
    parser.add_argument("--transcript_output", help="List of transcripts (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    vvb = VariantValidatorBatch()
    nomenclature_df = vvb.get_joined_dataframes(args.cgd_nomenclature, args.tfx_nomenclature)
    nomenclature_df = vvb.add_simple_g_dot(nomenclature_df)

    vvb.write_variants(args.variant_output, nomenclature_df)
    vvb.write_transcripts(args.transcript_output, nomenclature_df)

if __name__ == '__main__':
    main()        