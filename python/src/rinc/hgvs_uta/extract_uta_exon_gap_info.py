'''
Query UTA for all transcripts that have cigar strings indicating reference-refseq differences. 
Write CSV with each transcript and their cigar strings. 
Created on Jan 27, 2026

@author: pleyte
'''

import argparse
import logging.config
import os
from rinc.util.log_config import LogConfig
from rinc.util.uta_db import UtaDb
import pandas as pd
    
class ExtractUtaExonGapInfo(object):
    '''
    classdocs
    '''
    def __init__(self, uta_schema):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._uta_schema = uta_schema
        
        utaDbUrl = os.environ.get('UTA_DB_URL')
        if utaDbUrl:
            self._logger.info(f"Using UTA schema {self._uta_schema} at UTA_DB_URL={utaDbUrl}")
        else:
            self._logger.info(f"Using remote UTA db with schema {self._uta_schema}")
    
    def get_transcripts_with_gaps(self) -> dict:
        """
        Query UTA for all transcripts that have a cigar string that indicates a difference from reference genome. 
        If the whole exon matches the cigar string will look like "123=". Everything else indicates a mismatch.
        Returns a dataframe with one row for each accession, and a 'cigar' field that is a concatenated list of all cigar strings for the transcript.  
        """
        query = f"""
                SELECT * 
                  FROM {self._uta_schema}.tx_exon_aln_v 
                 WHERE alt_aln_method = 'splign'
                  -- Filter 1: Ensure we only pull hg19 rows for the final result
                  AND alt_ac IN ('NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 
                                 'NC_000005.9',  'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 
                                 'NC_000009.11', 'NC_000010.10', 'NC_000011.9',  'NC_000012.11', 
                                 'NC_000013.10', 'NC_000014.8',  'NC_000015.9',  'NC_000016.9', 
                                 'NC_000017.10', 'NC_000018.9',  'NC_000019.9',  'NC_000020.10', 
                                 'NC_000021.8',  'NC_000022.10', 'NC_000023.10', 'NC_000024.9')
                  AND tx_ac IN (
                      -- Subquery: Find transcripts that have a gap ON hg19
                      SELECT DISTINCT tx_ac 
                        FROM {self._uta_schema}.tx_exon_aln_v 
                       WHERE alt_aln_method = 'splign'
                         AND tx_ac LIKE 'NM_%'
                         AND cigar !~ '^[0-9]+=$'
                        -- Filter 2: Vital for speed; restricts the search space
                        AND alt_ac IN ('NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 
                                       'NC_000005.9',  'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 
                                       'NC_000009.11', 'NC_000010.10', 'NC_000011.9',  'NC_000012.11', 
                                       'NC_000013.10', 'NC_000014.8',  'NC_000015.9',  'NC_000016.9', 
                                       'NC_000017.10', 'NC_000018.9',  'NC_000019.9',  'NC_000020.10', 
                                       'NC_000021.8',  'NC_000022.10', 'NC_000023.10', 'NC_000024.9')
                  )
                ORDER BY tx_ac, ord;
            """

        with UtaDb() as uta_db:
            df = pd.read_sql(query, uta_db._hdp._conn)
            self._logger.info(f"Found {df.shape} transcript exons with cigar strings indicating mismatch")            
            concatenated_df = (
                df.groupby('tx_ac')['cigar']
                  .agg(' '.join)
                  .reset_index()
                  .sort_values('tx_ac')
                  .rename(columns={'tx_ac': 'accession'})
            )

        self._logger.info(f"Found {concatenated_df.shape} distinct transcripts with cigar strings indicating mismatch")
        return concatenated_df
    
    def write(self, output_file: str, transcript_gap_df):        
        transcript_gap_df.to_csv(output_file, index=False)
        self._logger.info(f"Wrote {transcript_gap_df.shape[0]} rows to {output_file}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Query UTA for alignment differences')
    parser.add_argument("--version", action="version", version="0.0.1")    
    parser.add_argument("--uta_schema", help="UTA db schema", required=True)
    parser.add_argument("--out_csv", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args
    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    uta_gap = ExtractUtaExonGapInfo(args.uta_schema)
    transcript_gap_df = uta_gap.get_transcripts_with_gaps()    
    uta_gap.write(args.out_csv, transcript_gap_df)
    
if __name__ == '__main__':
    main()