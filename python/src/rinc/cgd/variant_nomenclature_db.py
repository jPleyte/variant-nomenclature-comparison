'''
Created on Jan 18, 2026

@author: pleyte
'''

import argparse
import logging.config
import pandas as pd
from rinc.util.log_config import LogConfig
from pathlib import Path
from collections import Counter

class VariantNomenclatureDatabase(object):
    '''
    classdocs
    '''
    def __init__(self, db_file):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._vn_file = db_file
        self._schema = {'genomic_variant': pd.Series(dtype='Int64'), 
                        'chromosome': pd.Series(dtype='str'), 
                        'position_start': pd.Series(dtype='Int64'), 
                        'reference_base': pd.Series(dtype='str'), 
                        'variant_base': pd.Series(dtype='str'), 
                        'cdna_transcript': pd.Series(dtype='str'), 
                        'genotype_cdna': pd.Series(dtype='str'), 
                        'base_pair_position': pd.Series(dtype='str'), 
                        'protein_variant_type': pd.Series(dtype='str'),
                        'protein_transcript': pd.Series(dtype='str'),
                        'genotype_amino_acid_onel': pd.Series(dtype='str'),
                        'genotype_amino_acid_threel': pd.Series(dtype='str'),
                        'amino_acid_position': pd.Series(dtype='str'), 
                        'cdna_gene': pd.Series(dtype='str'), 
                        'exon': pd.Series(dtype='str'), 
                        'genomic_region_type': pd.Series(dtype='str'),                       
                        'splice_site': pd.Series(dtype='str')
                        }
        self._index_cols = ['chromosome', 'position_start', 'reference_base', 'variant_base', 'cdna_transcript']
        self._vn_df = self._get_or_create_db()
        self._vn_size = self._vn_df.shape[0] 
        
    
    def _get_or_create_db(self):
        """
        Load the database or create a new one
        """
        db_file = Path(self._vn_file)
        if db_file.exists():
            df = pd.read_parquet(self._vn_file).sort_index()
            self._logger.info(f"Read {df.shape[0]} rows from {self._vn_file}")
            return df
        else:            
            df = pd.DataFrame(self._schema)
            df.set_index(self._index_cols, drop=False, inplace=True)
            return df
            
    def _values_differ(self, val1, val2):
        # If both are null (None, NaN, etc), they do NOT differ
        if pd.isna(val1) and pd.isna(val2):
            return False
        # If one is null and the other isn't, they DO differ
        if pd.isna(val1) or pd.isna(val2):
            return True
        # Otherwise, do a standard comparison
        return val1 != val2
    
    def add_variant_transcripts_to_db(self, import_file: str, overwrite):
        """
        """
        incoming_data_counter = Counter()
        incoming_data_counter['initial_size'] = self._vn_df.shape[0]
        
        type_map = {col: series.dtype for col, series in self._schema.items()}
        incoming_df = pd.read_csv(import_file, 
                                  dtype=type_map)
        
        incoming_df.set_index(self._index_cols, drop=False, inplace=True)
        
        #Identify Conflicts: Intersection finds keys that exist in both DataFrames
        conflicts = incoming_df.index.intersection(self._vn_df.index)
        
        incoming_data_counter['match_count'] = len(conflicts)
        
        if not conflicts.empty:
            self._logger.info(f"Found {len(conflicts)} conflicting rows. Overwrite set to: {overwrite}")

            # Log the details of the conflicts
            for key in conflicts:
                old_row = self._vn_df.loc[key]
                new_row = incoming_df.loc[key]
                
                # Log the conflicting rows 
                c_diff = self._values_differ(old_row.get('genotype_cdna'), new_row.get('genotype_cdna'))
                p_diff = self._values_differ(old_row.get('genotype_amino_acid_onel'), new_row.get('genotype_amino_acid_onel'))
    
                if c_diff or p_diff:
                    incoming_data_counter['variant_transcript_match_nomenclature_mismatch'] += 1
                    msg = (f"Conflict for {key}\n"
                           f"  Existing: [c.: {old_row.get('genotype_cdna', 'N/A')}, p.: {old_row.get('genotype_amino_acid_onel', 'N/A')}]\n"
                           f"  Incoming: [c.: {new_row.get('genotype_cdna', 'N/A')}, p.: {new_row.get('genotype_amino_acid_onel', 'N/A')}]")
                    self._logger.info(msg)
            
            if overwrite:
                # Remove existing conflicts and then append all incoming
                self._vn_df = self._vn_df.drop(index=conflicts)
                self._vn_df = pd.concat([self._vn_df, incoming_df])
                incoming_data_counter['added_count'] = len(incoming_df)
            else:
                # SKIP: Only add rows that do NOT exist in the current index
                new_entries = incoming_df.drop(index=conflicts)
                self._vn_df = pd.concat([self._vn_df, new_entries])
                incoming_data_counter['added_count'] = len(new_entries)
        
        self._vn_df.sort_index(inplace=True)
            
        self._logger.info(f"Imported {import_file}")
        self._logger.info(f"results: {incoming_data_counter}")
    
        
    
    def _is_matching_nomenclature(self, incoming_row, row_key):
        """
        """
        existing_c_dot = self._vn_df.loc[row_key, 'genotype_cdna']
        incoming_c_dot = incoming_row['genotype_cdna']
        
        if not (existing_c_dot == incoming_c_dot):
            self._logger.info("Incoming row matches variant and transcript but has different nomenclature: ")
            return False
        else:
            return True
    
    def save_db(self):
        """
        Write the updated df out to Parquet file
        """        
        self._vn_df.to_parquet(self._vn_file)
        
        new_size = self._vn_df.shape[0]
        self._logger.info(f"Wrote {self._vn_file}. Initial size={self._vn_size}, new size={new_size}")
    
    def get_variant_transcripts(self, chromosome: str, position: int, reference: str, alt: str):
        """
        Query the variant transcript db for all all transcripts associated with a variant
        """
        chromosome = "chr"+ chromosome if not chromosome.startswith("chr") else chromosome
        
        try:
            match_df = self._vn_df.loc[(chromosome, position, reference, alt, slice(None)), :]
            return match_df
        except KeyError:
            return pd.DataFrame()
        
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Import variant transcripts and nomenclature into a db')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--db_file", help="Variant transcript db file to load or create (csv)", required=True)    
    parser.add_argument("--import", dest="import_file", help="SQL dump to import (csv)", required=True)
    parser.add_argument("--overwrite", action='store_true', help="Overwrite existing nomenclature values", required=False)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()
    
    vn_db = VariantNomenclatureDatabase(args.db_file)
    vn_db.add_variant_transcripts_to_db(args.import_file, args.overwrite)
    vn_db.save_db()
    
if __name__ == '__main__':
    main()                