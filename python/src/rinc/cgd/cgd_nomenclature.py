'''
Lookup variants in a saved export from the CGD database and write out all the variant
transcripts and nomenclature. 

Created on Jan 16, 2026

@author: pleyte
'''
import logging.config
import argparse
from rinc.io import variant_helper
from rinc.util.log_config import LogConfig
import pandas as pd
from rinc.variant_transcript import VariantTranscript
import csv
from rinc.cgd.variant_nomenclature_db import VariantNomenclatureDatabase

class CgdNomenclature(object):
    '''
    classdocs
    '''
    def __init__(self, cgd_db_file: str):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._cgd_db = VariantNomenclatureDatabase(cgd_db_file)
    
    def get_variants(self, variants_csv: str):
        """
        Read variants 
        """
        variants = variant_helper.get_variants(variants_csv)
        self._logger.info(f"Read {len(variants)} from {variants_csv}")
        return variants

    def _get_cgd_transcripts(self, v:VariantTranscript):
        """
        Return cgd transcripts for a variant 
        """
        return self._cgd_db.get_variant_transcripts(v.chromosome, v.position, v.reference, v.alt)
    
    def _get_cgd_to_variant_transcript(self, row):
        """
        Convert a CGD export row to a VariantTranscript object
        """            
        chromosome = row['chromosome'].replace('chr', '')
        
        vt = VariantTranscript(chromosome, row['position_start'], row['reference_base'], row['variant_base'], row['cdna_transcript'])
        vt.c_dot = None if pd.isna(row['genotype_cdna']) else row['genotype_cdna']
        vt.exon = None if pd.isna(row['exon']) else row['exon']
        vt.gene = None if pd.isna(row['cdna_gene']) else row['cdna_gene'] 
        vt.genomic_region_type = None if pd.isna(row['genomic_region_type']) else row['genomic_region_type'] 
        vt.p_dot1 = None if pd.isna(row['genotype_amino_acid_onel']) else row['genotype_amino_acid_onel']
        vt.p_dot3 =  None if pd.isna(row['genotype_amino_acid_threel']) else row['genotype_amino_acid_threel']
        vt.protein_transcript = None if pd.isna(row['protein_transcript']) else row['protein_transcript']
        vt.protein_variant_type = None if pd.isna(row['protein_variant_type']) else row['protein_variant_type']
        vt.additional_fields['splicing'] = None if pd.isna(row['splice_site']) else 'splicing'
        vt.additional_fields['genomic_variant_id'] = None if pd.isna(row['genomic_variant']) else row['genomic_variant']
        return vt
        
    def get_variant_transcripts(self, variants: list[VariantTranscript]):
        """
        Lookup the list of variants in the CGD db and return nomenclature for matches
        """
        variant_transcripts = []
        
        for v in variants:
            cgg_transcripts_df = self._get_cgd_transcripts(v)
            if not cgg_transcripts_df.empty:
                for _, row in cgg_transcripts_df.iterrows():
                    variant_transcripts.append(self._get_cgd_to_variant_transcript(row))
        
        return variant_transcripts

    def write(self, out_filename, variant_transcripts: list[VariantTranscript]):
        """
        Write variant transcripts to csv 
        Doesn't include g_dot because CGD doesn't have g_dot
        """
        # variant_helper.write_variant_transcripts(out_filename, variant_transcripts, [], 'cgd')
        key_headers = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript' ] 
        nomenclature_headers = ['c_dot', 'exon', 'gene', 'p_dot1', 'p_dot3', 'protein_transcript', 'protein_variant_type', 'splicing', 'genomic_variant_id']
        all_headers = key_headers + nomenclature_headers 
    
        rows = 0
        with open(out_filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(all_headers)
            
            for v in variant_transcripts:
                rows += 1
                row = [v.chromosome,
                       v.position,
                       v.reference,
                       v.alt,
                       v.cdna_transcript,
                       v.c_dot,
                       v.exon,
                       v.gene,
                       v.p_dot1,
                       v.p_dot3,
                       v.protein_transcript,
                       v.protein_variant_type,
                       v.additional_fields['splicing'],
                       v.additional_fields['genomic_variant_id']]
                
                writer.writerow(row)
        
        self._logger.info(f"Wrote {len(variant_transcripts)} variant transcripts to {out_filename}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Use hgvs to determine protein changes for a list of variants')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--cgd_db", help="File with variants and nomenclature exported from CGD (parquet)", required=True)
    parser.add_argument("--variants_input", help="Variant list", required=True)
    parser.add_argument("--out", help="CGD nomenclature output file (csv)", dest="output", required=True)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()

    cn = CgdNomenclature(args.cgd_db)
    
    variants = cn.get_variants(args.variants_input)
    variant_transcripts = cn.get_variant_transcripts(variants)
    cn.write(args.output, variant_transcripts)
    
if __name__ == '__main__':
    main()                
