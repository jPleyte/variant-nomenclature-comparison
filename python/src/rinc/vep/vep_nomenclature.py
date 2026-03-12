'''
Parse the file generate by VEP and write out a new csv 

@author: pleyte
'''

import logging.config

import argparse
from rinc.util.log_config import LogConfig
import pandas as pd
from rinc.variant_transcript import VariantTranscript
from collections import Counter
from rinc.util import chromosome_map
from rinc.util.pdot import PDot
import csv
import re
import urllib
from rinc.io import variant_helper

class VepNomenclature(object):
    '''
    Process VEP results
    '''
    def __init__(self, label):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._label = label
        self._p_dot_mapper = PDot()
        
        self._vep_df = None
        

    def read_vep_file(self, vep_file):
        """
        Read the VEP file.  The first 55 lines or so start with "##" and include meta data we can't use. 
        Then there's a line that starts with a single "#" which is the field names. 
        This function first reads the file to figure out how many lines to skip. 
        Then it reads the file using pandas, telling it to skip the lines starting with ##. 
        """
        with open(vep_file, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('#Uploaded_variation'):
                    header_line = i
                    break
        
        self._vep_df = pd.read_csv(vep_file, sep='\t', skiprows=header_line)
        self._logger.info(f"Read {self._vep_df.shape[0]} rows from {vep_file}")
    
    def _get_variant_values(self, row: pd.Series):
        """
        Return the variant and transcript parts. Pulls values from multiple fields and compares them to make sure they all match
        """
        if row['SOURCE'] == 'RefSeq':            
            transcript = row['Feature']            
        elif row['SOURCE'] == 'Ensembl' and row['CCDS'] != '-':
            transcript = row['CCDS']
        elif row['SOURCE'] == 'Ensembl' and row['Feature'].startswith('ENST'):
            transcript = row['Feature']
        else:
            raise ValueError(f"Unknown feature type: Source={row['SOURCE']}, Feature={row['Feature']}, CCDS={row['CCDS']}")
        
        chromosome, position, reference, alt = self._get_variant_values_from_gdot(row['#Uploaded_variation'])
        return chromosome, position, reference, alt, transcript
    
    def _get_variant_values_from_gdot(self, uploaded_variation: str):
        """
        Parse the variant values out of a Uploaded_variation field (eg 1_100908484_T/C) 
        """
        pattern = r"^(.+)_(\d+)_([ACGTN*.-]+)\/([ACGTN*.-]+)$"
        match = re.match(pattern, uploaded_variation)
                
        chrom, pos, ref, alt = match.groups()
        assert chrom and pos and ref and alt, "Variant parts could not be identified: "
        return chrom, pos, ref, alt
        
    def _get_c_dot(self, row):
        """
        Parse c. from the HGVSc field
        """
        if row['HGVSc'] == '-':
            return None

        transcript, c_dot = row['HGVSc'].split(':')
        
        if transcript.startswith('NM') and transcript != row['Feature']:
            raise ValueError(f"c. transcript does not match Feature: {transcript} != {row['Feature']}")
        elif not c_dot.startswith('c.'):
            raise ValueError(f"c. does not start with c.: {row['HGVSc']}")
        
        return c_dot
        
    def _get_g_dot(self, row):
        """
        Use HGVSg to return g.
        They use ncbi numbers for chromosome rather than refseq accessions so this function converts chromosome to refseq
        """
        chromosome, g_dot = row['HGVSg'].split(':')
        refseq_chromosome = chromosome_map.get_refseq(chromosome)
        return f"{refseq_chromosome}:{g_dot}"
     
    def _get_genomic_region_type(self, row) -> str:
        """
        Map the Consequence field to a one of our four genomic region types. 
        See https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
        """
        cons = row['Consequence'].lower()
        biotype = row['BIOTYPE'].lower()
        
        if 'upstream' in cons or 'downstream' in cons or 'intergenic' in cons:
            return 'intergenic'
        
        if "splice_acceptor_variant" in cons or "splice_donor_variant" in cons:
            return "splicing"
    
        if biotype == 'protein_coding':        
            exonic_terms = [ 'missense', 'synonymous', 'stop_gained', 'stop_lost', 
                            'frameshift', 'inframe_insertion', 'inframe_deletion', 
                            'coding_sequence_variant', 'protein_altering_variant', 
                            'start_lost', 'stop_retained_variant' ]
            
            if any(term in cons for term in exonic_terms):
                return 'exon'
        
            # UTR terms
            if 'utr_variant' in cons:
                return 'utr'
                
            # Intronic terms
            if 'intron_variant' in cons:
                # Special check: Essential Splice sites are intronic but critical
                if 'splice_donor' in cons or 'splice_acceptor' in cons:
                    return 'instron/splicing'
                return 'intron'
        
        
        # Non-coding or "Decay" transcripts
        # This handles biotypes like 'nonsense_mediated_decay' or 'antisense'
        if any(x in biotype for x in ['decay', 'antisense', 'rna', 'pseudogene', 'retained']):
            return 'non-coding'
    
        # Non-coding Exons (like in lincRNAs)
        if 'non_coding_transcript_exon_variant' in cons:
            return 'non-coding'
        
        raise ValueError(f"Unable to determine region type for {row['#Uploaded_variation']}: cons={cons}, biotype={biotype}")
    
    def _get_protein_variant_type(self, row):
        """
        Map VEP feilds to a protein variant type        
        """
        cons = row['Consequence'].lower()
        biotype = row['BIOTYPE'].lower()
    
        # 1. Essential Splice Loss (High Priority)
        if 'splice_donor' in cons or 'splice_acceptor' in cons:
            return 'Splice junction loss'
        
        # 2. Start/Stop Changes
        if 'start_lost' in cons: 
            return 'Start loss'
        if 'stop_gained' in cons: 
            return 'Stop gain'
        if 'stop_lost' in cons: 
            return 'Stop loss'
        
        # 3. Coding Changes
        if 'frameshift' in cons: 
            return 'Frameshift'
        
        # ITD Logic: Usually an inframe insertion where the sequence repeats
        # This is a placeholder; real ITD detection often needs the VCF ALT allele string
        if 'inframe_insertion' in cons:
            # If your data has a flag for duplication, use it here
            return 'In-frame'
        
        if 'inframe_deletion' in cons: 
            return 'In-frame'
        
        if 'missense' in cons:
            return 'Missense'
        
        if 'synonymous' in cons: 
            return 'Synonymous'
        
        if cons == 'stop_retained_variant':
            return 'Synonymous'
        
        # 4. MNV (Multi-nucleotide variant)
        # VEP often labels these as 'protein_altering_variant'
        if 'protein_altering' in cons: 
            return 'Multi nucleotide variant'
        
        # This isn't one of the CGD types but what else  to do with downstream_gene_variant?
        if 'downstream_gene_variant' in cons or 'intergenic' in cons:
            return 'Flanking'
    
        # 5. Promoter vs Flanking
        if 'upstream' in cons:
            # Standard convention: Promoter is within 2kb of the start
            # VEP provides 'DISTANCE' if you enable it
            return 'Flanking'
    
        # 6. Intron
        if 'intron' in cons: 
            return 'Intron'
        
        outside_terms = ['utr_variant', 'downstream', 'intergenic', 'non_coding_transcript_exon']
        if any(term in cons for term in outside_terms):
            return 'Flanking'
        
        # I don't know what to do with this one
        if cons == 'coding_sequence_variant' and biotype == 'protein_coding':
            return None     
        
        raise ValueError(f"Unable to determine protein variant type for {row['#Uploaded_variation']}: cons={cons}, biotype={biotype}")
        
        
    def _get_protein_change(self, row):
        """
        Return protein transcript and three letter p.
        """
        if row['HGVSp'] == '-':
            return None, None
        
        protein_transcript, p_dot3_raw = row['HGVSp'].split(':')
        
        # VEP uses "%3D" instead of "="
        p_dot3 = urllib.parse.unquote(p_dot3_raw)
        
        if row['SOURCE'] == 'RefSeq' and not protein_transcript.startswith("NP"):
            raise ValueError(f"Protein transcript does not start with NP: {row['HGVSp']}")
        elif not p_dot3.startswith("p."):
            raise ValueError(f"p. does not start with p.: {row['HGVSp']}")
        
        return protein_transcript, p_dot3
    
    def _get_exon(self, row):
        """
        Return the exon. VEP stores them in a string like '3/9'
        """
        if row['EXON'] == '-' and row['INTRON'] == '-':
            return None
        elif row['INTRON'] == '-':
            value = row['EXON']
        elif row['EXON'] == '-':
            value = row['INTRON']
        elif '/' in row['EXON'] and '/' in row['INTRON']:
            value = row['EXON']
        else: 
            raise ValueError(f"Unknown situation with exon/intron: {row['EXON']} {row['INTRON']}")
        
        return value.split('/')[0]
    
    def get_variant_transcripts(self) -> list[VariantTranscript]:
        """
        Parse the VEP dataframe        
        """
        transcript_counter = Counter()
        variant_transcripts = []
        
        for index, row in self._vep_df.iterrows():
            chromosome, position, reference, alt, transcript = self._get_variant_values(row)

            if transcript.startswith('NR'):
                transcript_counter['skipped_NR_transcript']
                continue
            elif transcript.startswith('ENST'):
                transcript_counter['skipped_ENST_transcript']
                continue
            
                         
            vt = VariantTranscript(chromosome, position, reference, alt, transcript)
            vt.c_dot = self._get_c_dot(row)
            vt.g_dot = self._get_g_dot(row)
            
            protein_transcript, p_dot3 = self._get_protein_change(row)
            vt.protein_transcript = protein_transcript
            vt.p_dot3 = p_dot3            
            if vt.p_dot3:
                vt.p_dot1 = self._p_dot_mapper.get_p_dot1(protein_transcript, p_dot3)
            
            vt.gene = row['SYMBOL'] if row['SYMBOL'] != '-' else None
            
            vt.genomic_region_type = self._get_genomic_region_type(row)
            vt.protein_variant_type = self._get_protein_variant_type(row)
            
            vt.exon = self._get_exon(row)
            
            vt.additional_fields['STRAND'] = row['STRAND']
            
            vt.additional_fields['REFSEQ_MATCH'] = row['REFSEQ_MATCH']
            if 'BAM_EDIT' in row: 
                vt.additional_fields['BAM_EDIT'] = row['BAM_EDIT']
            if 'GIVEN_REF' in row:
                vt.additional_fields['GIVEN_REF'] = row['GIVEN_REF']
                vt.additional_fields['USED_REF'] = row['USED_REF']
            
            transcript_counter['valid'] += 1
            variant_transcripts.append(vt)
            
        
        self._logger.info(f"Vep results: {transcript_counter}")
        return variant_transcripts
            
    def write(self, output_filename: str, variant_transcripts: list[VariantTranscript]):
        """
        """
        headers = ['chromosome', 'position', 'reference', 'alt',
                   'cdna_transcript', f'vep.{self._label}protein_transcript', 
                   f'vep.{self._label}.exon', f'vep.{self._label}.gene', f'vep.{self._label}.grt', f'vep.{self._label}.pvt',
                   f'vep.{self._label}.g_dot', f'vep.{self._label}.c_dot', f'vep.{self._label}.p_dot1', f'vep.{self._label}.p_dot3',
                   f'vep.{self._label}.REFSEQ_MATCH', f"vep.{self._label}.STRAND"]
        
        if variant_transcripts and 'BAM_EDIT' in variant_transcripts[0].additional_fields:
            headers.append(f'vep.{self._label}.BAM_EDIT')
        if variant_transcripts and 'GIVEN_REF' in variant_transcripts[0].additional_fields:
            headers.append(f'vep.{self._label}.GIVEN_REF')
            headers.append(f'vep.{self._label}.USED_REF')
        
        with open(output_filename, 'w', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(headers)

            unique_variant_transcripts = variant_helper._get_unique_varaint_transcripts(variant_transcripts)
            
            for v in unique_variant_transcripts:
                row = [v.chromosome, v.position, v.reference, v.alt, 
                       v.cdna_transcript, v.protein_transcript, 
                       v.exon, v.gene, v.genomic_region_type, v.protein_variant_type,
                       v.g_dot, v.c_dot, v.p_dot1, v.p_dot3,
                       v.additional_fields['REFSEQ_MATCH'], 
                       v.additional_fields['STRAND']]
                
                if 'BAM_EDIT' in variant_transcripts[0].additional_fields:
                    row.append(v.additional_fields['BAM_EDIT'])
                if 'GIVEN_REF' in variant_transcripts[0].additional_fields:
                    row.append(v.additional_fields['GIVEN_REF'])
                    row.append(v.additional_fields['USED_REF'])
                
                    
                writer.writerow(row)
            
        self._logger.info(f"Wrote {len(variant_transcripts)} variant transcripts to {output_filename}")
 
    
def _parse_args():
    parser = argparse.ArgumentParser(description='Read transcript nomenclature from vep output and write to csv')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--label", help="Label to include in all field names (eg 'refseq' or 'hg19')", required=True)
    parser.add_argument("--vep_results", help="File generated by VEP (tsv)", required=True)
    parser.add_argument("--out", help="VEP nomenclature output file (csv)", dest="output", required=True)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()

    vn = VepNomenclature(args.label)
    
    vn.read_vep_file(args.vep_results)
    variant_transcripts = vn.get_variant_transcripts()
    vn.write(args.output, variant_transcripts)
    
if __name__ == '__main__':
    main()                
        