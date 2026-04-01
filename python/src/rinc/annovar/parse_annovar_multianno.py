'''
Parse transcript details from an Annovar multianno file and write results to a csv.

A problem with Annovar's multianno format is that multiple transcripts are placed 
in the same row and it is not possible to determine which values from the "Func" 
column's list of functional classifications (eg, exonic, intronic, splicing) are
mean to be assigned to which transcript. This script accepts one or more "--annovar_variant_function"
paramters which are used to assocaite functional classification with transcript. 
  

it is assumed that Annovar was run using table_annovar.pl with the parameters below. 
The ``--remove`` parameter should not be used because this script parses both the 
multianno and variant_function files. 
```bash
$ANNOVAR_HOME/table_annovar.pl annovar.avinput \
    $ANNOVAR_HOME/humandb/ \
    --buildver hg19 \
    --out annovar \
    --protocol refGeneWithVer,ccdsGene \
    --operation g,g \
    --nastring . \
    --polish \
    --argument '--splicing_threshold 5 --exonicsplicing --transcript_function _update_variant_function,--splicing_threshold 5 --exonicsplicing --transcript_function --separate'
```

And that the resulting annovar.hg19_multianno.txt has the following fields:

| RefSeq                    | CCDS                | Explanation                                                             | Delimiters      |
| ------------------------- | ------------------- | ----------------------------------------------------------------------- | --------------- |
| Func.refGeneWithVer       | Func.ccdsGene       | Functional classifications variant (eg exonic or intronic)              | ; then ,        |
| Gene.refGeneWithVer       | Gene.ccdsGene       | List of transcript accessions                                           | ; then ,        |
| GeneDetail.refGeneWithVer | GeneDetail.ccdsGene | c. nomenclature for for variants located in non-coding regions          | ; then , then : |
| ExonicFunc.refGeneWithVer | ExonicFunc.ccdsGene | predicted biological impact of a variant on the protein-coding sequence | ;               |
| AAChange.refGeneWithVer   | AAChange.ccdsGene   | c. and p. nomenclature for variants located in coding regions           | ; then , then : |

     
Created on Jan 5, 2026

@author: pleyte
'''
import argparse
import csv
import logging.config
from rinc.util.log_config import LogConfig
from rinc.variant_transcript import VariantTranscript
import re
from rinc.annovar.parse_annovar_variant_function import ParseAnnovarVariantFunction
from rinc.annovar import annovar_util


class AnnovarMultianno(object):
    '''
    classdocs
    '''
    def __init__(self, variant_transcript_annotations: dict):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._variant_transcript_annotations = variant_transcript_annotations
    
    def get_variant_transcripts(self, annovar_multianno_file: str):
        """
        Reach each line of Annovar's crazy multianno tsv and extract each refseq and ccds transcript
        """
        variant_transcripts = [] 
        with open(annovar_multianno_file, mode='r', encoding='utf-8') as f:        
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                # Extract refseq transcripts                
                refseq_transcripts = self._parse_primary_blocks(row['Chr'], row['Start'], row['Ref'], row['Alt'], 
                                                                row['Func.refGeneWithVer'], 
                                                                row['Gene.refGeneWithVer'], 
                                                                row['GeneDetail.refGeneWithVer'], 
                                                                row['ExonicFunc.refGeneWithVer'], 
                                                                row['AAChange.refGeneWithVer'])

                variant_transcripts.extend(refseq_transcripts)

                # Extract ccds transcripts 
                ccds_transcripts = self._parse_primary_blocks(row['Chr'], row['Start'], row['Ref'], row['Alt'], 
                                                              row['Func.ccdsGene'], 
                                                              row['Gene.ccdsGene'], 
                                                              row['GeneDetail.ccdsGene'], 
                                                              row['ExonicFunc.ccdsGene'], 
                                                              row['AAChange.ccdsGene'])

                variant_transcripts.extend(ccds_transcripts)                

        return variant_transcripts
    
    def _parse_primary_blocks(self, chromosome, start, ref, alt, 
                              multianno_func, 
                              multianno_gene,                 
                              multianno_gene_detail,   
                              multianno_exonic_func,     
                              multianno_aa_change) -> list[VariantTranscript]:
        """
        Split the annovar fields at the top-most level using ';' and parse out the variant transcript definitions 
        """        
        variant_transcripts = []
       
        # Parse out all the non-coding transcripts
        non_coding_changes = []
        if multianno_gene_detail != '.': 
            for x in multianno_gene_detail.split(';'):
                results = self._parse_non_coding_region_changes(chromosome, start, ref, alt, x)
                if results:
                    non_coding_changes.extend(results)
                 
        variant_transcripts.extend(non_coding_changes)
        
        # Parse out all the coding transcripts
        coding_changes = []
        if multianno_aa_change != '.':
            for x in multianno_aa_change.split(';'):
                results = self._parse_coding_region_changes(chromosome, start, ref, alt, multianno_exonic_func, x)
                if results:
                    coding_changes.extend(results)
        
        variant_transcripts.extend(coding_changes)
                
        return variant_transcripts
            
    def _parse_non_coding_region_changes(self, chromosome, start, ref, alt, gene_detail):
        """
        Parse the comma delimited annovar GeneDetail field for the c. and transcript that do not affect the protein  
        """
        non_coding_change_transcripts = [] 
        for x in gene_detail.split(","):
            transcript = self._get_non_coding_transcript(chromosome, start, ref, alt, x)
            if transcript:
                non_coding_change_transcripts.append(transcript)
        
        return non_coding_change_transcripts

    def _get_non_coding_transcript(self, chromosome, start, ref, alt, gene_detail):
        """
        Parse transcript, c. and nearest exon out of the gene_detail field and return a VariantTranscript
        - . 
        - dist=44
        - NM_001345965.2:c.-38G>C
        - NM_058197.4:exon3:UTR3;
        - NM_001363763.2:exon3:c.305-3C>A;
        """
        # No use for intergenic
        if gene_detail.startswith("dist"):
            return None 
        
        c_dot = None
        exon = None
        
        parts = gene_detail.split(':')
        accession = parts[0]
        
        for x in parts:
            if x.startswith('c.'):
                c_dot = x
            elif x.startswith('exon'):
                exon = x.replace('exon', '')
        
        if c_dot:
            v = VariantTranscript(chromosome, start, ref, alt, accession)
            v.c_dot = self._get_normalized_c_dot(c_dot)
            v.exon = exon
            self._update_variant_function(v)            
            return v
        else:
            return None
    
    def _parse_coding_region_changes(self, chromosome, start, ref, alt, exonic_func, aa_change):
        """
        Parse the comma-delimited annovar field AAChange for the coding changes 
        """
        coding_change_transcripts = []
         
        for x in aa_change.split(","):
            transcript = self._get_coding_transcript(chromosome, start, ref, alt, exonic_func, x)
            if transcript:
                coding_change_transcripts.append(transcript)
        
        return coding_change_transcripts

    def _get_coding_transcript(self, chromosome: str, start: str, ref: str, alt: str, exonic_func: str, aa_change: str):
        """
        Parse transcript, c., p. and nearest exon out of the AAChange field and return a VariantTranscript
        Usually has p. but sometimes doesn't 
        - DDR2:NM_001014796.3:exon14:c.1722_1723delinsT:p.E576Rfs*20
        - TNFRSF14:NM_001297605.1:exon1:c.-13_1delinsC        
        """
        if aa_change == 'UNKNOWN':
            return None
            
        gene, accession, exon, c_dot, *rest = aa_change.split(':')
        
        if rest and len(rest) > 1:
            raise ValueError(f"annovar string has more than five parts: {aa_change}")
        elif rest:
            p_dot = rest[0]
        else:
            p_dot = None
        
        assert c_dot.startswith('c.'), f"AAChange nor formatted as expected: {aa_change}"
        
        v = VariantTranscript(chromosome, start, ref, alt, accession)
        v.c_dot = self._get_normalized_c_dot(c_dot)
        v.exon = exon.replace('exon', '')
        v.gene = gene
        v.p_dot1 = self._get_normalized_p_dot1(p_dot)        
        v.protein_variant_type = exonic_func
        v.additional_fields['p_raw_dot1'] = p_dot
        self._update_variant_function(v)        
        
        return v

    def _update_variant_function(self, vt: VariantTranscript):
        '''
        Update the VariantTranscript with information from the variant_function file 
        '''
        variant_transcript_annotation = annovar_util.get_variant_transcript_annotation(vt.chromosome, 
                                                                                       vt.position, 
                                                                                       vt.reference, 
                                                                                       vt.alt,
                                                                                       vt.cdna_transcript, 
                                                                                       self._variant_transcript_annotations)

        if variant_transcript_annotation:
            vt.genomic_region_type = annovar_util.get_variant_function(variant_transcript_annotation)
            vt.additional_fields['splicing'] = annovar_util.get_is_splicing(variant_transcript_annotation)
        
        
    def _get_normalized_p_dot1(self, p_dot_raw: str) -> str:
        """
        Normalize annovar's p.
        1) Convert X to *
        2) Use = instead of bases for equality (eg p.T73= rather than p.T73T)
        """
        if not p_dot_raw:
            return None
        
        # Convert X to * 
        normalized = re.sub(r'(\b|(?<=\w))X', '*', p_dot_raw)
        
        # Convert p.T73T to p.T73=
        pattern = r'^p\.([A-Z*])(\d+)\1$'
        normalized = re.sub(pattern, r'p.\1\2=', normalized)
        
        if normalized != p_dot_raw:
            self._logger.info(f"Normalized annovar pDot1 from {p_dot_raw} to {normalized}")
        return normalized
        
        
    def _get_normalized_c_dot(self, c_dot_raw) -> str:
        """
        Attempt to normalize annovar's c.
        Supported formats are:
        * Substitutions: c.C1307T --> c.1307C>T        
        """        
        substitution_pattern = r"c\.(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z])"
        match = re.search(substitution_pattern, c_dot_raw)
        if match:
            ref_base = match.group("ref")
            position = int(match.group("pos"))
            alt_base = match.group("alt")
            return f"c.{position}{ref_base}>{alt_base}"
        
        return c_dot_raw
    
    def write(self, output_file, variant_transcripts: list[VariantTranscript]):
        """
        Write file.
        Doesn't include g_dot or proteni transcript because annovar doesn't have them
        """
        key_headers = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript' ] 
        nomenclature_headers = ['c_dot', 'exon', 'gene', 'p_dot1', 'protein_variant_type', 'genomic_region_type', 'splicing', 'p_raw_dot1']
    
        all_headers = key_headers + nomenclature_headers 
    
        rows = 0
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
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
                       v.protein_variant_type,
                       v.genomic_region_type,
                       'splicing' if v.additional_fields.get('splicing') else None,
                       v.additional_fields['p_raw_dot1'] if 'p_raw_dot1' in v.additional_fields else None 
                       ]

                writer.writerow(row)
        
        self._logger.info(f"Wrote {len(variant_transcripts)} to {output_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract refseq and ccds transcripts and write to new csv')

    parser.add_argument('--annovar_multianno',
                        help='Multianno output from annovar (tsv)',
                        required=True)

    parser.add_argument('--annovar_variant_function',
                        help='Annovar variant_function file (tsv)',
                        action="append",
                        required=True)

    parser.add_argument("--out", help="output file (csv)", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()

    # Use variant_function file to get a list of variant-transcripts and their variant function (intronic, exonic, etc) and whether it is at a splice site
    pavf = ParseAnnovarVariantFunction()
    variant_transcript_annotations = pavf.get_variant_transcript_annotations(args.annovar_variant_function)

    # Read additiona variant-transcript information from the multianno file and combine it with values from the variant_function file    
    am = AnnovarMultianno(variant_transcript_annotations)
    variant_transcripts = am.get_variant_transcripts(args.annovar_multianno)

    am.write(args.out, variant_transcripts)


if __name__ == '__main__':
    main()
