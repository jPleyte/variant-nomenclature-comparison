'''
Convert csv to tab delimited SnpEff input file
Created on Jan 9, 2026

@author: pleyte
'''
import argparse
import csv
import logging.config
import pysam

import functools
from rinc.util.log_config import LogConfig
from rinc.util import chromosome_map

VERSION = '0.0.1'

@functools.total_ordering
class SimpleVariant:
    # Official genomic order for hg19/hg38
    _CHROM_ORDER = {str(i): i for i in range(1, 23)}
    _CHROM_ORDER.update({f"chr{i}": i for i in range(1, 23)})
    _CHROM_ORDER.update({"X": 23, "chrX": 23, "Y": 24, "chrY": 24, "MT": 25, "M": 25, "chrM": 25})

    def __init__(self, chromosome, position, reference, alt):
        self.chromosome = str(chromosome).replace('chr', '')
        self.position = int(position)
        self.reference = str(reference).upper()
        self.alt = str(alt).upper()

    def _sort_key(self):
        """Internal helper to return a (Rank, Position) tuple."""
        # Use 999 as a fallback for unknown scaffolds/contigs
        rank = self._CHROM_ORDER.get(self.chromosome, 999)
        return (rank, self.position, self.reference, self.alt)

    def __hash__(self):
        # Hash the tuple of all fields to ensure uniqueness in a set
        return hash((self.chromosome, self.position, self.reference, self.alt))

    def __eq__(self, other):
        if not isinstance(other, SimpleVariant):
            return NotImplemented
        return (self.chromosome, self.position, self.reference, self.alt) == \
               (other.chromosome, other.position, other.reference, other.alt)

    def __lt__(self, other):
        if not isinstance(other, SimpleVariant):
            return NotImplemented
        return self._sort_key() < other._sort_key()

    def __repr__(self):
        return f"SimpleVariant({self.chromosome}:{self.position} {self.reference}>{self.alt})"
    
class CsvToVcf(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def read(self, in_file_csv):
        """
        Read the csv and return list of variants 
        """         
        variant_set = set()
        
        with open(in_file_csv, mode='r', newline='', encoding='utf-8') as csvfile:
            for row in csv.DictReader(csvfile):
                v = SimpleVariant(row['chromosome'], row['position'], row['reference'], row['alt'])
                variant_set.add(v)
                
        return list(variant_set)
    
    def write(self, variants: list, out_file_vcf):
        """
        Write list of variants to VCF file 
        """
        header = pysam.VariantHeader()
        header.add_line("##fileformat=VCFv4.2")        
        header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
        
        for x in chromosome_map.refseq_to_ncbi.values():            
            header.contigs.add(x)

        with pysam.VariantFile(out_file_vcf, "w", header=header) as vcf_out:
            for x in variants:
                # pysam expects zero based position, and then adds one
                start_pos = x.position - 1
                rec = vcf_out.new_record(
                    contig=x.chromosome,
                    start=start_pos, 
                    stop=start_pos + len(x.reference),
                    alleles=(x.reference, x.alt),
                    id=".",
                    qual=None,
                    filter="PASS")

                vcf_out.write(rec)

def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a csv formatted file and write out a VCF')

    parser.add_argument('--in',
                dest="input",
                help='csv list of variants',
                required=True)

    parser.add_argument('--out',
                dest="output",
                help='Vcf out file',
                required=True)

    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    logger = logging.getLogger("rinc.etl.csv_to_vcf")

    args = _parse_args()

    c2v = CsvToVcf()
    variants = c2v.read(args.input)
    
    # This sort is not effective
    sorted_variants = sorted(variants)
    
    c2v.write(sorted_variants, args.output)

    logger.info(f"Wrote {len(variants)} variants to {args.output}")


if __name__ == '__main__':
    main()
