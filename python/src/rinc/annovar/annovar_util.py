'''
Created on Apr 1, 2026

@author: pleyte
'''
from collections import namedtuple

VariantTranscriptKey = namedtuple("VariantTranscriptKey", ["variant", "transcript"])

def get_variant_string(chromosome, position, ref, alt) -> str:
    '''
    Return a string represent the variant defined on the variant_function line
    '''
    return f"{chromosome}-{position}-{ref}-{alt}"

def get_variant_transcript_annotation(chromosome, position, ref, alt, cdna_transcript, variant_transcript_annotations):
    '''
    Return the variant function annotations from variant_transcript_annotations or None if none exist. 
    '''
    variant = get_variant_string(chromosome, position, ref, alt)
    variant_transcript_key = VariantTranscriptKey(variant, cdna_transcript)    
    return variant_transcript_annotations.get(variant_transcript_key)
    
def get_variant_function(variant_transcript_annotation: dict):
    '''
    Return the variant_function value (eg intronic, exonic) that was read by 
    the ParseAnnovarVariantFunction class.
    '''
    # The variant_function field is required 
    return variant_transcript_annotation['variant_function']

def get_is_splicing(variant_transcript_annotation):
    '''
    Return the 'splicing' value if it is present. 
    '''
    return variant_transcript_annotation.get('splicing')
    