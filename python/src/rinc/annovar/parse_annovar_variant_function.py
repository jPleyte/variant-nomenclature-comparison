'''
A class that parses the Annovar variant_function file

The ``get_variant_transcript_annotations`` function returns a dict where the key 
is a namedtuple object composed of a variant and transcript accession; and the 
value is information about the transcript. 

At the moment the only information about the transcript that is returned is:
* variant_function    eg exonic, intronic, UTR5, etc
* splicing            True when Annovar indicates that the transcript is in at a splice site. 

Created on Apr 1, 2026

@author: pleyte
'''
import logging
from collections import namedtuple
import re
from _collections import defaultdict
from rinc.annovar.annovar_util import VariantTranscriptKey
from rinc.annovar import annovar_util



class ParseAnnovarVariantFunction(object):
    '''
    Parse the Annovar variant_function file 
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def get_variant_transcript_annotations(self, variant_transcript_files: list[str]) -> dict:
        '''
        Parse one or more of Annovar's variant_transcript files to get the region type and a splicing indicator.         
        '''
        variant_transcript_annotations = defaultdict(dict)
        
        for variant_transcript_file in variant_transcript_files:
            with open(variant_transcript_file, 'r') as f:
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    variant = annovar_util.get_variant_string(fields[2], fields[3], fields[5], fields[6])
                    for transcript in self._get_variant_transcript_transcripts(fields):                        
                        variant_transcript_key = VariantTranscriptKey(variant, transcript)
                        variant_transcript_annotation = self._get_variant_transcript_annotation(fields)
                        variant_transcript_annotations[variant_transcript_key].update(variant_transcript_annotation)

        return variant_transcript_annotations    
        
    def _get_variant_transcript_transcripts(self, variant_function_fields: list):
        '''
        Return the list of transcripts defined on this line
        
        A variant_transcript line has 8 fields:
        0: function (eg exonic, splicing, intronic, downstream, UTR5)
        1: transcript info, can be in one of several formats:
            one transcript: NM_002524.5
            multiple transcript: NM_001172411.2,NM_001172412.1,NM_138959.3
            one transcript with info: NM_002524.5(NM_002524.5:exon3:c.112-2A>G)
            one or more transcripts with dist: NM_001007553.3,NM_001130523.3,NM_001242891.1,NM_001242892.2,NM_001242893.2,NM_007158.6(dist=865)
            multiple transcripts with info: NM_001170687.2(NM_001170687.2:exon8:c.997-1G>A),NM_001170688.1(NM_001170688.1:exon7:c.1015-1G>A),NM_001170689.2(NM_001170689.2:exon7:c.670-1G>A),NM_080875.3(NM_080875.3:exon8:c.1039-1G>A)
        '''
        # Split on commas that are NOT inside parentheses        
        gene_field = variant_function_fields[1]
        
        # Split on commas that are NOT inside parentheses
        entries = re.split(r',(?![^(]*\))', gene_field)
        
        transcripts = []         
        for entry in entries:
            match = re.match(r'([^(]+)(?:\(([^)]+)\))?', entry)
            if not match: 
                raise ValueError(f"Unable to parse line: {gene_field}")
            
            transcript = match.group(1).strip()
            transcripts.append(transcript)
        
        return transcripts
        
    def _get_variant_transcript_annotation(self, variant_function_fields: list):
        '''
        '''
        func = variant_function_fields[0]
        
        if func == 'splicing':
            return { 'splicing': True }
        else:
            return { 'variant_function': func }
