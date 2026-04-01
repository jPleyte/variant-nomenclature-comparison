#!/usr/bin/env nextflow

/*
 * Run ANNOVAR on the provided avinput file.
 */
process runAnnovar {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path annovar_avinput
    
    output:
    path "annovar.hg19_multianno.txt", emit: multianno    
    path "annovar.refGeneWithVer.variant_function", emit: refseq_variant_function
    path "annovar.refGeneWithVer.exonic_variant_function", emit: refseq_exonic_variant_function
    path "annovar.ccdsGene.variant_function", emit: ccds_variant_function
    path "annovar.ccdsGene.exonic_variant_function", emit: ccds_exonic_variant_function

    script:

    """    
    \$ANNOVAR_HOME/table_annovar.pl ${annovar_avinput} \
    \$ANNOVAR_HOME/humandb/ \
    --buildver hg19 \
    --out annovar \
    --protocol refGeneWithVer,ccdsGene \
    --operation g,g \
    --nastring . \
    --polish \
    --argument '--splicing_threshold 5 --exonicsplicing --transcript_function --separate,--splicing_threshold 5 --exonicsplicing --transcript_function --separate'
    """
}
