#!/usr/bin/env nextflow

/*
 * Convert csv variant list to vcf
 */
process csvToVcf {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path variant_csv
    
    output:
    path "variants.vcf", emit: vcf

    script:
    """
    python -m rinc.etl.csv_to_vcf \
           --in ${variant_csv} \
           --out variants_unsorted.vcf
    
    vcf-sort variants_unsorted.vcf > variants.vcf
    """
}
