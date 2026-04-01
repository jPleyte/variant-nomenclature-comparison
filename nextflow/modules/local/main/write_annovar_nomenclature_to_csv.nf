#!/usr/bin/env nextflow

/*
 * Convert csv variant list to annovar avinput format
 */
process writeAnnovarNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path multianno
    path refseq_variant_function
    path ccds_variant_function

    output:
    path "annovar_nomenclature.csv", emit: annovar_nomenclature

    script:
    """
    python -m rinc.annovar.parse_annovar_multianno               \
           --annovar_multianno ${multianno}                      \
           --annovar_variant_function ${refseq_variant_function} \
           --annovar_variant_function ${ccds_variant_function}   \
           --out annovar_nomenclature.csv
    """
}
