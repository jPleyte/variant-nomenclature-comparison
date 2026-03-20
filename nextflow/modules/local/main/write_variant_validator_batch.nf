/*
 * Generate variant and transcript lists for batch submission to VariantValidator. 
 * We only request variant-transcripts that have matching c. values in the cgd and tfx nomenclature files 
 */
process writeVariantValidatorBatch {
    publishDir "${params.outdir}/batch/variant_validator", mode: 'symlink'

    input:
    path cgd_nomenclature
    path tfx_nomenclature

    output:
    path "variant_validator_variants.txt", emit: variant_validator_variants
    path "variant_validator_transcripts.txt", emit: variant_validator_transcripts

    script:
    """
    python -m rinc.vv.variant_validator_batch           \
        --cgd_nomenclature ${cgd_nomenclature}          \
        --tfx_nomenclature ${tfx_nomenclature}          \
        --variant_output variant_validator_variants.txt \
        --transcript_output variant_validator_transcripts.txt
    """
}    
