#!/usr/bin/env nextflow

/*
 * Gather gap, hgvs, and annovar information; perform comparison, and write results to csv.
 */
process joinAndCompare {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    // Nomenclature dataframes    
    path annovar_nomenclature
    path snpeff_nomenclature
    path vep_refseq_nomenclautre
    path vep_hg19_nomenclature
    path tfx_nomenclature
    path cgd_nomenclature
    path variant_validator_nomenclature
    path mutalyzer_nomenclature
    
    // Additional inputs 
    path gff_and_uta_exon_gap_info
    path preferred_transcripts

    output:
    path "nomenclature_comparison.xlsx", emit: nomenclature_comparison_xlsx

    script:
    def tfx_nomenclature_arg = tfx_nomenclature ? "--tfx_nomenclature $tfx_nomenclature" : ""
    def cgd_nomenclature_arg = cgd_nomenclature ? "--cgd_nomenclature $cgd_nomenclature" : ""
    def vv_nomenclature_arg = variant_validator_nomenclature ? "--variant_validator_nomenclature $variant_validator_nomenclature" : ""
    def mutalyzer_nomenclature_arg = mutalyzer_nomenclature ? "--mutalyzer_nomenclature $mutalyzer_nomenclature" : ""
    def snpeff_nomenclature_arg = snpeff_nomenclature ? "--snpeff_nomenclature ${snpeff_nomenclature}" : ""
    def preferred_transcripts_arg = preferred_transcripts ? "--preferred_transcripts $preferred_transcripts" : ""
    """
    python -m rinc.join_and_compare                              \
        ${tfx_nomenclature_arg}                                  \
        ${cgd_nomenclature_arg}                                  \
        ${vv_nomenclature_arg}                                   \
        ${mutalyzer_nomenclature_arg}                            \
        ${snpeff_nomenclature_arg}                               \
        --annovar_nomenclature ${annovar_nomenclature}           \
        --vep_refseq_nomenclautre ${vep_refseq_nomenclautre}     \
        --vep_hg19_nomenclature ${vep_hg19_nomenclature}         \
        --gff_and_uta_exon_gap_info ${gff_and_uta_exon_gap_info} \
        ${preferred_transcripts_arg}                             \
        --out nomenclature_comparison.xlsx
    """
}
