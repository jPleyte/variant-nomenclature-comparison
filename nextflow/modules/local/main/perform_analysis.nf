#!/usr/bin/env nextflow

/*
 * Analyze differences in nomenclature between different tools 
 */
process performAnalysis {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path hgvs_nomenclature    
    path annovar_nomenclature    
    path snpeff_nomenclature
    path vep_refseq_nomenclautre
    path vep_hg19_nomenclature
    path tfx_nomenclature
    path cgd_nomenclature
    path mutalyzer_nomenclature

    output:
    path "analysis_pairwise_equality.csv", emit: analysis_pairwise_equality
    path "analysis_pairwise_equality.txt", emit: analysis_pairwise_equality_txt
    path "analysis_one_source_divergence.xlsx", emit: analysis_one_source_divergence
    
    script:
    def tfx_nomenclature_arg = tfx_nomenclature ? "--nomenclature tfx $tfx_nomenclature" : ""
    def hgvs_nomenclature_arg = hgvs_nomenclature ? "--nomenclature hgvs $hgvs_nomenclature" : ""
    def cgd_nomenclature_arg = cgd_nomenclature ? "--nomenclature cgd $cgd_nomenclature" : ""
    def mutalyzer_nomenclature_arg = mutalyzer_nomenclature ? "--nomenclature mutalyzer $mutalyzer_nomenclature" : ""
    """
    #python -m rinc.analysis.pairwise_equality \
    #    ${hgvs_nomenclature_arg} \
    #    ${tfx_nomenclature_arg} \
    #    ${cgd_nomenclature_arg} \
    #    ${mutalyzer_nomenclature_arg} \
    #    --nomenclature annovar ${annovar_nomenclature} \
    #    --nomenclature snpeff ${snpeff_nomenclature} \
    #    --nomenclature vepRefseq ${vep_refseq_nomenclautre} \
    #    --nomenclature vepHg19 ${vep_hg19_nomenclature} \
    #    --out analysis_pairwise_equality.csv
    
    python -c "import pandas as pd; from tabulate import tabulate; df=pd.read_csv('analysis_pairwise_equality.csv'); print(tabulate(df, headers='keys', tablefmt='psql', missingval=''))" > analysis_pairwise_equality.txt

    python -m rinc.analysis.one_source_divergence \
        ${hgvs_nomenclature_arg} \
        ${tfx_nomenclature_arg} \
        ${cgd_nomenclature_arg} \
        --nomenclature annovar ${annovar_nomenclature} \
        --nomenclature vepRefseq ${vep_refseq_nomenclautre} \
        --nomenclature vepHg19 ${vep_hg19_nomenclature} \
        --out analysis_one_source_divergence.xlsx
    """
}
