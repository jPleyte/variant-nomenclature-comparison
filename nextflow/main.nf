#!/usr/bin/env nextflow

/*
This workflow gathers nomnenclature data from multiple tools, converts their proprietary formats to csv so they all use the same fields, and then compares them and writes
the results to an excel spreadsheet.

# To Do
- [ ] Delete "writeExonDetail" (around line 164) because exon/gap analysis has been moved to the main spreadsheet.  
- [ ] The "ch_all_labeled" channel is a mix of all data, but is not currently used. 
- [ ] Right now VariantValidator is limitd to the rows where cgd and tfx have the same c. value (see writeVariantValidatorBatch). Why not fetch all variant transcripts? 
- [x] You remove the filter code but still have a "ch_variant_transcript_filter" channel and paremter that aren't being used.
- [ ] The performAnalysis process pairwise_equality.py script is no longer working so i have commented it out in perform_analysis.nf
*/

// main processes
include { validateParameters; paramsSummaryLog } from 'plugin/nf-validation'
include { getTfxVariants } from './modules/local/variants/get_tfx_variants.nf'
include { csvToAvinput } from './modules/local/main/csv_to_avinput.nf'
include { runAnnovar } from './modules/local/main/run_annovar.nf'
include { writeAnnovarNomenclatureToCsv } from './modules/local/main/write_annovar_nomenclature_to_csv.nf'
include { csvToVcf } from './modules/local/main/csv_to_vcf.nf'
include { runSnpEff } from './modules/local/main/run_snpeff.nf'
include { writeSnpEffNomenclatureToCsv } from './modules/local/main/write_snpeff_nomenclature_to_csv.nf'
include { runVep as runVepRefseq } from './modules/local/main/run_vep.nf'
include { runVep as runVepHg19 } from './modules/local/main/run_vep.nf'
include { writeVepNomenclatureToCsv as writeVepRefseqNomenclatureToCsv } from './modules/local/main/write_vep_nomenclature_to_csv.nf'
include { writeVepNomenclatureToCsv as writeVepHg19NomenclatureToCsv } from './modules/local/main/write_vep_nomenclature_to_csv.nf'
include { writeTfxNomenclatureToCsv } from './modules/local/main/write_tfx_nomenclature_to_csv.nf'
include { writeCgdNomenclatureToCsv } from './modules/local/main/write_cgd_nomenclature_to_csv.nf'
include { joinAndCompare } from './modules/local/main/join_and_compare.nf'
include { performAnalysis } from './modules/local/main/perform_analysis.nf'
include { writeVariantValidatorNomenclature } from './modules/local/main/write_variant_validator_nomenclature.nf'
include { writeVariantValidatorBatch } from './modules/local/main/write_variant_validator_batch.nf'

// workflows 
include  { EXTRACT_EXON_GAP_INFO } from './subworkflows/local/extract_exon_gap_info.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    ncbi_refseq_gff_db = channel.fromPath(params.ncbi_refseq_gff_db, checkIfExists: true)
    preferred_transcripts = channel.fromPath(params.preferred_transcripts, checkIfExists: true)

    def vep_sequence_modes = [
        refseq: "--use_transcript_ref",
        reference: "--use_given_ref"        
    ]
    
    validateParameters()

    println "Using variant source: ${params.variant_source}"
    if (params.variant_source != 'gap_query') {
        println "Using variant source file: ${params.variant_source_file}"
    } 

    // Optional filter - list of variants and transcripts to keep 
    def filter_file_obj = params.variant_transcript_filter ? file(params.variant_transcript_filter, checkIfExists: true) : []
    // Print the status message
    if (filter_file_obj) {
        println "Using Variant Transcript Filter: ${filter_file_obj.name}"
    } else {
        println "No Variant Transcript Filter provided."
    }

    // Parameter validation
    if (params.variant_source_file != null && params.variant_source == 'gap_query') {
        error "You provided an input file but chose gap_query as the input method, which does not use a file."
    } else if (params.variant_source_file == null && (params.variant_source == 'csv' || params.variant_source == 'tfx')) {
        error "You indicated ${params.variant_source} as the variant_source, but did not provide a variant_source_file"
    }

    // Read variants 
    def ch_variants
    if (params.variant_source == 'csv') {
        ch_variants = params.variant_source_file
    }
    else if (params.variant_source == 'tfx') {
        ch_variants = getTfxVariants(params.variant_source_file)
    }
    else {
        error "Unknown variant source: ${params.variant_source}"
    }

	/* Annovar */     
    // Convert variant list to annovar avinput file 
    csvToAvinput(ch_variants)

    // Run annovar on the avinput file
    runAnnovar(csvToAvinput.out.annovar_avinput)

    // Extract annovar nomenclature and write to new csv
    def ch_annovar_nomenclature = writeAnnovarNomenclatureToCsv(runAnnovar.out.multianno)

    // Convert variant list to vcf to be used by SnpEff and VEP
    csvToVcf(ch_variants)

    // Run SnpEfff on vcf
    def ch_snpeff_nomenclature = channel.empty() 
    if(params.enable_snpeff) {
        runSnpEff(csvToVcf.out.vcf)

        // Extract SnpEff nomenclature and write to new csv
        ch_snpeff_nomenclature = writeSnpEffNomenclatureToCsv(runSnpEff.out.snpeff_tsv)
    }

    /* VEP */     
    // Run VEP onece using coding sequence for reference andusing hg19 for reference 
    runVepRefseq(csvToVcf.out.vcf, params.vep_fasta, 'refseq', vep_sequence_modes.refseq)
    runVepHg19(csvToVcf.out.vcf, params.vep_fasta, 'hg19', vep_sequence_modes.reference)
    
    // extract VEP nomenclature and write to new csv
    def ch_vepRefSeq_nomenclature = writeVepRefseqNomenclatureToCsv(runVepRefseq.out.vep_output, 'refseq')
    def ch_vepHg19_nomenclature = writeVepHg19NomenclatureToCsv(runVepHg19.out.vep_output, 'hg19')
    
    /* Transcript Effects */ 
    def ch_tfx_nomenclature = channel.empty()
    if (params.variant_source == 'tfx') {
        ch_tfx_nomenclature = writeTfxNomenclatureToCsv(fasta_ch, params.variant_source_file)
    }

    /* CGD */ 
    def ch_cgd_nomenclature = channel.empty()
    if (params.cgd_export_df != null) {
        ch_cgd_nomenclature = writeCgdNomenclatureToCsv(params.cgd_export_df, ch_variants)
    }

    /* Variant Validator */
    def ch_variant_validator_nomenclature = channel.empty()
    if (params.enable_variant_validator) {
        ch_variant_validator_batch_results = channel.fromPath(params.variant_validator_batch_results, checkIfExists: true)
        ch_variant_validator_nomenclature = writeVariantValidatorNomenclature(ch_variant_validator_batch_results)
    }

    /* Mutalyzer */ 
    def ch_mutalyzer_nomenclature = channel.empty()
    if (params.enable_mutalyzer) {
        // Run 
        if (params.mutalyzer_nomenclature_file) {
            ch_mutalyzer_nomenclature = channel.fromPath(params.mutalyzer_nomenclature_file, checkIfExists: true)
        }
    }

    def ch_all_labeled = ch_annovar_nomenclature.map { file -> ["annovar", file] }
        .mix( ch_snpeff_nomenclature.map  { file -> ["snpeff", file] } )
        .mix( ch_vepRefSeq_nomenclature.map { file -> ["vep_refseq", file] } )
        .mix( ch_vepHg19_nomenclature.map   { file -> ["vep_hg19", file] } )
        .mix( ch_tfx_nomenclature.map       { file -> ["tfx", file] } )
        .mix( ch_cgd_nomenclature.map               { file -> ["cgd", file] } )
        .mix( ch_variant_validator_nomenclature.map { file -> ["vv", file] } )
        .mix( ch_mutalyzer_nomenclature.map         { file -> ["mut", file] } )

    def ch_final_tool_outputs = ch_all_labeled.toList()

    // Pirint out the list of nomenclature files
    ch_final_tool_outputs.view { list -> 
        "Collected Outputs:\n" + list.collect { tuple -> "  - $tuple" }.join("\n") 
    }

    // Use a gff and the UTA db to create a file with all known transcripts that have reference gaps)
    def gff_and_uta_exon_gap_info = EXTRACT_EXON_GAP_INFO(ncbi_refseq_gff_db, uta_schema).gff_and_uta_exon_gap_info

    // jDebug: Delete this
    // Write out exon position and cigar strings for every transcript    
    // writeExonDetail(ncbi_refseq_gff_db, ncbi_refseq_gff_accession_index_df, ch_final_tool_outputs)    
    // error "STOPPING WORKFLOW for debuging"

    // Compare hgvs and annovar, join hgvs, annovar, and gaps file into final output
    joinAndCompare(ch_annovar_nomenclature,
                   ch_snpeff_nomenclature.ifEmpty([]),
                   ch_vepRefSeq_nomenclature,
                   ch_vepHg19_nomenclature,
                   ch_tfx_nomenclature.ifEmpty([]),
                   ch_cgd_nomenclature.ifEmpty([]),
                   ch_variant_validator_nomenclature.ifEmpty([]),
                   ch_mutalyzer_nomenclature.ifEmpty([]),
                   gff_and_uta_exon_gap_info,
                   preferred_transcripts)

    /*
    performAnalysis(ch_hgvs_nomenclature.ifEmpty([]),
                    ch_annovar_nomenclature,
                    ch_snpeff_nomenclature,
                    ch_vepRefSeq_nomenclature,
                    ch_vepHg19_nomenclature,
                    ch_tfx_nomenclature.ifEmpty([]),
                    ch_cgd_nomenclature.ifEmpty([]),
                    ch_mutalyzer_nomenclature.ifEmpty([])) */

    // Create a batch file of variants to submit to Variant Validator 
    if (params.generate_variant_validator_batch) {
        writeVariantValidatorBatch(ch_cgd_nomenclature, ch_tfx_nomenclature)
    }
}