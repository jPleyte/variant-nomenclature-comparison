# variant-nomenclature-comparison

This project uses multipe tools to generates cDNA and protein nomenclature and produces a spreadsheet comparing the values from each. 

## Versions

* Python 3.14.3
* Nextflow 25.10.4 build 11173
* UTA uta_20241220
* SeqRepo 2024-12-20

## Dependencies

* vcf-sort from vcftools

## Running the workflow

### Running the workflow in Github Codespaces

This repsoitory has a ``.devcontainer`` directory that makes it able to be deployed as a GitHub Codespace. When the Codespace is first launched the ``.devcontainer/setup.sh`` script will install necessary dependencies including the hg19 fasta and annovar databases. 

To run the workflow change to the nextflow directory and launch the worfklow:
```bash
nextflow run .
```

### Variant list 

The workflow depends on a variant list. Each tool identifies transcripts and nomenclature for the variants. When the ``variant_source`` parameter is set to 'tfx' the variant list will be extracted from the variants fround in the output of the Transcript Effects tool.  


## Running the workflow locally 

To run the workflow on your local environment edit the ``nextflow/nextflow.config`` file and make the following changes to the ``mac`` profile:
* Set the ``params.fasta`` parameter to the full path to the hg19 fasta
* Set the ``params.uta_schema`` parameter to the schema in your UTA database
* Set the env ``ANNOVAR_HOME`` to directory where you installed annovar. 

Example command line for launching the workflow specifying the ``mac`` profile

```bash
nextflow run main.nf -profile mac main \ 
    --variant_source tfx  \
    --variant_source_file /tmp/transcript_effects/out.json \
    --generate_variant_validator_batch true \
    --mutalyzer_nomenclature_file ../data/mutalyzer/mutalyzer_nomenclature.csv \
    --variant_validator_batch_results ../data/variant_validator/variant_validator_result.tsv
```
