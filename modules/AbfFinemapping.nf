#!/usr/bin/env nextflow

process FinemapAbf {

    tag "${name}"
    
    // Define inputs
    input:
        tuple file(sumstats), file(reference)

    // Define outputs
    output:
        tuple path("*_finemapping_results.txt.gz"), path("*_finemapping_results_summary.txt")

    // Script to execute
    script:
    """
    pheno=\$(basename "$sumstats" | sed 's/^results_concat_\\(.*\\)\\.parquet\\.snappy\$/\\1/')
    echo \${pheno}

    Rscript ${baseDir}/bin/FinemapAbf.R \
    --pheno \${pheno} \
    --input ${sumstats} \
    --reference ${reference} \
    --outputpip \${pheno}_finemapping_results.txt.gz \
    --outputsummary \${pheno}_finemapping_results_summary.txt
    """
}
