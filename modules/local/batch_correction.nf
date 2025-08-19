process BATCH_CORRECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-sva=3.50.0 bioconda::bioconductor-limma=3.58.1 conda-forge::r-base=4.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(counts_matrix)
    path(sample_info)

    output:
    tuple val(meta), path("*_batch_corrected.txt"), emit: corrected_counts
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(sva)
    library(limma)

    # Read count matrix and sample information
    counts <- read.table("${counts_matrix}", header=TRUE, row.names=1, sep="\\t")
    sample_info <- read.table("${sample_info}", header=TRUE, sep="\\t")

    # Ensure sample order matches
    sample_info <- sample_info[match(colnames(counts), sample_info\$sample_id), ]

    # Create model matrix
    mod <- model.matrix(~ condition, data=sample_info)
    mod0 <- model.matrix(~ 1, data=sample_info)

    # Estimate batch effects using ComBat-seq
    if ("batch" %in% colnames(sample_info)) {
        # Use ComBat-seq for batch correction
        corrected_counts <- ComBat_seq(as.matrix(counts), 
                                     batch=sample_info\$batch, 
                                     group=sample_info\$condition)
    } else {
        # Use SVA to identify surrogate variables
        svobj <- sva(as.matrix(log2(counts + 1)), mod, mod0)
        corrected_counts <- removeBatchEffect(log2(counts + 1), 
                                            covariates=svobj\$sv, 
                                            design=mod)
        corrected_counts <- 2^corrected_counts - 1
    }

    # Write corrected counts
    write.table(corrected_counts, 
                file="${prefix}_batch_corrected.txt", 
                sep="\\t", 
                quote=FALSE, 
                row.names=TRUE, 
                col.names=NA)

    # Write versions
    writeLines(c('"${task.process}":',
                 paste('    sva:', packageVersion("sva")),
                 paste('    limma:', packageVersion("limma"))),
               "versions.yml")
    """
}