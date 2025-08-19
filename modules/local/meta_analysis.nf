process META_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-limma=3.58.1 bioconda::bioconductor-edger=4.0.2 conda-forge::r-base=4.3.1 conda-forge::r-metafor=4.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-87ca1a2dd8b2a8b7b9b6e1e5b5b5b5b5b5b5b5b5:latest' :
        'biocontainers/mulled-v2-87ca1a2dd8b2a8b7b9b6e1e5b5b5b5b5b5b5b5b5:latest' }"

    input:
    tuple val(meta), path(study_results, stageAs: "study_*.txt")
    path(sample_info)

    output:
    tuple val(meta), path("*_meta_analysis_results.txt"), emit: results
    tuple val(meta), path("*_meta_analysis_plots.pdf") , emit: plots
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(limma)
    library(edgeR)
    library(metafor)

    # Read study results
    study_files <- list.files(pattern="study_.*\\\\.txt", full.names=TRUE)
    study_results <- lapply(study_files, function(x) {
        read.table(x, header=TRUE, sep="\\t", row.names=1)
    })
    names(study_results) <- gsub("study_|\\\\.txt", "", basename(study_files))

    # Combine results for meta-analysis
    all_genes <- Reduce(intersect, lapply(study_results, rownames))
    
    # Extract effect sizes and standard errors
    effect_sizes <- sapply(study_results, function(x) x[all_genes, "logFC"])
    std_errors <- sapply(study_results, function(x) x[all_genes, "logFC"] / x[all_genes, "t"])

    # Perform meta-analysis for each gene
    meta_results <- data.frame(
        gene_id = all_genes,
        meta_logFC = numeric(length(all_genes)),
        meta_se = numeric(length(all_genes)),
        meta_pval = numeric(length(all_genes)),
        meta_padj = numeric(length(all_genes)),
        heterogeneity_pval = numeric(length(all_genes)),
        stringsAsFactors = FALSE
    )

    pdf("${prefix}_meta_analysis_plots.pdf", width=12, height=8)
    
    for (i in 1:length(all_genes)) {
        gene <- all_genes[i]
        
        # Meta-analysis using random effects model
        tryCatch({
            res <- rma(yi = effect_sizes[i, ], 
                      sei = abs(std_errors[i, ]), 
                      method = "REML")
            
            meta_results\$meta_logFC[i] <- res\$beta[1]
            meta_results\$meta_se[i] <- res\$se
            meta_results\$meta_pval[i] <- res\$pval
            meta_results\$heterogeneity_pval[i] <- res\$QEp
            
            # Create forest plot for top genes
            if (i <= 20) {
                forest(res, main = paste("Meta-analysis:", gene))
            }
        }, error = function(e) {
            meta_results\$meta_logFC[i] <<- NA
            meta_results\$meta_se[i] <<- NA
            meta_results\$meta_pval[i] <<- NA
            meta_results\$heterogeneity_pval[i] <<- NA
        })
    }
    
    dev.off()

    # Adjust p-values
    meta_results\$meta_padj <- p.adjust(meta_results\$meta_pval, method = "BH")

    # Sort by p-value
    meta_results <- meta_results[order(meta_results\$meta_pval), ]

    # Write results
    write.table(meta_results, 
                file = "${prefix}_meta_analysis_results.txt", 
                sep = "\\t", 
                quote = FALSE, 
                row.names = FALSE)

    # Write versions
    writeLines(c('"${task.process}":',
                 paste('    limma:', packageVersion("limma")),
                 paste('    edgeR:', packageVersion("edgeR")),
                 paste('    metafor:', packageVersion("metafor"))),
               "versions.yml")
    """
}