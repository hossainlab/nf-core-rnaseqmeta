process DIFFERENTIAL_EXPRESSION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-limma=3.58.1 bioconda::bioconductor-edger=4.0.2 conda-forge::r-base=4.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(counts_matrix)
    path(sample_info)

    output:
    tuple val(meta), path("*_de_results.txt")    , emit: results
    tuple val(meta), path("*_de_plots.pdf")      , emit: plots
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(limma)
    library(edgeR)

    # Read count matrix and sample information
    counts <- read.table("${counts_matrix}", header=TRUE, row.names=1, sep="\\t")
    sample_info <- read.table("${sample_info}", header=TRUE, sep="\\t")

    # Ensure sample order matches
    sample_info <- sample_info[match(colnames(counts), sample_info\$sample_id), ]

    # Create DGEList object
    dge <- DGEList(counts = counts, samples = sample_info)

    # Filter low-expressed genes
    keep <- filterByExpr(dge, group = sample_info\$condition)
    dge <- dge[keep, , keep.lib.sizes = FALSE]

    # Normalize
    dge <- calcNormFactors(dge)

    # Create design matrix
    design <- model.matrix(~ 0 + condition, data = sample_info)
    colnames(design) <- levels(factor(sample_info\$condition))

    # Estimate dispersion
    dge <- estimateDisp(dge, design)

    # Fit model
    fit <- glmQLFit(dge, design)

    # Define contrasts (assuming case vs control)
    contrast_matrix <- makeContrasts(
        case_vs_control = case - control,
        levels = design
    )

    # Test for differential expression
    qlf <- glmQLFTest(fit, contrast = contrast_matrix)
    results <- topTags(qlf, n = Inf, sort.by = "PValue")

    # Create plots
    pdf("${prefix}_de_plots.pdf", width = 12, height = 8)
    
    # MA plot
    plotMD(qlf, main = "MA Plot")
    
    # Volcano plot
    with(results\$table, plot(logFC, -log10(PValue), 
                             pch = 20, 
                             main = "Volcano Plot",
                             xlab = "log2 Fold Change",
                             ylab = "-log10 P-value"))
    
    # Add significance thresholds
    abline(h = -log10(0.05), col = "red", lty = 2)
    abline(v = c(-1, 1), col = "red", lty = 2)
    
    # PCA plot
    logCPM <- cpm(dge, log = TRUE)
    pca <- prcomp(t(logCPM))
    plot(pca\$x[,1], pca\$x[,2], 
         col = as.factor(sample_info\$condition),
         pch = 19,
         main = "PCA Plot",
         xlab = paste0("PC1 (", round(summary(pca)\$importance[2,1]*100, 1), "%)"),
         ylab = paste0("PC2 (", round(summary(pca)\$importance[2,2]*100, 1), "%)"))
    legend("topright", 
           legend = levels(as.factor(sample_info\$condition)),
           col = 1:length(levels(as.factor(sample_info\$condition))),
           pch = 19)
    
    dev.off()

    # Write results
    write.table(results\$table, 
                file = "${prefix}_de_results.txt", 
                sep = "\\t", 
                quote = FALSE, 
                row.names = TRUE, 
                col.names = NA)

    # Write versions
    writeLines(c('"${task.process}":',
                 paste('    limma:', packageVersion("limma")),
                 paste('    edgeR:', packageVersion("edgeR"))),
               "versions.yml")
    """
}