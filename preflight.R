## OncoSeq preflight -- source() this, or paste the whole block at once.
## Wrapped in local({...}) so it either runs completely or not at all.
local({
    root <- if (file.exists("DESCRIPTION")) "." else
            if (file.exists(file.path("OncoSeq", "DESCRIPTION"))) "OncoSeq" else NA

    if (is.na(root)) {
        cat("NOT at the package root.\n")
        cat("  working directory: ", getwd(), "\n", sep = "")
        cat("  Fix: open OncoSeq.Rproj, or setwd() to the folder holding DESCRIPTION.\n")
        return(invisible(NULL))
    }

    dcf <- read.dcf(file.path(root, "DESCRIPTION"))
    cat("package: ", dcf[, "Package"], "  version: ", dcf[, "Version"],
        "\n", sep = "")
    cat("R:       ", as.character(getRversion()), "\n", sep = "")
    if (getRversion() < "4.4.0")
        cat("  WARNING: OncoSeq needs R >= 4.4.0\n")

    deps <- c("SummarizedExperiment", "S4Vectors", "DESeq2", "edgeR", "limma",
              "fgsea", "BiocParallel", "survival", "ggplot2", "ggrepel",
              "matrixStats", "httr2", "jsonlite")
    tools <- c("roxygen2", "testthat", "knitr", "rmarkdown", "BiocStyle",
               "BiocCheck", "devtools")

    installed <- rownames(installed.packages())
    missDeps <- setdiff(deps, installed)
    missTools <- setdiff(tools, installed)

    cat("dependencies: ", length(deps) - length(missDeps), "/", length(deps),
        "   verification tools: ", length(tools) - length(missTools), "/",
        length(tools), "\n", sep = "")

    allMiss <- c(missDeps, missTools)
    if (length(allMiss) == 0L) {
        cat("\ntoolchain complete - ready to verify\n")
    } else {
        cat("\nmissing: ", paste(allMiss, collapse = ", "), "\n", sep = "")
        cat("\nRun this to install them:\n\n")
        cat('  if (!requireNamespace("BiocManager", quietly = TRUE))\n')
        cat('      install.packages("BiocManager")\n')
        cat("  BiocManager::install(c(",
            paste0('"', allMiss, '"', collapse = ", "), "))\n\n", sep = "")
    }
    invisible(NULL)
})
