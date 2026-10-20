## Build data/simExample.rda, the simulated cohort returned by
## exampleOncoSeq().
##
## The simulation lives here rather than in R/ so that the shipped object
## is fixed once and the package calls no random number generator at run
## time: exampleOncoSeq() simply loads this object and subsets it.
##
## Run with:  Rscript data-raw/make-simExample.R

library(OncoSeq)

set.seed(2026L)

nPerGroup <- 20L
nDE <- 40L

drivers <- c("TP53", "RB1", "FGFR3", "ERBB2", "PIK3CA", "KDM6A",
             "ARID1A", "KMT2D", "STAG2", "CDKN2A", "MDM2", "EGFR",
             "CCND1", "E2F3", "KRT5", "KRT20", "PPARG", "GATA3",
             "FOXA1", "MKI67")
genes <- c(drivers, sprintf("GENE%03d", seq_len(280L)))
nGenes <- length(genes)

samples <- c(sprintf("T%02d", seq_len(nPerGroup)),
             sprintf("N%02d", seq_len(nPerGroup)))
condition <- factor(rep(c("tumour", "normal"), each = nPerGroup),
                    levels = c("normal", "tumour"))
isTumour <- condition == "tumour"

baseMu <- 2^runif(nGenes, 3, 11)
lfc <- rep(0, nGenes)
deIdx <- sample(seq_len(nGenes), nDE)
lfc[deIdx] <- rnorm(nDE, 0, 1.6)

mu <- outer(baseMu, rep(1, length(samples)))
mu[, isTumour] <- mu[, isTumour] * 2^lfc
dispersion <- 0.1 + 2 / sqrt(baseMu)
counts <- matrix(
    rnbinom(length(mu), mu = as.vector(mu),
            size = 1 / rep(dispersion, length(samples))),
    nrow = nGenes, dimnames = list(genes, samples))

mutGenes <- c("TP53", "KMT2D", "ARID1A", "KDM6A", "PIK3CA", "RB1",
              "FGFR3", "STAG2")
mutFreq <- c(0.48, 0.28, 0.25, 0.24, 0.20, 0.17, 0.14, 0.11)
classes <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
             "Splice_Site", "In_Frame_Del")
tumours <- samples[isTumour]
mutRows <- list()
for (k in seq_along(mutGenes)) {
    hit <- tumours[runif(length(tumours)) < mutFreq[k]]
    if (length(hit) == 0L) next
    mutRows[[k]] <- data.frame(
        sample = hit,
        gene = mutGenes[k],
        variant_class = sample(classes, length(hit), replace = TRUE,
                               prob = c(0.6, 0.15, 0.12, 0.08, 0.05)),
        chromosome = paste0("chr", ((k - 1L) %% 22L) + 1L),
        protein_change = paste0("p.X", sample(50:900, length(hit),
                                              replace = TRUE), "Y"),
        stringsAsFactors = FALSE)
}
mutations <- do.call(rbind, mutRows)

cnvGenes <- c("CDKN2A", "MDM2", "ERBB2", "E2F3", "CCND1", "EGFR",
              "PPARG", "RB1", "FGFR3", "TP53")
cnvMat <- matrix(0L, nrow = length(cnvGenes), ncol = length(samples),
                 dimnames = list(cnvGenes, samples))
for (g in cnvGenes) {
    direction <- if (g %in% c("CDKN2A", "RB1", "TP53")) -1L else 1L
    hit <- tumours[runif(length(tumours)) < 0.3]
    cnvMat[g, hit] <- direction *
        sample(1:2, length(hit), replace = TRUE, prob = c(0.4, 0.6))
}

tp53Mutant <- unique(mutations$sample[mutations$gene == "TP53"])
hazard <- ifelse(samples %in% tp53Mutant, 1 / 500, 1 / 1100)
trueTime <- rexp(length(samples), rate = hazard)
censorTime <- rexp(length(samples), rate = 1 / 1400)

colData <- data.frame(
    condition = condition,
    age = round(rnorm(length(samples), 68, 9)),
    stage = factor(sample(c("II", "III", "IV"), length(samples),
                          replace = TRUE, prob = c(0.35, 0.4, 0.25)),
                   levels = c("II", "III", "IV")),
    os_time = round(pmin(trueTime, censorTime)),
    os_event = as.integer(trueTime <= censorTime),
    row.names = samples)

simExample <- OncoSeqExperiment(
    assays = list(counts = counts),
    colData = colData,
    mutations = mutations,
    cnv = cnvMat,
    provenance = list(source = "simulation",
                      script = "data-raw/make-simExample.R",
                      seed = 2026L, nDE = nDE, built_on = Sys.Date()))

message("simExample: ", nrow(simExample), " genes x ", ncol(simExample),
        " samples; ", nrow(mutations(simExample)), " variants")
save(simExample, file = "data/simExample.rda", compress = "xz",
     compression_level = 9)
message("written (", round(file.size("data/simExample.rda") / 1024, 1),
        " KB)")
