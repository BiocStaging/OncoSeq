#' @importFrom httr2 request req_url_path_append req_url_query req_body_json
#'   req_retry req_timeout req_user_agent req_perform resp_body_json
#'   resp_body_string req_method
#' @importFrom jsonlite fromJSON
NULL

CBIOPORTAL_URL <- "https://www.cbioportal.org/api"

.cbioRequest <- function(path, query = list(), body = NULL,
                         baseUrl = CBIOPORTAL_URL, timeout = 300) {
    req <- httr2::request(baseUrl)
    req <- httr2::req_url_path_append(req, path)
    if (length(query) > 0L) req <- do.call(httr2::req_url_query,
                                           c(list(req), query))
    req <- httr2::req_user_agent(req, "OncoSeq R package")
    req <- httr2::req_timeout(req, timeout)
    req <- httr2::req_retry(req, max_tries = 4L)
    if (!is.null(body)) {
        req <- httr2::req_body_json(req, body, auto_unbox = FALSE)
        req <- httr2::req_method(req, "POST")
    }
    resp <- httr2::req_perform(req)
    jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyDataFrame = TRUE)
}

#' List studies available from cBioPortal
#'
#' @description Queries a cBioPortal instance for the studies it serves,
#'   optionally filtered by a keyword matched against the study identifier,
#'   name and cancer type.
#'
#' @param keyword Optional case-insensitive search term, for example
#'   `"bladder"`.
#' @param baseUrl Base URL of the cBioPortal API.
#'
#' @return A `data.frame` with one row per study, including `studyId`,
#'   `name`, `cancerTypeId` and the sample counts.
#'
#' @examples
#' \dontrun{
#' bladder <- cbioStudies("bladder")
#' head(bladder[, c("studyId", "name", "allSampleCount")])
#' }
#'
#' @seealso [cbioProfiles()], [importCBioPortal()]
#' @export
cbioStudies <- function(keyword = NULL, baseUrl = CBIOPORTAL_URL) {
    studies <- .cbioRequest("studies", baseUrl = baseUrl)
    studies <- as.data.frame(studies)
    if (!is.null(keyword)) {
        hay <- tolower(paste(studies$studyId, studies$name,
                             studies$cancerTypeId))
        studies <- studies[grepl(tolower(keyword), hay, fixed = TRUE), ,
                           drop = FALSE]
    }
    studies
}

#' List the molecular profiles of a cBioPortal study
#'
#' @param studyId cBioPortal study identifier, for example
#'   `"blca_tcga_pan_can_atlas_2018"`.
#' @param baseUrl Base URL of the cBioPortal API.
#'
#' @return A `data.frame` of molecular profiles with their
#'   `molecularProfileId`, `molecularAlterationType` and `datatype`.
#'
#' @examples
#' \dontrun{
#' cbioProfiles("blca_tcga_pan_can_atlas_2018")[, c("molecularProfileId",
#'                                                  "molecularAlterationType")]
#' }
#'
#' @export
cbioProfiles <- function(studyId, baseUrl = CBIOPORTAL_URL) {
    as.data.frame(.cbioRequest(paste0("studies/", studyId,
                                      "/molecular-profiles"),
                               baseUrl = baseUrl))
}

#' Import a cBioPortal study
#'
#' @description
#' Assembles a curated cBioPortal study into an
#' \linkS4class{OncoSeqExperiment}: an mRNA expression matrix, the somatic
#' mutation table, discrete (GISTIC-style) gene-level copy-number calls and
#' the clinical table, all restricted to samples shared by the requested
#' layers.
#'
#' @details
#' cBioPortal serves processed expression values rather than raw counts: for
#' TCGA PanCancer Atlas studies the default profile holds RSEM-normalised
#' values. The resulting assay is therefore named `expression`, not `counts`,
#' and count-based differential-expression engines will refuse to run on it;
#' use `engine = "limma"` on log-transformed values, or [importGDC()] when
#' raw counts are required.
#'
#' Clinical attributes arrive in long form and are pivoted to one column per
#' attribute. Overall survival is harmonised into `os_time` and `os_event`
#' when the study carries `OS_MONTHS`/`OS_STATUS`.
#'
#' @param studyId cBioPortal study identifier; see [cbioStudies()].
#' @param genes Optional character vector of HUGO symbols to restrict the
#'   expression and copy-number layers to. `NULL` fetches every gene in the
#'   profile, which for a whole study is a large download.
#' @param expressionProfile Molecular profile identifier for expression.
#'   When `NULL` the first continuous mRNA profile of the study is used.
#' @param mutationProfile Molecular profile identifier for mutations, or
#'   `NULL` to use the study's extended mutation profile. Set `mutations` to
#'   `FALSE` to skip.
#' @param cnvProfile Molecular profile identifier for discrete copy number,
#'   or `NULL` to use the study's GISTIC profile. Set `cnv` to `FALSE` to
#'   skip.
#' @param mutations,cnv,clinical Logical flags controlling which optional
#'   layers are retrieved.
#' @param chunkSize Number of genes per API request.
#' @param baseUrl Base URL of the cBioPortal API.
#'
#' @return An \linkS4class{OncoSeqExperiment} with an `expression` assay.
#'
#' @examples
#' \dontrun{
#' ose <- importCBioPortal("blca_tcga_pan_can_atlas_2018",
#'                         genes = c("TP53", "RB1", "FGFR3", "ERBB2"))
#' ose
#' }
#'
#' @seealso [importGDC()], [cbioStudies()]
#' @export
importCBioPortal <- function(studyId, genes = NULL,
                             expressionProfile = NULL,
                             mutationProfile = NULL, cnvProfile = NULL,
                             mutations = TRUE, cnv = TRUE, clinical = TRUE,
                             chunkSize = 500L, baseUrl = CBIOPORTAL_URL) {
    profiles <- cbioProfiles(studyId, baseUrl = baseUrl)
    if (nrow(profiles) == 0L) {
        stop("study ", sQuote(studyId), " has no molecular profiles; check ",
             "the identifier with cbioStudies()", call. = FALSE)
    }

    expressionProfile <- expressionProfile %||% .pickProfile(profiles,
        "MRNA_EXPRESSION", "CONTINUOUS", required = TRUE,
        what = "expression")

    geneMap <- .cbioGenes(genes, baseUrl = baseUrl)
    sampleIds <- .cbioSampleIds(studyId, expressionProfile, baseUrl = baseUrl)

    exprMat <- .cbioMolecularMatrix(expressionProfile, geneMap, sampleIds,
                                    chunkSize = chunkSize, baseUrl = baseUrl)
    if (nrow(exprMat) == 0L) {
        stop("no expression values were returned for study ", sQuote(studyId),
             call. = FALSE)
    }
    sampleIds <- colnames(exprMat)

    mutDF <- NULL
    if (isTRUE(mutations)) {
        mp <- mutationProfile %||% .pickProfile(profiles, "MUTATION_EXTENDED",
            NULL, required = FALSE, what = "mutation")
        if (!is.null(mp)) {
            mutDF <- .cbioMutations(mp, geneMap, sampleIds,
                                    chunkSize = chunkSize, baseUrl = baseUrl)
        }
    }

    cnvMat <- NULL
    if (isTRUE(cnv)) {
        cp <- cnvProfile %||% .pickProfile(profiles, "COPY_NUMBER_ALTERATION",
            "DISCRETE", required = FALSE, what = "copy number")
        if (!is.null(cp)) {
            cnvMat <- .cbioMolecularMatrix(cp, geneMap, sampleIds,
                                           chunkSize = chunkSize,
                                           baseUrl = baseUrl)
            if (nrow(cnvMat) == 0L) cnvMat <- NULL
        }
    }

    cd <- data.frame(sampleId = sampleIds,
                     patientId = sub("-[0-9]+$", "", sampleIds),
                     row.names = sampleIds, stringsAsFactors = FALSE)
    if (isTRUE(clinical)) {
        clinDF <- .cbioClinical(studyId, baseUrl = baseUrl)
        if (!is.null(clinDF)) {
            idx <- match(cd$patientId, clinDF$patientId)
            keep <- setdiff(colnames(clinDF), c(colnames(cd), "patientId"))
            cd <- cbind(cd, clinDF[idx, keep, drop = FALSE])
            rownames(cd) <- sampleIds
        }
    }
    cd <- .harmoniseCBioSurvival(cd)

    OncoSeqExperiment(
        assays = list(expression = exprMat),
        colData = cd,
        mutations = mutDF,
        cnv = cnvMat,
        provenance = list(source = "cBioPortal", studyId = studyId,
                          expressionProfile = expressionProfile,
                          baseUrl = baseUrl, retrieved = Sys.time()))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.pickProfile <- function(profiles, alterationType, datatype, required,
                         what) {
    ok <- profiles$molecularAlterationType == alterationType
    if (!is.null(datatype)) ok <- ok & profiles$datatype == datatype
    hits <- profiles$molecularProfileId[ok]
    if (length(hits) == 0L) {
        if (required) {
            stop("the study has no ", what, " profile of type ",
                 alterationType, call. = FALSE)
        }
        message("no ", what, " profile found; that layer will be empty")
        return(NULL)
    }
    hits[1L]
}

.cbioGenes <- function(genes, baseUrl) {
    if (is.null(genes)) {
        all_genes <- .cbioRequest("genes",
            query = list(projection = "SUMMARY", pageSize = 100000L),
            baseUrl = baseUrl)
        all_genes <- as.data.frame(all_genes)
        all_genes <- all_genes[all_genes$type == "protein-coding", ,
                               drop = FALSE]
        return(all_genes)
    }
    res <- .cbioRequest("genes/fetch",
        query = list(geneIdType = "HUGO_GENE_SYMBOL"),
        body = as.character(genes), baseUrl = baseUrl)
    res <- as.data.frame(res)
    missing_genes <- setdiff(genes, res$hugoGeneSymbol)
    if (length(missing_genes) > 0L) {
        warning(length(missing_genes), " gene symbol(s) not recognised by ",
                "cBioPortal, e.g. ",
                paste(head(missing_genes, 3L), collapse = ", "),
                call. = FALSE)
    }
    res
}

.cbioSampleIds <- function(studyId, profileId, baseUrl) {
    listId <- paste0(profileId)
    res <- try(.cbioRequest(paste0("sample-lists/", listId),
                            baseUrl = baseUrl), silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res$sampleIds)) {
        return(as.character(unlist(res$sampleIds)))
    }
    res <- .cbioRequest(paste0("sample-lists/", studyId, "_all"),
                        baseUrl = baseUrl)
    as.character(unlist(res$sampleIds))
}

.cbioMolecularMatrix <- function(profileId, geneMap, sampleIds, chunkSize,
                                 baseUrl) {
    ids <- geneMap$entrezGeneId
    chunks <- split(ids, ceiling(seq_along(ids) / chunkSize))
    pieces <- vector("list", length(chunks))
    for (k in seq_along(chunks)) {
        res <- .cbioRequest(
            paste0("molecular-profiles/", profileId, "/molecular-data/fetch"),
            query = list(projection = "SUMMARY"),
            body = list(entrezGeneIds = as.integer(chunks[[k]]),
                        sampleIds = as.character(sampleIds)),
            baseUrl = baseUrl)
        if (length(res) > 0L && is.data.frame(res) && nrow(res) > 0L) {
            pieces[[k]] <- res[, c("entrezGeneId", "sampleId", "value")]
        }
    }
    long <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
    if (is.null(long) || nrow(long) == 0L) {
        return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }
    symbol <- geneMap$hugoGeneSymbol[match(long$entrezGeneId,
                                           geneMap$entrezGeneId)]
    .longToMatrix(symbol, long$sampleId, suppressWarnings(
        as.numeric(long$value)))
}

.longToMatrix <- function(rowKey, colKey, value) {
    rowLevels <- sort(unique(rowKey))
    colLevels <- sort(unique(colKey))
    mat <- matrix(NA_real_, nrow = length(rowLevels),
                  ncol = length(colLevels),
                  dimnames = list(rowLevels, colLevels))
    idx <- cbind(match(rowKey, rowLevels), match(colKey, colLevels))
    mat[idx] <- value
    mat
}

.cbioMutations <- function(profileId, geneMap, sampleIds, chunkSize,
                           baseUrl) {
    ids <- geneMap$entrezGeneId
    chunks <- split(ids, ceiling(seq_along(ids) / chunkSize))
    pieces <- vector("list", length(chunks))
    wanted <- c("sampleId", "patientId", "entrezGeneId", "mutationType",
                "variantType", "chr", "startPosition", "endPosition",
                "referenceAllele", "variantAllele", "proteinChange",
                "mutationStatus", "tumorAltCount", "tumorRefCount")
    for (k in seq_along(chunks)) {
        res <- .cbioRequest(
            paste0("molecular-profiles/", profileId, "/mutations/fetch"),
            query = list(projection = "DETAILED"),
            body = list(entrezGeneIds = as.integer(chunks[[k]]),
                        sampleIds = as.character(sampleIds)),
            baseUrl = baseUrl)
        if (is.data.frame(res) && nrow(res) > 0L) {
            keep <- intersect(wanted, colnames(res))
            df <- res[, keep, drop = FALSE]
            df$gene <- geneMap$hugoGeneSymbol[match(res$entrezGeneId,
                                                    geneMap$entrezGeneId)]
            pieces[[k]] <- df
        }
    }
    long <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
    if (is.null(long) || nrow(long) == 0L) return(NULL)
    names(long)[names(long) == "sampleId"] <- "sample"
    .normaliseMutationTable(long)
}

.cbioClinical <- function(studyId, baseUrl, pageSize = 100000L) {
    out <- list()
    for (type in c("PATIENT", "SAMPLE")) {
        page <- 0L
        repeat {
            res <- .cbioRequest(paste0("studies/", studyId, "/clinical-data"),
                query = list(clinicalDataType = type, projection = "SUMMARY",
                             pageSize = pageSize, pageNumber = page),
                baseUrl = baseUrl)
            if (!is.data.frame(res) || nrow(res) == 0L) break
            out[[length(out) + 1L]] <- data.frame(
                patientId = res$patientId,
                attribute = res$clinicalAttributeId,
                value = res$value, stringsAsFactors = FALSE)
            if (nrow(res) < pageSize) break
            page <- page + 1L
        }
    }
    long <- do.call(rbind, out)
    if (is.null(long) || nrow(long) == 0L) return(NULL)
    long <- long[!duplicated(long[, c("patientId", "attribute")]), ,
                 drop = FALSE]
    wide <- stats::reshape(long, idvar = "patientId", timevar = "attribute",
                           direction = "wide")
    colnames(wide) <- sub("^value\\.", "", colnames(wide))
    rownames(wide) <- NULL
    .coerceNumericColumns(wide)
}

## Convert columns that are entirely numeric strings to numeric.
.coerceNumericColumns <- function(df) {
    for (nm in setdiff(colnames(df), "patientId")) {
        v <- df[[nm]]
        if (!is.character(v)) next
        num <- suppressWarnings(as.numeric(v))
        if (!all(is.na(num[!is.na(v) & nzchar(v)]))) {
            if (all(is.na(num) == (is.na(v) | !nzchar(v)))) df[[nm]] <- num
        }
    }
    df
}

## cBioPortal encodes status as "1:DECEASED" / "0:LIVING".
.harmoniseCBioSurvival <- function(cd) {
    cd <- as.data.frame(cd)
    if ("OS_STATUS" %in% colnames(cd)) {
        st <- as.character(cd$OS_STATUS)
        cd$os_event <- suppressWarnings(as.integer(sub(":.*$", "", st)))
    }
    if ("OS_MONTHS" %in% colnames(cd)) {
        cd$os_time <- suppressWarnings(as.numeric(cd$OS_MONTHS))
        cd$os_time_unit <- "months"
    }
    cd
}
