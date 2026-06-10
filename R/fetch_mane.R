#' Fetch MANE Select transcripts
#'
#' Retrieves MANE Select transcripts using AnnotationHub and ensembldb.
#' MANE Select transcripts are a single representative transcript per
#' protein-coding gene, jointly curated by NCBI and EMBL-EBI.
#'
#' @param txdb An \code{EnsDb} or \code{TxDb} object. If \code{NULL},
#'   the function fetches one from AnnotationHub.
#' @param ah_id Optional AnnotationHub ID (e.g. \code{"AH113665"}) for
#'   a specific EnsDb release. Overrides automatic lookup.
#' @param species Species identifier. Default \code{"Homo sapiens"}.
#' @param assembly Genome assembly name. Default \code{"GRCh38"}.
#'
#' @return An \code{EnsDb} object containing only MANE Select transcripts,
#'   or the original \code{TxDb}/\code{EnsDb} if filtering is not
#'   possible.
#'
#' @importFrom ensembldb transcripts filter
#' @importFrom GenomicFeatures isActiveSeq
#' @export
#'
#' @examples
#' # Auto-fetch requires AnnotationHub (internet):
#' if (requireNamespace("AnnotationHub", quietly = TRUE) &&
#'     requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
#'     cat("Ready for MANE transcript fetching")
#' }
fetch_mane_txdb <- function(txdb = NULL, ah_id = NULL,
                            species = "Homo sapiens",
                            assembly = "GRCh38") {

    if (is.null(txdb)) {
        need_ah <- requireNamespace("AnnotationHub", quietly = TRUE)
        need_edb <- requireNamespace("ensembldb", quietly = TRUE)

        if (!need_ah) {
            stop("AnnotationHub is required to auto-fetch a TxDb. ",
                 "Install it with BiocManager::install('AnnotationHub').")
        }
        if (!need_edb) {
            stop("ensembldb is required for EnsDb support. ",
                 "Install it with BiocManager::install('ensembldb').")
        }

        ah <- AnnotationHub::AnnotationHub()

        if (!is.null(ah_id)) {
            txdb <- ah[[ah_id]]
        } else {
            # Query latest EnsDb for the species
            edb_records <- AnnotationHub::query(
                ah, c("EnsDb", species, assembly)
            )

            # Prefer current Ensembl release (> 100)
            rnums <- as.integer(sub(".*Ensembl ([0-9]+).*", "\\1",
                                    edb_records$title))
            rnums[is.na(rnums)] <- 0L
            latest_idx <- which.max(rnums)
            if (length(latest_idx) == 0 || rnums[latest_idx] == 0) {
                stop("No EnsDb found for ", species, " (", assembly, ")")
            }
            txdb <- edb_records[[names(edb_records)[latest_idx]]]
        }
    }

    # Try to filter to MANE Select transcripts
    txdb <- .filter_mane_select(txdb)
    txdb
}


#' @keywords internal
.filter_mane_select <- function(txdb) {
    # Check if this is an EnsDb (has filter support)
    if (!inherits(txdb, "EnsDb")) {
        warning("Input is not an EnsDb object. MANE Select filtering ",
                "requires EnsDb. Returning unfiltered TxDb.")
        return(txdb)
    }

    # Filter by MANE Select: look for the 'MANE Select' tag in
    # transcript support level (TSL) metadata.
    # In EnsDb, MANE Select is identified by tx_biotype or
    # a transcript tag. We try both approaches.

    # Approach 1: filter by MANE Select tag
    tx_filter <- ensembldb::TxidFilter(
        .get_mane_tx_ids(txdb)
    )

    if (is.null(tx_filter)) {
        warning("Could not identify MANE Select transcripts. ",
                "Returning all transcripts.")
        return(txdb)
    }

    n_before <- length(ensembldb::transcripts(txdb))
    txdb_filtered <- ensembldb::filter(txdb, tx_filter)
    n_after <- length(ensembldb::transcripts(txdb_filtered))

    message("Retained ", n_after, " MANE Select transcripts ",
            "(filtered from ", n_before, ").")
    txdb_filtered
}


#' @keywords internal
.get_mane_tx_ids <- function(txdb) {
    # Try known columns for MANE Select tags
    all_tx <- ensembldb::transcripts(
        txdb,
        columns = c("tx_id", "tx_biotype", "tx_support_level",
                    "gene_name", "uniprot_id")
    )

    mane_ids <- character(0)

    # Method 1: Check tx_biotype for MANE-related biotypes
    if ("tx_biotype" %in% colnames(all_tx)) {
        mane_biotypes <- grepl("MANE", all_tx$tx_biotype,
                               ignore.case = TRUE)
        if (any(mane_biotypes, na.rm = TRUE)) {
            mane_ids <- c(mane_ids,
                          all_tx$tx_id[mane_biotypes & !is.na(mane_biotypes)])
        }
    }

    # Method 2: Check if there is a tx_external_name or
    # transcript_name with 'NM_' prefix (RefSeq MANE Select)
    # MANE Select transcripts are always RefSeq NM_ accessions
    if ("tx_external_name" %in% colnames(all_tx)) {
        nm_tx <- grepl("^NM_", all_tx$tx_external_name)
        mane_ids <- c(mane_ids, all_tx$tx_id[nm_tx])
    }

    # Method 3: Try to get MANE tag from transcript metadata
    tx_cols <- ensembldb::listColumns(txdb)
    tag_cols <- grep("tag", names(tx_cols), ignore.case = TRUE,
                     value = TRUE)
    if (length(tag_cols) > 0) {
        for (col in tag_cols) {
            met <- ensembldb::transcripts(txdb,
                                          columns = c("tx_id", col))
            mane_tagged <- grepl("MANE.Select", met[[col]],
                                 ignore.case = TRUE)
            if (any(mane_tagged, na.rm = TRUE)) {
                mane_ids <- c(mane_ids,
                              met$tx_id[mane_tagged & !is.na(mane_tagged)])
            }
        }
    }

    mane_ids <- unique(mane_ids)

    if (length(mane_ids) == 0) {
        return(NULL)
    }

    mane_ids
}
