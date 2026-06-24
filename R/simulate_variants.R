#' Simulate all possible SNVs for MANE Select transcripts
#'
#' Main entry point for the variant simulation pipeline. Fetches
#' MANE (Matched Annotation from NCBI and EMBL-EBI) Select
#' transcripts, extracts transcript structure, generates all possible
#' single nucleotide variants (SNVs) across CDS, UTR, and canonical
#' splice site regions, and formats them in HGVS notation.
#'
#' @param txdb An \code{EnsDb} or \code{TxDb} object. If \code{NULL},
#'   fetches from AnnotationHub (requires internet).
#' @param bsgenome A \code{BSgenome} object providing the reference
#'   genome sequence (e.g.,
#'   \code{BSgenome.Hsapiens.UCSC.hg38}).
#' @param transcript_ids Optional character vector of transcript IDs
#'   to simulate. \code{NULL} processes all MANE Select transcripts.
#'   Limiting is recommended for interactive use due to the large
#'   number of variants (~10,000 per transcript).
#' @param regions Character vector. Which regions to simulate.
#'   Options: \code{"cds"} (default), \code{"five_utr"},
#'   \code{"three_utr"}, \code{"splice_site"}.
#' @param use_messages Logical. If \code{TRUE}, prints progress
#'   messages. Default: \code{TRUE}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{transcript_id}{Transcript identifier.}
#'   \item{region}{Region type.}
#'   \item{genomic_pos}{Genomic position (1-based).}
#'   \item{genomic_ref}{Reference allele.}
#'   \item{genomic_alt}{Alternative allele.}
#'   \item{cds_pos}{CDS-relative position (NA for splice sites).}
#'   \item{exon_boundary}{Exon boundary coordinate (splice sites only).}
#'   \item{offset}{Intronic offset (splice sites only).}
#'   \item{hgvs_c}{HGVS coding DNA notation.}
#'   \item{hgvs_g}{HGVS genomic notation.}
#' }
#'
#' @seealso
#' \code{\link{fetch_mane_txdb}} for obtaining MANE Select transcripts.
#' \code{\link{get_transcript_structure}} for transcript structure
#' extraction.
#' \code{\link{generate_variants}} for the variant generation step.
#' \code{\link{format_hgvs}} for HGVS formatting details.
#'
#' @export
#'
#' @examples
#' # Check validity of HGVS strings:
#' is_valid_hgvs("NM_000546.6:c.215C>G")
#'
#' # Simulate requires matching TxDb + BSgenome with same seqlevels style:
#' # library(EnsDb.Hsapiens.v86)
#' # library(BSgenome.Hsapiens.UCSC.hg38)
#' # simulate_variants(EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38,
#' #                    transcript_ids = "ENST00000357654")
simulate_variants <- function(txdb = NULL,
                               bsgenome = NULL,
                               transcript_ids = NULL,
                               regions = c("cds", "five_utr",
                                           "three_utr",
                                           "splice_site"),
                               use_messages = TRUE) {

    # --- Validate inputs ---
    if (is.null(bsgenome)) {
        stop("bsgenome is required. Provide a BSgenome object, e.g., ",
             "BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38")
    }

    # --- Step 1: Fetch MANE Select transcripts ---
    if (use_messages) message("Step 1/4: Fetching MANE Select transcripts...")
    txdb <- fetch_mane_txdb(txdb)

    # --- Step 2: Extract transcript structure ---
    if (use_messages) message("Step 2/4: Extracting transcript structure...")
    tx_struct <- get_transcript_structure(txdb, transcript_ids)

    # --- Step 3: Generate all SNVs ---
    if (use_messages) {
        n_tx <- length(names(tx_struct$cds))
        message("Step 3/4: Generating variants for ", n_tx,
                " transcript(s)...")
    }
    variants <- generate_variants(tx_struct, bsgenome,
                                   regions = regions,
                                   transcript_ids = transcript_ids)

    # --- Step 4: Format HGVS ---
    if (use_messages) {
        message("Step 4/4: Formatting ", nrow(variants),
                " variants in HGVS notation...")
    }
    result <- format_hgvs(variants)

    if (use_messages) {
        message("Done. Generated ", nrow(result), " variants across ",
                length(unique(result$transcript_id)), " transcript(s).")
    }

    result
}
