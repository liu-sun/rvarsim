#' Extract transcript structure for variant simulation
#'
#' Extracts the coding sequence (CDS), 5' UTR, 3' UTR, and canonical
#' splice site positions for one or more transcripts from a TxDb/EnsDb.
#'
#' Canonical splice sites are defined as the first and last two
#' positions of each intron (±1, ±2 from exon boundaries).
#'
#' @param txdb A \code{TxDb} or \code{EnsDb} object.
#' @param transcript_ids Optional character vector of transcript IDs
#'   to limit extraction. If \code{NULL}, all transcripts are used.
#'
#' @return A list with components:
#' \describe{
#'   \item{cds}{A \code{GRangesList} of CDS exons by transcript.}
#'   \item{five_utr}{A \code{GRangesList} of 5' UTR regions by transcript.}
#'   \item{three_utr}{A \code{GRangesList} of 3' UTR regions by transcript.}
#'   \item{exons}{A \code{GRangesList} of all exons by transcript.}
#'   \item{splice_sites}{A \code{GRanges} of splice site positions
#'     (donor +1,+2 and acceptor -1,-2).}
#'   \item{cds_start}{Named integer vector of CDS start positions
#'     (genomic coordinate of ATG).}
#'   \item{cds_end}{Named integer vector of CDS end positions
#'     (genomic coordinate of stop codon last base).}
#'   \item{strand}{Named character vector of strand ("+" or "-").}
#' }
#'
#' @importFrom GenomicFeatures cdsBy fiveUTRsByTranscript threeUTRsByTranscript exonsBy isActiveSeq
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors metadata mcols
#' @export
#'
#' @examples
#' \donttest{
#' library(EnsDb.Hsapiens.v86)
#' mane <- fetch_mane_txdb(EnsDb.Hsapiens.v86)
#' struct <- get_transcript_structure(mane, transcript_ids = NULL)
#' }
get_transcript_structure <- function(txdb, transcript_ids = NULL) {

    # If transcript_ids provided, restrict the TxDb
    if (!is.null(transcript_ids)) {
        if (inherits(txdb, "EnsDb")) {
            tx_filter <- ensembldb::TxIdFilter(transcript_ids)
            txdb <- ensembldb::filter(txdb, tx_filter)
        } else {
            # For TxDb, use isActiveSeq mechanism
            all_tx <- GenomicFeatures::transcripts(txdb)
            active <- tx_id %in% transcript_ids
            # Work with GRangesList subsetting instead
        }
    }

    # Extract regions by transcript
    cds_all <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
    five_utr_all <- GenomicFeatures::fiveUTRsByTranscript(txdb,
                                                          use.names = TRUE)
    three_utr_all <- GenomicFeatures::threeUTRsByTranscript(txdb,
                                                            use.names = TRUE)
    exons_all <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

    # Filter to transcripts that actually have CDS
    tx_with_cds <- intersect(names(cds_all), names(exons_all))
    if (length(tx_with_cds) == 0) {
        stop("No transcripts with CDS found in the provided TxDb.")
    }

    cds_all <- cds_all[tx_with_cds]
    exons_all <- exons_all[tx_with_cds]
    five_utr_all <- five_utr_all[intersect(names(five_utr_all), tx_with_cds)]
    three_utr_all <- three_utr_all[intersect(names(three_utr_all),
                                             tx_with_cds)]

    # Compute CDS start/end per transcript (genomic coordinates, 1-based)
    cds_starts <- vapply(cds_all, function(gr) {
        if (as.character(BiocGenerics::strand(gr)[1]) == "+") {
            min(BiocGenerics::start(gr))
        } else {
            max(BiocGenerics::end(gr))
        }
    }, numeric(1))

    cds_ends <- vapply(cds_all, function(gr) {
        if (as.character(BiocGenerics::strand(gr)[1]) == "+") {
            max(BiocGenerics::end(gr))
        } else {
            min(BiocGenerics::start(gr))
        }
    }, numeric(1))

    strands <- vapply(cds_all, function(gr) {
        as.character(BiocGenerics::strand(gr)[1])
    }, character(1))

    # Compute canonical splice sites (±1, ±2 from exon boundaries)
    splice_sites <- .compute_splice_sites(exons_all, cds_all, strands)

    structure(
        list(
            cds          = cds_all,
            five_utr     = five_utr_all,
            three_utr    = three_utr_all,
            exons        = exons_all,
            splice_sites = splice_sites,
            cds_start    = cds_starts,
            cds_end      = cds_ends,
            strand       = strands
        ),
        class = "transcript_structure"
    )
}


#' Compute canonical splice site positions
#'
#' Identifies ±1 and ±2 intronic positions at each exon-intron
#' boundary that lies within or adjacent to the CDS.
#'
#' @param exons GRangesList of exons by transcript.
#' @param cds GRangesList of CDS regions by transcript.
#' @param strands Named character vector of strand per transcript.
#'
#' @return A GRanges of splice site positions with mcols:
#'   \code{transcript_id}, \code{site_type} ("donor" or "acceptor"),
#'   \code{exon_boundary} (nearest coding coordinate),
#'   \code{offset} (+1, +2, -1, -2).
#'
#' @keywords internal
.compute_splice_sites <- function(exons, cds, strands) {

    all_sites <- GenomicRanges::GRanges()

    for (tx_id in names(exons)) {
        tx_exons <- GenomicRanges::reduce(exons[[tx_id]])
        tx_cds <- GenomicRanges::reduce(cds[[tx_id]])
        tx_strand <- strands[tx_id]

        if (length(tx_exons) < 2) next  # Need at least 2 exons for splicing

        # Sort exons by genomic position, then by strand direction
        tx_exons <- sort(tx_exons)
        if (tx_strand == "-") {
            tx_exons <- rev(tx_exons)
        }

        for (i in seq_len(length(tx_exons) - 1)) {
            exon_a <- tx_exons[i]
            exon_b <- tx_exons[i + 1]

            # Donor: first 2 bases of intron (after exon_a)
            # Acceptor: last 2 bases of intron (before exon_b)
            if (tx_strand == "+") {
                donor_start <- BiocGenerics::end(exon_a) + 1L
                donor_end   <- BiocGenerics::end(exon_a) + 2L
                acceptor_end   <- BiocGenerics::start(exon_b) - 1L
                acceptor_start <- BiocGenerics::start(exon_b) - 2L
            } else {
                donor_start <- BiocGenerics::start(exon_a) - 2L
                donor_end   <- BiocGenerics::start(exon_a) - 1L
                acceptor_end   <- BiocGenerics::end(exon_b) + 2L
                acceptor_start <- BiocGenerics::end(exon_b) + 1L
            }

            seqname <- as.character(GenomeInfoDb::seqnames(exon_a))

            # Create donor sites (+1, +2)
            if (donor_start > 0 && donor_end > 0) {
                donor_gr <- GenomicRanges::GRanges(
                    seqnames = seqname,
                    ranges = IRanges::IRanges(
                        start = donor_start,
                        end   = donor_end
                    ),
                    strand = tx_strand,
                    transcript_id = tx_id,
                    site_type     = "donor",
                    exon_idx      = i,
                    offset        = c(1L, 2L)
                )
                all_sites <- c(all_sites, donor_gr)
            }

            # Create acceptor sites (-1, -2)
            if (acceptor_start > 0 && acceptor_end > 0) {
                acceptor_gr <- GenomicRanges::GRanges(
                    seqnames = seqname,
                    ranges = IRanges::IRanges(
                        start = acceptor_start,
                        end   = acceptor_end
                    ),
                    strand = tx_strand,
                    transcript_id = tx_id,
                    site_type     = "acceptor",
                    exon_idx      = i + 1L,
                    offset        = c(-2L, -1L)
                )
                all_sites <- c(all_sites, acceptor_gr)
            }
        }
    }

    all_sites
}
