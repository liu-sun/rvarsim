#' Format variants in HGVS notation
#'
#' Converts a variants data.frame (from \code{\link{generate_variants}})
#' into HGVS-compliant coding DNA notation.
#'
#' HGVS conventions used:
#' \itemize{
#'   \item CDS: \code{NM_.*:c.[pos][ref]>[alt]} (e.g., c.215C>G)
#'   \item 5' UTR: \code{c.-[pos][ref]>[alt]} (e.g., c.-14G>A)
#'   \item 3' UTR: \code{c.*[pos][ref]>[alt]} (e.g., c.*32T>C)
#'   \item Splice donor: \code{c.[boundary]+[offset][ref]>[alt]}
#'     (e.g., c.453+1G>T)
#'   \item Splice acceptor: \code{c.[boundary]-[offset][ref]>[alt]}
#'     (e.g., c.612-2A>G)
#' }
#'
#' @param variants A \code{data.frame} from
#'   \code{\link{generate_variants}}.
#' @param use_transcript_version Logical. If \code{TRUE}, appends
#'   transcript version numbers (e.g., \code{NM_000546.6}) when the
#'   transcript ID is a RefSeq accession. Default: \code{FALSE}.
#'
#' @return A \code{data.frame} with all original columns plus:
#' \describe{
#'   \item{hgvs_c}{HGVS coding DNA notation string.}
#'   \item{hgvs_g}{HGVS genomic notation string
#'     (e.g., \code{NC_000001.11:g.123456A>G}).}
#' }
#'
#' @export
#'
#' @examples
#' # Self-contained format demo:
#' df <- data.frame(
#'     transcript_id = "NM_TEST",
#'     region = c("cds", "five_utr"),
#'     genomic_pos = c(215L, 30L),
#'     genomic_ref = c("C", "G"),
#'     genomic_alt = c("G", "A"),
#'     cds_pos = c(215L, -14L),
#'     exon_boundary = c(NA_integer_, NA_integer_),
#'     offset = c(NA_integer_, NA_integer_),
#'     stringsAsFactors = FALSE
#' )
#' format_hgvs(df)
format_hgvs <- function(variants, use_transcript_version = FALSE) {

    if (nrow(variants) == 0) {
        variants$hgvs_c <- character(0)
        variants$hgvs_g <- character(0)
        return(variants)
    }

    # Build HGVS c. notation for each variant
    variants$hgvs_c <- vapply(seq_len(nrow(variants)), function(i) {
        .build_hgvs_c(
            region      = variants$region[i],
            cds_pos     = variants$cds_pos[i],
            exon_boundary = variants$exon_boundary[i],
            offset      = variants$offset[i],
            ref         = variants$genomic_ref[i],
            alt         = variants$genomic_alt[i],
            transcript  = variants$transcript_id[i],
            use_version = use_transcript_version
        )
    }, character(1))

    # Build HGVS g. notation
    variants$hgvs_g <- vapply(seq_len(nrow(variants)), function(i) {
        .build_hgvs_g(
            chr         = as.character(variants$genomic_pos[i]),
            pos         = variants$genomic_pos[i],
            ref         = variants$genomic_ref[i],
            alt         = variants$genomic_alt[i]
        )
    }, character(1))

    variants
}


#' Build HGVS c. notation for a single variant
#'
#' @param region Region type string.
#' @param cds_pos CDS-relative position.
#' @param exon_boundary Exon boundary coordinate for splice sites.
#' @param offset Intronic offset for splice sites.
#' @param ref Reference allele.
#' @param alt Alternative allele.
#' @param transcript Transcript ID.
#' @param use_version Whether to include version number.
#'
#' @return HGVS c. notation string.
#'
#' @keywords internal
.build_hgvs_c <- function(region, cds_pos, exon_boundary, offset,
                           ref, alt, transcript,
                           use_version = FALSE) {

    tx_label <- if (use_version && grepl("^NM_", transcript)) {
        transcript  # version already in the name if from RefSeq
    } else {
        transcript
    }

    if (is.na(ref) || is.na(alt) || ref == alt) {
        return(NA_character_)
    }

    switch(
        region,
        cds = {
            if (is.na(cds_pos)) return(NA_character_)
            sprintf("%s:c.%d%s>%s", tx_label, cds_pos, ref, alt)
        },
        five_utr = {
            if (is.na(cds_pos)) return(NA_character_)
            # cds_pos is negative or zero
            sprintf("%s:c.%d%s>%s", tx_label, cds_pos, ref, alt)
        },
        three_utr = {
            if (is.na(cds_pos)) return(NA_character_)
            # cds_pos is positive (distance after stop codon)
            sprintf("%s:c.*%d%s>%s", tx_label, cds_pos, ref, alt)
        },
        splice_donor = {
            if (is.na(exon_boundary) || is.na(offset))
                return(NA_character_)
            sprintf("%s:c.%d+%d%s>%s", tx_label,
                    exon_boundary, abs(offset), ref, alt)
        },
        splice_acceptor = {
            if (is.na(exon_boundary) || is.na(offset))
                return(NA_character_)
            sprintf("%s:c.%d%d%s>%s", tx_label,
                    exon_boundary, offset, ref, alt)
        },
        NA_character_
    )
}


#' Build HGVS genomic (g.) notation for a single variant
#'
#' @param chr Chromosome name.
#' @param pos Genomic position (1-based).
#' @param ref Reference allele.
#' @param alt Alternative allele.
#'
#' @return HGVS g. notation string.
#'
#' @keywords internal
.build_hgvs_g <- function(chr, pos, ref, alt) {
    if (is.na(ref) || is.na(alt) || ref == alt) {
        return(NA_character_)
    }
    sprintf("%s:g.%d%s>%s", chr, pos, ref, alt)
}
