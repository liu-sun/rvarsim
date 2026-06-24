#' Generate all possible single nucleotide variants
#'
#' For each position across CDS, UTR, and splice site regions, extracts
#' the reference base from the genome and generates the three possible
#' alternative alleles (A, C, G, T minus the reference).
#'
#' @param tx_struct A \code{transcript_structure} object from
#'   \code{\link{get_transcript_structure}}.
#' @param bsgenome A \code{BSgenome} object (e.g.,
#'   \code{BSgenome.Hsapiens.UCSC.hg38}) providing the reference
#'   genome sequence.
#' @param regions Character vector specifying which region types to
#'   include. Options: \code{"cds"}, \code{"five_utr"},
#'   \code{"three_utr"}, \code{"splice_site"}. Default: all four.
#' @param transcript_ids Optional character vector of transcript IDs
#'   to limit generation. \code{NULL} means all transcripts in the
#'   structure.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{transcript_id}{Transcript identifier.}
#'   \item{region}{Region type: "cds", "five_utr", "three_utr",
#'     "splice_donor", "splice_acceptor".}
#'   \item{genomic_pos}{Genomic position (1-based).}
#'   \item{genomic_ref}{Reference allele at this position.}
#'   \item{genomic_alt}{Alternative allele.}
#'   \item{cds_pos}{Position relative to CDS start
#'     (1 = A of ATG; negative = 5'UTR; NA for splice sites
#'     described by boundary).}
#'   \item{exon_boundary}{For splice sites: the nearest coding
#'     coordinate (NA for other regions).}
#'   \item{offset}{For splice sites: the intronic offset
#'     (\code{+1}, \code{+2}, \code{-1}, \code{-2}; NA otherwise).}
#' }
#'
#' @importFrom BSgenome getSeq
#' @importFrom GenomeInfoDb seqnames seqinfo
#' @importFrom IRanges width
#' @export
#'
#' @examples
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     cat("Generate variants with BSgenome.Hsapiens.UCSC.hg38")
#' }
generate_variants <- function(tx_struct, bsgenome,
                               regions = c("cds", "five_utr",
                                           "three_utr",
                                           "splice_site"),
                               transcript_ids = NULL) {

    valid_regions <- c("cds", "five_utr", "three_utr", "splice_site")
    regions <- match.arg(regions, valid_regions, several.ok = TRUE)

    all_variants <- vector("list", 0)

    # Determine which transcripts to process
    tx_ids <- names(tx_struct$cds)
    if (!is.null(transcript_ids)) {
        tx_ids <- intersect(tx_ids, transcript_ids)
    }

    # Validate that the genome has the expected chromosomes
    seqlevels_ok <- .validate_seqlevels(tx_struct, bsgenome)
    if (!seqlevels_ok) {
        warning("Some chromosome names in the TxDb may not match ",
                "the BSgenome. Variants on unmatched chromosomes ",
                "will be skipped.")
    }

    for (tx_id in tx_ids) {
        tx_vars <- list()

        # --- CDS variants ---
        if ("cds" %in% regions && tx_id %in% names(tx_struct$cds)) {
            tx_vars$cds <- .generate_region_variants(
                tx_struct$cds[[tx_id]], bsgenome, tx_id,
                tx_struct$cds_start[tx_id],
                tx_struct$cds_end[tx_id],
                tx_struct$strand[tx_id],
                region_type = "cds"
            )
        }

        # --- 5' UTR variants ---
        if ("five_utr" %in% regions && tx_id %in% names(tx_struct$five_utr)) {
            tx_vars$five_utr <- .generate_region_variants(
                tx_struct$five_utr[[tx_id]], bsgenome, tx_id,
                tx_struct$cds_start[tx_id],
                tx_struct$cds_end[tx_id],
                tx_struct$strand[tx_id],
                region_type = "five_utr"
            )
        }

        # --- 3' UTR variants ---
        if ("three_utr" %in% regions &&
            tx_id %in% names(tx_struct$three_utr)) {
            tx_vars$three_utr <- .generate_region_variants(
                tx_struct$three_utr[[tx_id]], bsgenome, tx_id,
                tx_struct$cds_start[tx_id],
                tx_struct$cds_end[tx_id],
                tx_struct$strand[tx_id],
                region_type = "three_utr"
            )
        }

        # --- Splice site variants ---
        if ("splice_site" %in% regions) {
            tx_vars$splice <- .generate_splice_variants(
                tx_struct$splice_sites, bsgenome, tx_id,
                tx_struct$cds_start[tx_id],
                tx_struct$strand[tx_id]
            )
        }

        # Combine variants for this transcript
        if (length(tx_vars) > 0) {
            all_variants[[tx_id]] <- do.call(rbind, tx_vars)
        }
    }

    result <- do.call(rbind, all_variants)
    if (is.null(result)) {
        result <- data.frame(
            transcript_id   = character(0),
            region          = character(0),
            genomic_pos     = integer(0),
            genomic_ref     = character(0),
            genomic_alt     = character(0),
            cds_pos         = integer(0),
            exon_boundary   = integer(0),
            offset          = integer(0),
            stringsAsFactors = FALSE
        )
    }

    rownames(result) <- NULL
    result
}


#' Generate variants for a single region
#'
#' @param gr A GRanges of region positions.
#' @param bsgenome BSgenome object.
#' @param tx_id Transcript ID.
#' @param cds_start CDS start genomic position.
#' @param cds_end CDS end genomic position.
#' @param strand Strand.
#' @param region_type One of "cds", "five_utr", "three_utr".
#'
#' @return data.frame of variants.
#'
#' @keywords internal
.generate_region_variants <- function(gr, bsgenome, tx_id,
                                       cds_start, cds_end, strand,
                                       region_type) {

    # Expand each range into individual positions
    positions <- .range_to_positions(gr)

    if (length(positions) == 0) return(NULL)

    # Get reference sequence for each position
    ref_bases <- .get_ref_bases(gr, bsgenome)

    # Map position index back to genomic coordinate
    genomic_pos <- .range_to_positions(gr, as_genomic = TRUE)

    # Build data.frame: for each position, 3 alternative alleles
    result_list <- lapply(seq_along(ref_bases), function(i) {
        ref <- ref_bases[i]
        alts <- setdiff(c("A", "C", "G", "T"), ref)

        gpos <- genomic_pos[i]
        cpos <- .compute_cds_pos(gpos, cds_start, cds_end, strand,
                                 region_type)

        data.frame(
            transcript_id = tx_id,
            region        = region_type,
            genomic_pos   = gpos,
            genomic_ref   = ref,
            genomic_alt   = alts,
            cds_pos       = cpos,
            exon_boundary = NA_integer_,
            offset        = NA_integer_,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, result_list)
}


#' Generate splice site variants
#'
#' @param splice_gr GRanges of splice site positions from
#'   transcript_structure.
#' @param bsgenome BSgenome object.
#' @param tx_id Transcript ID to filter.
#' @param cds_start CDS start genomic position.
#' @param strand Strand.
#'
#' @return data.frame of splice site variants.
#'
#' @keywords internal
.generate_splice_variants <- function(splice_gr, bsgenome, tx_id,
                                       cds_start, strand) {

    tx_sites <- splice_gr[splice_gr$transcript_id == tx_id]
    if (length(tx_sites) == 0) return(NULL)

    # For each splice site position, get ref base and generate alts
    result_list <- lapply(seq_along(tx_sites), function(i) {
        site <- tx_sites[i]
        gpos <- BiocGenerics::start(site)

        # Get reference base
        seqname <- as.character(GenomeInfoDb::seqnames(site))
        seqlevels_style <- GenomeInfoDb::seqlevelsStyle(bsgenome)
        GenomeInfoDb::seqlevelsStyle(seqname) <- seqlevels_style[1]

        ref <- tryCatch({
            s <- BSgenome::getSeq(bsgenome, seqname,
                                  start = gpos, end = gpos)
            as.character(s)
        }, error = function(e) NA_character_)

        if (is.na(ref)) return(NULL)

        alts <- setdiff(c("A", "C", "G", "T"), ref)

        # Compute exon boundary coordinate (nearest coding position
        # at the exon edge)
        boundary <- .compute_exon_boundary(
            gpos, site$site_type, site$exon_idx,
            site$offset, cds_start, strand,
            tx_id
        )

        data.frame(
            transcript_id = tx_id,
            region        = paste0("splice_", site$site_type),
            genomic_pos   = gpos,
            genomic_ref   = ref,
            genomic_alt   = alts,
            cds_pos       = NA_integer_,
            exon_boundary = boundary,
            offset        = site$offset,
            stringsAsFactors = FALSE
        )
    })

    result_list <- result_list[!vapply(result_list, is.null, logical(1))]
    if (length(result_list) == 0) return(NULL)
    do.call(rbind, result_list)
}


#' Expand GRanges to a vector of individual positions
#'
#' @param gr A GRanges object.
#' @param as_genomic If TRUE, return genomic coordinates.
#'   If FALSE, return sequential indices.
#'
#' @return Integer vector of positions.
#'
#' @keywords internal
.range_to_positions <- function(gr, as_genomic = TRUE) {
    pos_list <- lapply(seq_along(gr), function(i) {
        rng <- gr[i]
        if (as_genomic) {
            seq(BiocGenerics::start(rng), BiocGenerics::end(rng))
        } else {
            seq_len(IRanges::width(rng))
        }
    })
    unlist(pos_list)
}


#' Get reference bases for each position in a GRanges
#'
#' Expands the GRanges to single-base positions and retrieves the
#' reference allele from the BSgenome for each.
#'
#' @param gr A GRanges object.
#' @param bsgenome A BSgenome object.
#'
#' @return Character vector of reference bases (same length as the
#'   expanded positions).
#'
#' @keywords internal
.get_ref_bases <- function(gr, bsgenome) {
    # For each range, get the sequence from the BSgenome
    bases <- character(0)
    for (i in seq_along(gr)) {
        rng <- gr[i]
        seqname <- as.character(GenomeInfoDb::seqnames(rng))
        seq_bases <- tryCatch({
            GenomeInfoDb::seqlevelsStyle(seqname) <-
                GenomeInfoDb::seqlevelsStyle(bsgenome)[1]
            s <- BSgenome::getSeq(bsgenome, seqname,
                                  start = BiocGenerics::start(rng),
                                  end   = BiocGenerics::end(rng))
            strsplit(as.character(s), "")[[1]]
        }, error = function(e) {
            rep(NA_character_, IRanges::width(rng))
        })
        bases <- c(bases, seq_bases)
    }
    bases
}


#' Compute CDS-relative position
#'
#' Maps a genomic coordinate to a CDS-relative position following
#' HGVS numbering conventions.
#'
#' @param gpos Genomic position (1-based).
#' @param cds_start Genomic coordinate of CDS start (ATG).
#' @param cds_end Genomic coordinate of CDS end (stop codon last base).
#' @param strand Strand ("+" or "-").
#' @param region_type Region type ("cds", "five_utr", "three_utr").
#'
#' @return Integer: positive for CDS, negative for 5'UTR,
#'   NA for 3'UTR (stored separately, * notation added by formatter).
#'
#' @keywords internal
.compute_cds_pos <- function(gpos, cds_start, cds_end, strand,
                              region_type) {
    if (is.na(cds_start) || is.na(cds_end)) return(NA_integer_)

    if (strand == "+") {
        if (region_type == "cds") {
            gpos - cds_start + 1L
        } else if (region_type == "five_utr") {
            gpos - cds_start  # negative or zero
        } else if (region_type == "three_utr") {
            gpos - cds_end     # positive, *-notation
        } else {
            NA_integer_
        }
    } else {
        if (region_type == "cds") {
            cds_start - gpos + 1L
        } else if (region_type == "five_utr") {
            cds_start - gpos   # negative or zero for 5'UTR
        } else if (region_type == "three_utr") {
            cds_end - gpos      # negative for 3'UTR on minus strand
        } else {
            NA_integer_
        }
    }
}


#' Compute exon boundary coordinate for splice site HGVS notation
#'
#' Determines the nearest coding coordinate at the exon-intron
#' boundary for building HGVS c. notation (e.g., c.453+1G>T).
#'
#' @param gpos Genomic position of the splice site.
#' @param site_type "donor" or "acceptor".
#' @param exon_idx Exon index in transcript.
#' @param offset Intronic offset (+1, +2, -1, -2).
#' @param cds_start CDS start genomic coordinate.
#' @param strand Strand.
#' @param tx_id Transcript ID (for error context).
#'
#' @return Integer: the nearest coding coordinate, or NA.
#'
#' @keywords internal
.compute_exon_boundary <- function(gpos, site_type, exon_idx, offset,
                                    cds_start, strand, tx_id) {
    if (site_type == "donor") {
        if (strand == "+") {
            # Donor: intron starts at exon_end + 1
            # gpos = exon_end + |offset|
            # exon_end = gpos - 1
            exon_end_in_genomic <- gpos - abs(offset)
            .compute_cds_pos(exon_end_in_genomic, cds_start,
                             NA_real_, strand, "cds")
        } else {
            # Donor on minus strand: intron starts at exon_start - 1
            exon_start_in_genomic <- gpos + abs(offset)
            .compute_cds_pos(exon_start_in_genomic, cds_start,
                             NA_real_, strand, "cds")
        }
    } else {
        # acceptor
        if (strand == "+") {
            # Acceptor: intron ends at exon_start - 1
            exon_start_in_genomic <- gpos + abs(offset) + 1L
            .compute_cds_pos(exon_start_in_genomic, cds_start,
                             NA_real_, strand, "cds")
        } else {
            # Acceptor on minus strand: intron ends at exon_end + 1
            exon_end_in_genomic <- gpos - abs(offset) - 1L
            .compute_cds_pos(exon_end_in_genomic, cds_start,
                             NA_real_, strand, "cds")
        }
    }
}


#' Validate that chromosome names in the structure exist in the
#' BSgenome.
#'
#' @param tx_struct Transcript structure object.
#' @param bsgenome BSgenome object.
#'
#' @return Logical: TRUE if all seqlevels are valid.
#'
#' @keywords internal
.validate_seqlevels <- function(tx_struct, bsgenome) {
    tx_seqs <- unique(c(
        as.character(GenomeInfoDb::seqnames(
            unlist(tx_struct$cds, use.names = FALSE))),
        as.character(GenomeInfoDb::seqnames(tx_struct$splice_sites))
    ))
    bs_seqs <- GenomeInfoDb::seqnames(bsgenome)
    all(tx_seqs %in% bs_seqs)
}
