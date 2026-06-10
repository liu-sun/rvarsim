#' Map variants between genome and transcript notation
#'
#' Converts HGVS variants between genomic (g.) and coding (c.) notation
#' using a TxDb/EnsDb to resolve exon/CDS structure. Handles intronic
#' offsets, UTR positions, and splice boundary coordinates.
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param txdb A \code{TxDb} or \code{EnsDb} object providing
#'   transcript structure.
#' @param bsgenome A \code{BSgenome} object for sequence context.
#' @param transcript_id Optional transcript ID to use when
#'   disambiguating c. notation (required when multiple transcripts
#'   share the same accession).
#' @param direction Direction of transcription: \code{"c_to_g"},
#'   \code{"g_to_c"}, or \code{"auto"} (auto-detect from notation).
#'   Default: \code{"auto"}.
#'
#' @return A list of HGVS variant strings in the target notation.
#'
#' @export
#'
#' @examples
#' # Transcribe needs a TxDb with matching BSgenome seqlevels:
#' if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     cat("Ready for transcription mapping")
#' }
transcribe_hgvs <- function(hgvs_strings, txdb, bsgenome,
                              transcript_id = NULL,
                              direction = c("c_to_g",
                                            "g_to_c",
                                            "auto")) {
    direction <- match.arg(direction)

    lapply(hgvs_strings, function(s) {
        parsed <- tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
        if (is.na(parsed[[1]])) return(NA_character_)

        detect <- if (direction == "auto") {
            if (parsed$notation == "c") "c_to_g" else "g_to_c"
        } else direction

        switch(detect,
               c_to_g = .transcribe_c_to_g(parsed, txdb, bsgenome,
                                            transcript_id),
               g_to_c = .transcribe_g_to_c(parsed, txdb, bsgenome,
                                            transcript_id),
               NA_character_)
    })
}


#' @rdname transcribe_hgvs
#' @export
c_to_g <- function(hgvs_strings, txdb, bsgenome,
                    transcript_id = NULL) {
    transcribe_hgvs(hgvs_strings, txdb, bsgenome, transcript_id,
                     direction = "c_to_g")
}


#' @rdname transcribe_hgvs
#' @export
g_to_c <- function(hgvs_strings, txdb, bsgenome,
                    transcript_id = NULL) {
    transcribe_hgvs(hgvs_strings, txdb, bsgenome, transcript_id,
                     direction = "g_to_c")
}


#' Map c. variant to g. notation
#'
#' Uses TxDb to resolve the CDS-relative position to a genomic
#' coordinate, considering strand, intronic offsets, and UTR
#' positions.
#'
#' @keywords internal
.transcribe_c_to_g <- function(parsed, txdb, bsgenome,
                                transcript_id) {
    # Get CDS structure for the transcript
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx",
                                         use.names = TRUE)

    # Find matching transcript
    tx_id <- .resolve_transcript(parsed$accession, cds_by_tx,
                                  transcript_id)
    if (is.null(tx_id) || !tx_id %in% names(cds_by_tx)) {
        warning("Transcript ", parsed$accession,
                " not found in TxDb. Cannot map.")
        return(NA_character_)
    }

    cds_gr <- cds_by_tx[[tx_id]]
    exon_gr <- GenomicFeatures::exonsBy(txdb, by = "tx",
                                         use.names = TRUE)[[tx_id]]

    tx_strand <- as.character(BiocGenerics::strand(cds_gr)[1])

    # Get CDS start and end in genomic coordinates
    if (tx_strand == "+") {
        cds_genomic_start <- min(BiocGenerics::start(cds_gr))
        cds_genomic_end   <- max(BiocGenerics::end(cds_gr))
    } else {
        cds_genomic_start <- max(BiocGenerics::end(cds_gr))
        cds_genomic_end   <- min(BiocGenerics::start(cds_gr))
    }

    pos <- parsed$position
    cds_pos <- pos$start

    # Compute genomic position based on strand and region
    if (pos$is_five_utr) {
        # 5'UTR: c.-14 → genomic position upstream of CDS start
        if (tx_strand == "+") {
            gpos <- cds_genomic_start + cds_pos
        } else {
            gpos <- cds_genomic_start - cds_pos
        }
    } else if (pos$is_three_utr) {
        # 3'UTR: c.*32 → genomic position downstream of CDS end
        if (tx_strand == "+") {
            gpos <- cds_genomic_end + cds_pos
        } else {
            gpos <- cds_genomic_end - cds_pos
        }
    } else if (pos$offset_start != 0) {
        # Intronic: c.453+1 → genomic position with intronic offset
        gpos <- .resolve_intronic_position(
            cds_pos, pos$offset_start, cds_gr, exon_gr, tx_strand
        )
    } else {
        # Standard CDS position
        if (tx_strand == "+") {
            gpos <- cds_genomic_start + cds_pos - 1L
        } else {
            gpos <- cds_genomic_start - cds_pos + 1L
        }
    }

    # Build g. notation
    chr <- as.character(GenomeInfoDb::seqnames(cds_gr)[1])
    hgvs_g <- sprintf("%s:g.%d%s>%s",
                       chr, gpos, parsed$reference, parsed$alternate)

    hgvs_g
}


#' Map g. variant to c. notation
#'
#' Uses TxDb to resolve a genomic position to CDS-relative coordinates.
#'
#' @keywords internal
.transcribe_g_to_c <- function(parsed, txdb, bsgenome,
                                transcript_id) {
    # Get CDS for candidate transcripts that overlap this genomic pos
    gpos <- parsed$position$start

    gr_point <- GenomicRanges::GRanges(
        seqnames = .extract_chromosome(parsed$accession),
        ranges   = IRanges::IRanges(start = gpos, end = gpos),
        strand   = "*"
    )

    # Find transcripts whose CDS overlaps this position
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx",
                                         use.names = TRUE)
    utr5_by_tx <- GenomicFeatures::fiveUTRsByTranscript(txdb,
                                                         use.names = TRUE)
    utr3_by_tx <- GenomicFeatures::threeUTRsByTranscript(txdb,
                                                          use.names = TRUE)
    exons_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx",
                                             use.names = TRUE)

    if (!is.null(transcript_id)) {
        tx_ids <- transcript_id
    } else {
        # Find all transcripts overlapping this position
        tx_hits <- unique(c(
            names(cds_by_tx)[IRanges::overlapsAny(cds_by_tx, gr_point)],
            names(utr5_by_tx)[IRanges::overlapsAny(utr5_by_tx, gr_point)],
            names(utr3_by_tx)[IRanges::overlapsAny(utr3_by_tx, gr_point)]
        ))
        tx_ids <- tx_hits
    }

    if (length(tx_ids) == 0) {
        warning("No transcript overlaps genomic position ", gpos)
        return(NA_character_)
    }

    # Use the first matching transcript
    tx_id <- tx_ids[1]
    cds_gr <- cds_by_tx[[tx_id]]
    tx_strand <- as.character(BiocGenerics::strand(cds_gr)[1])

    if (tx_strand == "+") {
        cds_start <- min(BiocGenerics::start(cds_gr))
        cds_end   <- max(BiocGenerics::end(cds_gr))
    } else {
        cds_start <- max(BiocGenerics::end(cds_gr))
        cds_end   <- min(BiocGenerics::start(cds_gr))
    }

    # Determine region and compute c. position
    is_utr5 <- FALSE
    is_utr3 <- FALSE
    cds_pos <- NA_integer_
    offset <- 0L

    # Check if position is in CDS
    cds_ranges <- GenomicRanges::reduce(cds_gr)
    if (IRanges::overlapsAny(gr_point, cds_ranges)) {
        if (tx_strand == "+") {
            cds_pos <- gpos - cds_start + 1L
        } else {
            cds_pos <- cds_start - gpos + 1L
        }
    } else if (tx_strand == "+") {
        if (gpos < cds_start) {
            # 5'UTR
            cds_pos <- gpos - cds_start
            is_utr5 <- TRUE
        } else if (gpos > cds_end) {
            # 3'UTR
            cds_pos <- gpos - cds_end
            is_utr3 <- TRUE
        }
    } else {
        if (gpos > cds_start) {
            # 5'UTR (upstream on minus strand = numerically higher)
            cds_pos <- cds_start - gpos
            is_utr5 <- TRUE
        } else if (gpos < cds_end) {
            # 3'UTR
            cds_pos <- cds_end - gpos
            is_utr3 <- TRUE
        }
    }

    # Also check intronic positions: between exons but within transcript
    if (is.na(cds_pos)) {
        result <- .compute_cds_from_intronic(
            gpos, parsed$accession, exons_by_tx, tx_id, tx_strand
        )
        if (!is.null(result)) {
            cds_pos <- result$cds_pos
            offset  <- result$offset
        }
    }

    if (is.na(cds_pos)) {
        warning("Cannot map genomic position ", gpos,
                " to a CDS coordinate for transcript ", tx_id)
        return(NA_character_)
    }

    # Build c. notation
    .build_c_notation(tx_id, cds_pos, offset, is_utr5, is_utr3,
                       parsed$reference, parsed$alternate)
}


#' Build a c. notation string from coordinates
#'
#' @keywords internal
.build_c_notation <- function(tx_id, cds_pos, offset,
                               is_utr5, is_utr3, ref, alt) {
    if (is_utr5) {
        sprintf("%s:c.%d%s>%s", tx_id, cds_pos, ref, alt)
    } else if (is_utr3) {
        sprintf("%s:c.*%d%s>%s", tx_id, cds_pos, ref, alt)
    } else if (offset > 0) {
        sprintf("%s:c.%d+%d%s>%s", tx_id, cds_pos, offset, ref, alt)
    } else if (offset < 0) {
        sprintf("%s:c.%d%d%s>%s", tx_id, cds_pos, offset, ref, alt)
    } else {
        sprintf("%s:c.%d%s>%s", tx_id, cds_pos, ref, alt)
    }
}


#' Resolve intronic CDS position to genomic coordinates
#'
#' Given a CDS exon boundary position and intronic offset,
#' compute the genomic position.
#'
#' @keywords internal
.resolve_intronic_position <- function(cds_pos, offset,
                                        cds_gr, exon_gr, strand) {
    # For a CDS position with intronic offset, find which exon
    # this position falls in/near and compute the genomic position
    if (strand == "+") {
        cds_start <- min(BiocGenerics::start(cds_gr))
        gpos <- cds_start + cds_pos - 1L + offset
    } else {
        cds_start <- max(BiocGenerics::end(cds_gr))
        gpos <- cds_start - cds_pos + 1L - offset
    }
    gpos
}


#' Compute CDS-relative position and offset for intronic genomic pos
#'
#' @keywords internal
.compute_cds_from_intronic <- function(gpos, chr, exons_by_tx,
                                        tx_id, strand) {
    # Check if gpos falls between exons (intronic)
    exons <- exons_by_tx[[tx_id]]
    exons <- sort(exons)

    # Find the nearest exon boundary
    for (i in seq_len(length(exons))) {
        exon <- exons[i]
        if (strand == "+") {
            # Check if position is just after this exon (donor)
            if (gpos > BiocGenerics::end(exon) &&
                (i == length(exons) ||
                 gpos < BiocGenerics::start(exons[i + 1]))) {
                # Intronic: compute offset from exon end
                exon_end_cds <- .compute_cds_pos(
                    BiocGenerics::end(exon),
                    min(BiocGenerics::start(exons)), NA_real_,
                    strand, "cds"
                )
                offset <- gpos - BiocGenerics::end(exon)
                return(list(cds_pos = exon_end_cds, offset = offset))
            }
        } else {
            # Minus strand
            if (gpos < BiocGenerics::start(exon) &&
                (i == 1 ||
                 gpos > BiocGenerics::end(exons[i - 1]))) {
                exon_start_genomic <- BiocGenerics::start(exon)
                exon_start_cds <- .compute_cds_pos(
                    exon_start_genomic,
                    max(BiocGenerics::end(exons)), NA_real_,
                    strand, "cds"
                )
                offset <- exon_start_genomic - gpos
                return(list(cds_pos = exon_start_cds,
                            offset = -offset))
            }
        }
    }

    NULL
}


#' Resolve a transcript accession to a transcript ID in the TxDb
#'
#' @keywords internal
.resolve_transcript <- function(accession, cds_by_tx, transcript_id) {
    if (!is.null(transcript_id)) return(transcript_id)

    # Try exact match
    if (accession %in% names(cds_by_tx)) return(accession)

    # Try matching transcript names (for RefSeq NM_ accessions)
    all_tx <- names(cds_by_tx)
    matches <- grep(accession, all_tx, fixed = TRUE)
    if (length(matches) > 0) return(all_tx[matches[1]])

    NULL
}


#' Extract chromosome identifier from an accession
#'
#' @keywords internal
.extract_chromosome <- function(accession) {
    # NC_000001.11 → 1, or just return the accession for chr naming
    if (grepl("^NC_", accession)) {
        chr_num <- as.integer(substr(accession, 4, 9)) - 1000L
        paste0("chr", chr_num)
    } else {
        accession
    }
}
