#' Convert HGVS to and from other variant formats
#'
#' Converts between HGVS notation and common variant representation
#' formats: VCF (Variant Call Format), SPDI (Sequence Position
#' Deletion Insertion), and VRS (Variation Representation
#' Specification, basic support).
#'
#' @section VCF format:
#' VCF uses 1-based positions: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
#' Conversion from HGVS c. notation requires transcript-to-genomic mapping.
#'
#' @section SPDI format:
#' NCBI's SPDI uses 0-based interbase coordinates:
#' sequence_id:position:deleted_sequence:inserted_sequence
#'
#' @return A \code{data.frame} or character vector, depending on the
#'   conversion function used.
#'
#' @name hgvs_conversion
NULL


#' Convert HGVS to VCF data.frame
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param txdb Optional \code{TxDb} for mapping c. to genomic
#'   coordinates. Required for transcript-level variants.
#' @param bsgenome Optional \code{BSgenome} for chromosome naming.
#'
#' @return A \code{data.frame} with VCF columns: CHROM, POS, ID, REF,
#'   ALT, QUAL, FILTER, INFO.
#'
#' @export
#'
#' @examples
#' # VCF conversion:
#' vcf <- hgvs_to_vcf("NC_000001.11:g.123456A>G")
#'
#' # SPDI round-trip:
#' hgvs_to_spdi("NC_000001.11:g.123456A>G")
hgvs_to_vcf <- function(hgvs_strings, txdb = NULL, bsgenome = NULL) {
    result <- lapply(hgvs_strings, function(s) {
        parsed <- tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
        if (is.na(parsed[[1]])) {
            return(data.frame(
                CHROM = NA_character_, POS = NA_integer_,
                ID = ".", REF = NA_character_, ALT = NA_character_,
                QUAL = ".", FILTER = ".", INFO = ".",
                stringsAsFactors = FALSE
            ))
        }

        .hgvs_to_vcf_single(parsed, txdb, bsgenome)
    })
    do.call(rbind, result)
}


#' Convert a single parsed HGVS to VCF
#'
#' @keywords internal
.hgvs_to_vcf_single <- function(parsed, txdb, bsgenome) {
    # If g. notation, extract chromosome and position directly
    if (parsed$notation == "g") {
        chr <- sub("\\..*$", "", parsed$accession)
        pos <- parsed$position$start
        ref <- if (is.na(parsed$reference) || nchar(parsed$reference) == 0)
            "" else parsed$reference
        alt <- if (is.na(parsed$alternate) || nchar(parsed$alternate) == 0)
            "" else parsed$alternate

        # VCF requires REF to be non-empty; for insertions, include
        # the preceding base
        if (ref == "" && alt != "") {
            # Need upstream reference context
        }
    } else if (parsed$notation == "c") {
        # Map c. to g. via TxDb
        if (!is.null(txdb) && !is.null(bsgenome)) {
            genomic <- .cds_to_genomic(parsed, txdb, bsgenome)
            if (!is.null(genomic)) {
                return(genomic)
            }
        }
        return(data.frame(
            CHROM = NA_character_, POS = NA_integer_,
            ID = ".", REF = parsed$reference, ALT = parsed$alternate,
            QUAL = ".", FILTER = ".", INFO = ".",
            stringsAsFactors = FALSE
        ))
    }

    data.frame(
        CHROM = if (exists("chr")) chr else NA_character_,
        POS   = if (exists("pos")) pos else NA_integer_,
        ID    = ".",
        REF   = if (exists("ref")) ref else parsed$reference,
        ALT   = if (exists("alt")) alt else parsed$alternate,
        QUAL  = ".",
        FILTER = ".",
        INFO  = paste0("HGVS=", parsed$raw),
        stringsAsFactors = FALSE
    )
}


#' Convert VCF row to HGVS notation
#'
#' @param vcf_df A \code{data.frame} with VCF columns: CHROM, POS, REF,
#'   ALT. May also contain ID, QUAL, FILTER.
#' @param assembly Genome assembly name for creating NC_ accession
#'   (e.g., "GRCh38").
#'
#' @return Character vector of HGVS g. notation strings.
#'
#' @export
#'
#' @examples
#' vcf <- data.frame(CHROM = "1", POS = 123456L, REF = "A", ALT = "G")
#' vcf_to_hgvs(vcf, "GRCh38")
vcf_to_hgvs <- function(vcf_df, assembly = "GRCh38") {
    vapply(seq_len(nrow(vcf_df)), function(i) {
        chr <- vcf_df$CHROM[i]
        pos <- vcf_df$POS[i]
        ref <- vcf_df$REF[i]
        alt <- vcf_df$ALT[i]

        # Determine variant type from REF/ALT
        if (nchar(ref) == 1 && nchar(alt) == 1) {
            # SNV
            sprintf("%s:g.%d%s>%s",
                    .get_chr_accession(chr, assembly), pos, ref, alt)
        } else if (nchar(ref) > nchar(alt) && alt == "") {
            # Deletion
            end_pos <- pos + nchar(ref) - 1L
            if (end_pos == pos) {
                sprintf("%s:g.%ddel%s",
                        .get_chr_accession(chr, assembly), pos, ref)
            } else {
                sprintf("%s:g.%d_%ddel%s",
                        .get_chr_accession(chr, assembly),
                        pos, end_pos, ref)
            }
        } else if (nchar(ref) < nchar(alt) && ref == "") {
            # Insertion
            sprintf("%s:g.%d_%dins%s",
                    .get_chr_accession(chr, assembly),
                    pos - 1L, pos, alt)
        } else if (nchar(ref) > 0 && nchar(alt) > 0) {
            # delins
            end_pos <- pos + nchar(ref) - 1L
            if (end_pos == pos) {
                sprintf("%s:g.%ddelins%s",
                        .get_chr_accession(chr, assembly), pos, alt)
            } else {
                sprintf("%s:g.%d_%ddelins%s",
                        .get_chr_accession(chr, assembly),
                        pos, end_pos, alt)
            }
        } else {
            NA_character_
        }
    }, character(1), USE.NAMES = FALSE)
}


#' Get standard chromosome accession from assembly
#'
#' @keywords internal
.get_chr_accession <- function(chr, assembly) {
    chr_num <- gsub("chr", "", chr, ignore.case = TRUE)
    if (assembly == "GRCh38" || assembly == "hg38") {
        sprintf("NC_%06d.11", as.integer(chr_num) + 1000L)
    } else if (assembly == "GRCh37" || assembly == "hg19") {
        sprintf("NC_%06d.10", as.integer(chr_num) + 1000L - 1L)
    } else {
        # Fallback
        paste0("NC_", chr)
    }
}


#' Convert HGVS to SPDI format
#'
#' SPDI uses 0-based interbase coordinates:
#' seqid:position:deletion:insertion
#'
#' @param hgvs_strings Character vector of HGVS g. notation strings.
#' @param assembly Genome assembly name.
#'
#' @return Character vector of SPDI strings.
#'
#' @export
#'
#' @examples
#' hgvs_to_spdi("NC_000001.11:g.123456A>G")
hgvs_to_spdi <- function(hgvs_strings, assembly = "GRCh38") {
    vapply(hgvs_strings, function(s) {
        parsed <- tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
        if (is.na(parsed[[1]])) return(NA_character_)

        # SPDI uses 0-based interbase: position = 1-based HGVS pos - 1
        spdi_pos <- parsed$position$start - 1L
        seqid <- parsed$accession

        deletion <- if (is.na(parsed$reference) ||
                        nchar(parsed$reference) == 0) {
            ""
        } else {
            parsed$reference
        }

        insertion <- if (is.na(parsed$alternate) ||
                         nchar(parsed$alternate) == 0) {
            ""
        } else {
            parsed$alternate
        }

        sprintf("%s:%d:%s:%s", seqid, spdi_pos, deletion, insertion)
    }, character(1), USE.NAMES = FALSE)
}


#' Convert SPDI to HGVS g. notation
#'
#' @param spdi_strings Character vector of SPDI strings.
#'
#' @return Character vector of HGVS g. notation strings.
#'
#' @examples
#' spdi_to_hgvs("NC_000001.11:123455:A:G")
#'
#' @export
spdi_to_hgvs <- function(spdi_strings) {
    vapply(spdi_strings, function(s) {
        parts <- strsplit(s, ":", fixed = TRUE)[[1]]
        if (length(parts) != 4) return(NA_character_)

        seqid    <- parts[1]
        position <- as.integer(parts[2])  # 0-based interbase
        deletion <- parts[3]
        insertion <- parts[4]

        # Convert to 1-based
        hgvs_pos <- position + 1L

        if (nchar(deletion) == 1 && nchar(insertion) == 1) {
            sprintf("%s:g.%d%s>%s", seqid, hgvs_pos, deletion, insertion)
        } else if (nchar(deletion) > 0 && nchar(insertion) == 0) {
            if (nchar(deletion) == 1) {
                sprintf("%s:g.%ddel%s", seqid, hgvs_pos, deletion)
            } else {
                end_pos <- hgvs_pos + nchar(deletion) - 1L
                sprintf("%s:g.%d_%ddel%s", seqid, hgvs_pos, end_pos,
                        deletion)
            }
        } else if (nchar(deletion) == 0 && nchar(insertion) > 0) {
            sprintf("%s:g.%d_%dins%s", seqid, hgvs_pos - 1L, hgvs_pos,
                    insertion)
        } else {
            if (nchar(deletion) == 1) {
                sprintf("%s:g.%ddelins%s", seqid, hgvs_pos, insertion)
            } else {
                end_pos <- hgvs_pos + nchar(deletion) - 1L
                sprintf("%s:g.%d_%ddelins%s", seqid, hgvs_pos, end_pos,
                        insertion)
            }
        }
    }, character(1), USE.NAMES = FALSE)
}


#' Map CDS position to genomic coordinates (internal helper)
#'
#' Requires a TxDb and BSgenome. Used by hgvs_to_vcf for c. variants.
#'
#' @keywords internal
.cds_to_genomic <- function(parsed, txdb, bsgenome) {
    # This requires transcript-to-genome mapping which is implemented
    # in hgvs_transcribe.R
    NULL
}
