#' Infer transcript variants from protein consequences
#'
#' Given an observed protein (p.) variant, infers all possible
#' coding (c.) variants that could produce it by enumerating
#' possible codon changes in the genetic code.
#'
#' For example, p.Arg215Gly could be caused by any of the 6 AGA/Gly
#' codon pair changes. This function enumerates all possibilities.
#'
#' @param hgvs_p Character vector of HGVS p. notation strings
#'   (e.g., \code{"p.Arg215Gly"}).
#' @param txdb A \code{TxDb} or \code{EnsDb} object.
#' @param bsgenome A \code{BSgenome} object for retrieving the
#'   reference coding sequence.
#' @param genetic_code Genetic code table (default: standard code).
#'
#' @return A list of character vectors, where each element is a
#'   vector of possible HGVS c. notation strings for the
#'   corresponding protein variant.
#'
#' @export
#'
#' @examples
#' # The backtranslate function requires a TxDb and BSgenome with
#' # matching seqlevels (use same reference for both):
#' if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     cat("Ready for backtranslation")
#' }
backtranslate_hgvs <- function(hgvs_p, txdb, bsgenome,
                                genetic_code = NULL) {

    if (is.null(genetic_code)) {
        genetic_code <- .standard_genetic_code()
    }

    lapply(hgvs_p, function(s) {
        parsed <- .parse_protein_variant(s)
        if (is.null(parsed)) return(NA_character_)

        .enumerate_coding_variants(parsed, txdb, bsgenome,
                                    genetic_code)
    })
}


#' Parse a p. HGVS string (simplified parser for protein notation)
#'
#' @keywords internal
.parse_protein_variant <- function(p_str) {
    # Parse patterns like:
    # p.Arg215Gly     (missense)
    # p.Arg215Ter     (nonsense)
    # p.(=)           (silent)
    # p.Arg215LeufsTer5 (frameshift)
    # p.Lys215_Gly217del (deletion)
    # p.Lys215_Gly216insAla (insertion)

    if (!grepl("^p\\.", p_str)) {
        warning("Not a protein variant: ", p_str)
        return(NULL)
    }

    body <- sub("^p\\.", "", p_str)

    # Silent
    if (body == "(=)") {
        return(list(type = "silent"))
    }

    # Frameshift
    if (grepl("fs", body, fixed = TRUE)) {
        m <- regexec("^([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})fsTer([0-9]+)$",
                     body)
        parts <- regmatches(body, m)[[1]]
        if (length(parts) >= 5) {
            return(list(
                type    = "frameshift",
                ref_aa  = .aa_one_letter(parts[2]),
                pos     = as.integer(parts[3]),
                new_aa  = parts[4],
                ter_pos = as.integer(parts[5])
            ))
        }
    }

    # Nonsense: p.Arg215Ter or p.Arg215*
    if (grepl("Ter$", body) || grepl("\\*$", body)) {
        m <- regexec("^([A-Z][a-z]{2})([0-9]+)(Ter|\\*)$", body)
        parts <- regmatches(body, m)[[1]]
        if (length(parts) >= 4) {
            return(list(
                type   = "nonsense",
                ref_aa = .aa_one_letter(parts[2]),
                pos    = as.integer(parts[3])
            ))
        }
    }

    # Deletion: p.Lys215_Gly217del or p.Lys215del
    if (grepl("del$", body)) {
        m <- regexec(
            "^([A-Z][a-z]{2})([0-9]+)(?:_([A-Z][a-z]{2})([0-9]+))?del$",
            body, perl = TRUE
        )
        parts <- regmatches(body, m)[[1]]
        if (length(parts) >= 3) {
            end_aa <- if (nchar(parts[4]) > 0) parts[4] else parts[2]
            end_pos <- if (nchar(parts[5]) > 0) as.integer(parts[5])
                else as.integer(parts[3])
            return(list(
                type     = "deletion",
                ref_aa   = .aa_one_letter(parts[2]),
                pos      = as.integer(parts[3]),
                end_aa   = .aa_one_letter(end_aa),
                end_pos  = end_pos
            ))
        }
    }

    # Insertion: p.Lys215_Gly216insAla
    if (grepl("ins", body, fixed = TRUE)) {
        m <- regexec(
            paste0(
                "^([A-Z][a-z]{2})([0-9]+)",
                "_([A-Z][a-z]{2})([0-9]+)",
                "ins([A-Z][a-z]{3,})$"
            ),
            body
        )
        parts <- regmatches(body, m)[[1]]
        if (length(parts) >= 6) {
            return(list(
                type    = "insertion",
                ref_aa1 = .aa_one_letter(parts[2]),
                pos1    = as.integer(parts[3]),
                ref_aa2 = .aa_one_letter(parts[4]),
                pos2    = as.integer(parts[5]),
                ins_aa  = .parse_aa_three_letter_string(parts[6])
            ))
        }
    }

    # Missense: p.Arg215Gly
    m <- regexec("^([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})$", body)
    parts <- regmatches(body, m)[[1]]
    if (length(parts) >= 4) {
        return(list(
            type   = "missense",
            ref_aa = .aa_one_letter(parts[2]),
            pos    = as.integer(parts[3]),
            alt_aa = .aa_one_letter(parts[4])
        ))
    }

    warning("Could not parse protein variant: ", p_str)
    NULL
}


#' Parse a string of concatenated three-letter amino acid codes
#'
#' @keywords internal
.parse_aa_three_letter_string <- function(s) {
    # "AlaGlySer" → c("A", "G", "S")
    m <- gregexec("[A-Z][a-z]{2}", s)[[1]]
    if (length(m) == 0) return(character(0))
    codes <- regmatches(s, list(m))[[1]]
    vapply(codes, .aa_one_letter, character(1), USE.NAMES = FALSE)
}


#' Convert three-letter AA code to single-letter
#'
#' @keywords internal
.aa_one_letter <- function(aa3) {
    codes <- c(
        Ala = "A", Cys = "C", Asp = "D", Glu = "E",
        Phe = "F", Gly = "G", His = "H", Ile = "I",
        Lys = "K", Leu = "L", Met = "M", Asn = "N",
        Pro = "P", Gln = "Q", Arg = "R", Ser = "S",
        Thr = "T", Val = "V", Trp = "W", Tyr = "Y",
        Ter = "*", Xaa = "X"
    )
    if (aa3 %in% names(codes)) unname(codes[aa3]) else "X"
}


#' Enumerate all possible coding variants for a protein change
#'
#' @keywords internal
.enumerate_coding_variants <- function(protein_var, txdb,
                                        bsgenome, genetic_code) {

    # Get CDS for the transcript (need transcript ID — extract from
    # context or use first MANE Select transcript)
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx",
                                         use.names = TRUE)
    if (length(cds_by_tx) == 0) return(NA_character_)

    tx_id <- names(cds_by_tx)[1]
    cds_gr <- cds_by_tx[[tx_id]]
    tx_strand <- as.character(BiocGenerics::strand(cds_gr)[1])

    ref_cds <- .get_cds_sequence(cds_gr, bsgenome, tx_strand)
    if (length(ref_cds) != 1 || is.na(ref_cds) ||
        nchar(ref_cds) == 0) return(character(0))
    ref_protein <- .translate_sequence(ref_cds, genetic_code)

    variants <- switch(
        protein_var$type,
        missense = .enumerate_missense(protein_var, ref_cds,
                                        ref_protein, genetic_code,
                                        tx_id),
        nonsense = .enumerate_nonsense(protein_var, ref_cds,
                                        ref_protein, genetic_code,
                                        tx_id),
        frameshift = character(0),  # Too many possibilities
        deletion = character(0),
        insertion = character(0),
        silent = character(0),
        character(0)
    )

    variants
}


#' Enumerate missense codon changes
#'
#' @keywords internal
.enumerate_missense <- function(pv, ref_cds, ref_protein,
                                 genetic_code, tx_id) {
    aa_pos <- pv$pos
    codon_start <- (aa_pos - 1) * 3 + 1
    codon_end <- codon_start + 2

    if (is.na(ref_cds) || nchar(ref_cds) == 0) return(character(0))
    if (codon_end > nchar(ref_cds)) return(character(0))

    ref_codon <- substr(ref_cds, codon_start, codon_end)
    expected_aa <- genetic_code[ref_codon]
    if (is.na(expected_aa) || expected_aa != pv$ref_aa) {
        warning("Reference amino acid mismatch at position ", aa_pos)
    }

    # Find all codons that encode the target amino acid
    target_codons <- names(genetic_code)[genetic_code == pv$alt_aa]

    variants <- character(0)

    for (target_codon in target_codons) {
        # Find which positions differ between ref_codon and target_codon
        ref_bases <- strsplit(ref_codon, "")[[1]]
        alt_bases <- strsplit(target_codon, "")[[1]]
        diff_positions <- which(ref_bases != alt_bases)

        if (length(diff_positions) == 1) {
            # Single nucleotide change
            i <- diff_positions[1]
            cds_pos <- codon_start + i - 1
            variants <- c(variants,
                          sprintf("%s:c.%d%s>%s", tx_id,
                                  cds_pos,
                                  ref_bases[i], alt_bases[i]))
        } else if (length(diff_positions) == 2) {
            # Two nucleotide changes — may be represented as single
            # or two separate variants
            i1 <- diff_positions[1]
            i2 <- diff_positions[2]
            cds_pos1 <- codon_start + i1 - 1
            cds_pos2 <- codon_start + i2 - 1
            variants <- c(variants,
                          sprintf("%s:c.%d%s>%s (%d%s>%s)", tx_id,
                                  cds_pos1, ref_bases[i1],
                                  alt_bases[i1],
                                  cds_pos2, ref_bases[i2],
                                  alt_bases[i2]))
        } else if (length(diff_positions) == 3) {
            # Three changes: codon substitution with delins
            variants <- c(variants,
                          sprintf("%s:c.%d_%ddelins%s", tx_id,
                                  codon_start, codon_end,
                                  target_codon))
        }
    }

    unique(variants)
}


#' Enumerate nonsense codon changes
#'
#' @keywords internal
.enumerate_nonsense <- function(pv, ref_cds, ref_protein,
                                  genetic_code, tx_id) {
    aa_pos <- pv$pos
    codon_start <- (aa_pos - 1) * 3 + 1
    codon_end <- codon_start + 2

    if (is.na(ref_cds) || nchar(ref_cds) == 0) return(character(0))
    if (codon_end > nchar(ref_cds)) return(character(0))

    ref_codon <- substr(ref_cds, codon_start, codon_end)

    # Find all codons that are stop codons (*)
    stop_codons <- names(genetic_code)[genetic_code == "*"]

    variants <- character(0)

    for (stop_codon in stop_codons) {
        ref_bases <- strsplit(ref_codon, "")[[1]]
        alt_bases <- strsplit(stop_codon, "")[[1]]
        diff_positions <- which(ref_bases != alt_bases)

        for (i in diff_positions) {
            cds_pos <- codon_start + i - 1
            variants <- c(variants,
                          sprintf("%s:c.%d%s>%s", tx_id,
                                  cds_pos,
                                  ref_bases[i], alt_bases[i]))
        }

        if (length(diff_positions) == 2) {
            cds_pos1 <- codon_start + diff_positions[1] - 1
            cds_pos2 <- codon_start + diff_positions[2] - 1
            variants <- c(variants,
                          sprintf("%s:c.%d%s>%s [%d%s>%s]", tx_id,
                                  cds_pos1,
                                  ref_bases[diff_positions[1]],
                                  alt_bases[diff_positions[1]],
                                  cds_pos2,
                                  ref_bases[diff_positions[2]],
                                  alt_bases[diff_positions[2]]))
        }
    }

    unique(variants)
}
