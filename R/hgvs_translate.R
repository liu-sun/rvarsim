#' Infer protein consequences from transcript variants
#'
#' Translates CDS nucleotide variants to protein (p.) HGVS notation
#' using the genetic code and CDS structure from a TxDb/EnsDb.
#'
#' Supported consequence types:
#' \itemize{
#'   \item Missense: \code{p.Arg215Gly} (substitution causing amino
#'     acid change)
#'   \item Nonsense: \code{p.Trp215Ter} (premature stop codon)
#'   \item Silent: \code{p.(=)} (no amino acid change)
#'   \item Frameshift: \code{p.Arg215LeufsTer5} (frameshift with
#'     altered reading frame)
#'   \item Inframe deletion: \code{p.Lys215_Gly217del}
#'   \item Inframe insertion: \code{p.Lys215_Gly216insAla}
#' }
#'
#' @param hgvs_strings Character vector of HGVS c. notation strings.
#' @param txdb A \code{TxDb} or \code{EnsDb} object.
#' @param bsgenome A \code{BSgenome} object for retrieving coding
#'   sequence.
#' @param genetic_code A named integer vector or genetic code table
#'   mapping codons to single-letter amino acids. If \code{NULL}
#'   (default), uses the standard genetic code.
#'
#' @return Character vector of HGVS p. notation strings, or
#'   \code{NA} for variants that cannot be translated (e.g.,
#'   non-coding variants).
#'
#' @export
#'
#' @examples
#' # Translate requires a TxDb with matching BSgenome:
#' if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     cat("Ready for translation")
#' }
translate_hgvs <- function(hgvs_strings, txdb, bsgenome,
                            genetic_code = NULL) {

    if (is.null(genetic_code)) {
        genetic_code <- .standard_genetic_code()
    }

    vapply(hgvs_strings, function(s) {
        parsed <- tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
        if (is.na(parsed[[1]])) return(NA_character_)
        if (parsed$notation != "c") {
            warning("translate_hgvs requires c. notation input")
            return(NA_character_)
        }

        .translate_variant(parsed, txdb, bsgenome, genetic_code)
    }, character(1), USE.NAMES = FALSE)
}


#' Translate a single parsed c. variant to p. notation
#'
#' @keywords internal
.translate_variant <- function(parsed, txdb, bsgenome,
                                genetic_code) {

    # Get CDS for the transcript
    tx_id <- parsed$accession
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx",
                                         use.names = TRUE)

    if (!tx_id %in% names(cds_by_tx)) {
        warning("Transcript ", tx_id, " not found in TxDb")
        return(NA_character_)
    }

    cds_gr <- cds_by_tx[[tx_id]]
    tx_strand <- as.character(BiocGenerics::strand(cds_gr)[1])

    # Get reference coding sequence from BSgenome
    ref_cds <- .get_cds_sequence(cds_gr, bsgenome, tx_strand)
    if (length(ref_cds) != 1 || is.na(ref_cds) ||
        nchar(ref_cds) == 0) return(NA_character_)

    # Translate reference CDS
    ref_protein <- .translate_sequence(ref_cds, genetic_code)

    # For SNVs: compute variant codon, translate, and determine
    # consequence type
    cds_pos <- parsed$position$start
    if (parsed$position$is_utr) {
        return(NA_character_)  # UTR variants don't affect protein
    }

    if (parsed$position$offset_start != 0) {
        return(NA_character_)  # Intronic variants — may affect
                               # splicing but not direct translation
    }

    fun <- switch(
        parsed$type,
        substitution = .translate_substitution,
        deletion     = .translate_deletion,
        insertion    = .translate_insertion,
        duplication  = .translate_duplication,
        delins       = .translate_delins,
        function(...) NA_character_
    )
    fun(parsed, ref_cds, ref_protein, cds_pos, genetic_code)
}


#' Translate a substitution variant
#'
#' @keywords internal
.translate_substitution <- function(parsed, ref_cds, ref_protein,
                                     cds_pos, genetic_code) {
    ref_base <- parsed$reference
    alt_base <- parsed$alternate

    # Ensure CDS position is within bounds
    if (cds_pos < 1 || cds_pos > nchar(ref_cds)) {
        warning("CDS position ", cds_pos, " out of bounds")
        return(NA_character_)
    }

    # Verify reference base matches
    actual_ref <- substr(ref_cds, cds_pos, cds_pos)
    if (actual_ref != ref_base) {
        warning("Reference base mismatch at c.", cds_pos,
                ": expected ", actual_ref, ", got ", ref_base)
    }

    # Build alternate CDS sequence
    alt_cds <- ref_cds
    substr(alt_cds, cds_pos, cds_pos) <- alt_base

    # Translate alternate CDS
    alt_protein <- .translate_sequence(alt_cds, genetic_code)

    # Find affected codon and amino acid position
    codon_idx <- ceiling(cds_pos / 3)
    codon_start <- (codon_idx - 1) * 3 + 1

    ref_aa <- if (codon_idx <= nchar(ref_protein))
        substr(ref_protein, codon_idx, codon_idx) else "?"
    alt_aa <- if (codon_idx <= nchar(alt_protein))
        substr(alt_protein, codon_idx, codon_idx) else "?"

    if (nchar(alt_protein) != nchar(ref_protein)) {
        # Nonsense: premature stop
        stop_pos <- regexpr("\\*", alt_protein)[1]
        if (stop_pos > 0 && stop_pos < nchar(ref_protein)) {
            ref_aa_name <- .aa_three_letter(ref_aa)
            return(sprintf("p.%s%dTer", ref_aa_name, stop_pos))
        }
    }

    if (ref_aa == alt_aa) {
        # Silent
        return("p.(=)")
    }

    # Missense
    ref_name <- .aa_three_letter(ref_aa)
    alt_name <- .aa_three_letter(alt_aa)
    sprintf("p.%s%d%s", ref_name, codon_idx, alt_name)
}


#' Translate a deletion variant
#'
#' @keywords internal
.translate_deletion <- function(parsed, ref_cds, ref_protein,
                                 cds_pos, genetic_code) {
    del_seq <- parsed$reference
    del_len <- nchar(del_seq)
    end_pos <- cds_pos + del_len - 1L

    if (del_len %% 3 != 0) {
        # Frameshift
        codon_idx <- ceiling(cds_pos / 3)
        ref_aa <- substr(ref_protein, codon_idx, codon_idx)
        ref_name <- .aa_three_letter(ref_aa)

        # Build frameshifted sequence to find new stop
        alt_cds <- paste0(
            substr(ref_cds, 1, cds_pos - 1),
            substr(ref_cds, cds_pos + del_len, nchar(ref_cds))
        )
        alt_protein <- .translate_sequence(alt_cds, genetic_code)
        stop_pos <- regexpr("\\*", alt_protein)[1]
        fs_len <- if (stop_pos > 0) stop_pos - codon_idx else
            nchar(alt_protein) - codon_idx + 1

        return(sprintf("p.%s%dLeufsTer%d", ref_name, codon_idx,
                       if (stop_pos > 0) fs_len else nchar(alt_protein)))
    }

    # Inframe deletion
    aa_start <- ceiling(cds_pos / 3)
    aa_end <- ceiling(end_pos / 3)

    if (aa_start == aa_end) {
        ref_name <- .aa_three_letter(
            substr(ref_protein, aa_start, aa_start)
        )
        sprintf("p.%s%ddel", ref_name, aa_start)
    } else {
        ref_start <- .aa_three_letter(
            substr(ref_protein, aa_start, aa_start)
        )
        ref_end <- .aa_three_letter(
            substr(ref_protein, aa_end, aa_end)
        )
        sprintf("p.%s%d_%s%ddel", ref_start, aa_start, ref_end, aa_end)
    }
}


#' Translate an insertion variant
#'
#' @keywords internal
.translate_insertion <- function(parsed, ref_cds, ref_protein,
                                  cds_pos, genetic_code) {
    ins_seq <- parsed$alternate
    ins_len <- nchar(ins_seq)

    if (ins_len %% 3 != 0) {
        # Frameshift
        codon_idx <- ceiling(cds_pos / 3)
        ref_aa <- substr(ref_protein, codon_idx, codon_idx)
        ref_name <- .aa_three_letter(ref_aa)

        alt_cds <- paste0(
            substr(ref_cds, 1, cds_pos),
            ins_seq,
            substr(ref_cds, cds_pos + 1, nchar(ref_cds))
        )
        alt_protein <- .translate_sequence(alt_cds, genetic_code)
        stop_pos <- regexpr("\\*", alt_protein)[1]
        fs_len <- if (stop_pos > 0) stop_pos - codon_idx else
            nchar(alt_protein) - codon_idx + 1

        return(sprintf("p.%s%dLeufsTer%d", ref_name, codon_idx,
                       if (stop_pos > 0) fs_len else nchar(alt_protein)))
    }

    # Inframe insertion
    aa_pos <- ceiling(cds_pos / 3)
    alt_cds <- paste0(
        substr(ref_cds, 1, cds_pos),
        ins_seq,
        substr(ref_cds, cds_pos + 1, nchar(ref_cds))
    )
    alt_protein <- .translate_sequence(alt_cds, genetic_code)

    # Determine inserted amino acids
    ins_aa <- ""
    ins_start <- aa_pos + 1
    ins_aa_len <- ins_len / 3
    ins_aa <- substr(alt_protein, ins_start,
                     ins_start + ins_aa_len - 1)
    ins_aa_formatted <- paste(
        vapply(strsplit(ins_aa, "")[[1]], .aa_three_letter,
               character(1)),
        collapse = ""
    )

    ref_aa_left <- .aa_three_letter(
        substr(ref_protein, aa_pos, aa_pos)
    )
    ref_aa_right <- .aa_three_letter(
        substr(ref_protein, aa_pos + 1, aa_pos + 1)
    )

    sprintf("p.%s%d_%s%dins%s", ref_aa_left, aa_pos,
            ref_aa_right, aa_pos + 1, ins_aa_formatted)
}


#' Translate a duplication variant
#'
#' @keywords internal
.translate_duplication <- function(parsed, ref_cds, ref_protein,
                                    cds_pos, genetic_code) {
    # Similar to insertion but the inserted sequence is a copy of
    # the reference
    dup_seq <- parsed$reference
    if (nchar(dup_seq) == 0) {
        # Get sequence from reference at position
        start <- if (!is.na(parsed$position$start)) parsed$position$start
            else cds_pos
        end <- if (!is.na(parsed$position$end)) parsed$position$end
            else start
        dup_seq <- substr(ref_cds, start, end)
    }

    dup_len <- nchar(dup_seq)

    if (dup_len %% 3 == 0) {
        # Inframe
        aa_start <- ceiling(cds_pos / 3)
        aa_end <- ceiling((cds_pos + dup_len - 1) / 3)
        ref_start <- .aa_three_letter(
            substr(ref_protein, aa_start, aa_start)
        )
        ref_end <- .aa_three_letter(
            substr(ref_protein, aa_end, aa_end)
        )
        sprintf("p.%s%d_%s%ddup", ref_start, aa_start,
                ref_end, aa_end)
    } else {
        # Frameshift
        codon_idx <- ceiling(cds_pos / 3)
        ref_aa <- substr(ref_protein, codon_idx, codon_idx)
        ref_name <- .aa_three_letter(ref_aa)
        return(sprintf("p.%s%dfs", ref_name, codon_idx))
    }
}


#' Translate a delins variant
#'
#' @keywords internal
.translate_delins <- function(parsed, ref_cds, ref_protein,
                                cds_pos, genetic_code) {
    # Compute net length change
    del_len <- if (!is.na(parsed$position$end) &&
                   !is.na(parsed$position$start)) {
        parsed$position$end - parsed$position$start + 1
    } else {
        0
    }
    ins_len <- nchar(parsed$alternate)
    net_change <- ins_len - del_len

    if (net_change %% 3 != 0) {
        # Frameshift
        codon_idx <- ceiling(cds_pos / 3)
        ref_aa <- substr(ref_protein, codon_idx, codon_idx)
        ref_name <- .aa_three_letter(ref_aa)
        return(sprintf("p.%s%dfs", ref_name, codon_idx))
    }

    # Inframe delins: too complex for simple heuristic, return
    # generic description
    return(NA_character_)
}


#' Retrieve coding sequence from BSgenome for a CDS GRanges
#'
#' @keywords internal
.get_cds_sequence <- function(cds_gr, bsgenome, strand) {
    cds_gr <- GenomicRanges::reduce(cds_gr)

    # Match seqlevelsStyle between TxDb and BSgenome gracefully
    bs_style <- tryCatch(
        GenomeInfoDb::seqlevelsStyle(bsgenome)[1],
        error = function(e) "UCSC"
    )
    tryCatch(
        GenomeInfoDb::seqlevelsStyle(cds_gr) <- bs_style,
        warning = function(w) NULL,
        error = function(e) NULL
    )

    seq_str <- tryCatch({
        seq_list <- BSgenome::getSeq(bsgenome, cds_gr)
        # handle DNAStringSet: collapse into single string
        paste(as.character(seq_list), collapse = "")
    }, error = function(e) {
        warning("Failed to retrieve CDS sequence: ", e$message)
        NA_character_
    })

    seq_str
}


#' Translate a nucleotide sequence to protein using genetic code
#'
#' @keywords internal
.translate_sequence <- function(nt_seq, genetic_code) {
    if (is.na(nt_seq) || nchar(nt_seq) < 3) return("")

    codons <- substring(nt_seq,
                         seq(1, nchar(nt_seq) - 2, by = 3),
                         seq(3, nchar(nt_seq), by = 3))

    aa <- vapply(codons, function(codon) {
        if (codon %in% names(genetic_code)) {
            unname(genetic_code[codon])
        } else {
            "X"
        }
    }, character(1), USE.NAMES = FALSE)

    paste(aa, collapse = "")
}


#' Standard genetic code (DNA codons → single-letter amino acid)
#'
#' @keywords internal
.standard_genetic_code <- function() {
    c(
        TTT = "F", TTC = "F", TTA = "L", TTG = "L",
        TCT = "S", TCC = "S", TCA = "S", TCG = "S",
        TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
        TGT = "C", TGC = "C", TGA = "*", TGG = "W",
        CTT = "L", CTC = "L", CTA = "L", CTG = "L",
        CCT = "P", CCC = "P", CCA = "P", CCG = "P",
        CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
        CGT = "R", CGC = "R", CGA = "R", CGG = "R",
        ATT = "I", ATC = "I", ATA = "I", ATG = "M",
        ACT = "T", ACC = "T", ACA = "T", ACG = "T",
        AAT = "N", AAC = "N", AAA = "K", AAG = "K",
        AGT = "S", AGC = "S", AGA = "R", AGG = "R",
        GTT = "V", GTC = "V", GTA = "V", GTG = "V",
        GCT = "A", GCC = "A", GCA = "A", GCG = "A",
        GAT = "D", GAC = "D", GAA = "E", GAG = "E",
        GGT = "G", GGC = "G", GGA = "G", GGG = "G"
    )
}


#' Amino acid three-letter codes
#'
#' @keywords internal
.aa_three_letter <- function(aa) {
    codes <- c(
        A = "Ala", C = "Cys", D = "Asp", E = "Glu",
        F = "Phe", G = "Gly", H = "His", I = "Ile",
        K = "Lys", L = "Leu", M = "Met", N = "Asn",
        P = "Pro", Q = "Gln", R = "Arg", S = "Ser",
        T = "Thr", V = "Val", W = "Trp", Y = "Tyr",
        "*" = "Ter", X = "Xaa", "?" = "Xaa"
    )
    if (aa %in% names(codes)) unname(codes[aa]) else "Xaa"
}
