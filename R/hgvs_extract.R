#' Generate HGVS descriptions from reference and observed sequences
#'
#' Aligns a reference and observed sequence, identifies differences
#' (substitutions, insertions, deletions), and generates HGVS
#' variant descriptions.
#'
#' @param ref_seq Character string: reference sequence.
#' @param obs_seq Character string: observed (variant) sequence.
#' @param accession Accession number for the HGVS output
#'   (e.g., \code{"NM_000546.6"}).
#' @param notation HGVS notation prefix. Default \code{"c"}.
#' @param start_pos Starting position offset. For c. notation,
#'   1 = first base of CDS. For g. notation, the genomic position
#'   of the first base in the reference sequence.
#' @param algorithm Alignment algorithm: \code{"needleman"}
#'   (global) or \code{"smith"} (local). Default \code{"needleman"}.
#'
#' @return A character vector of HGVS variant description strings, one
#'   per detected variant. Returns \code{"="} if sequences are
#'   identical.
#'
#' @export
#'
#' @examples
#' ref <- "ATGCGTACGTAG"
#' obs <- "ATGCATACCTAG"
#' extract_hgvs(ref, obs, "NM_000546.6", "c", 1)
extract_hgvs <- function(ref_seq, obs_seq, accession,
                           notation = "c", start_pos = 1L,
                           algorithm = c("needleman", "smith")) {

    algorithm <- match.arg(algorithm)

    if (ref_seq == obs_seq) {
        return(sprintf("%s:%s.=", accession, notation))
    }

    # Pairwise alignment
    alignment <- .pairwise_align(ref_seq, obs_seq, algorithm)

    # Extract differences from alignment
    variants <- .extract_variants_from_alignment(
        alignment, accession, notation, start_pos
    )

    if (length(variants) == 0) {
        return(sprintf("%s:%s.=", accession, notation))
    }

    variants
}


#' Simple pairwise sequence alignment (Needleman-Wunsch or
#' Smith-Waterman)
#'
#' @keywords internal
.pairwise_align <- function(ref, obs, algorithm) {
    if (algorithm == "smith") {
        .smith_waterman(ref, obs)
    } else {
        .needleman_wunsch(ref, obs)
    }
}


#' Needleman-Wunsch global alignment
#'
#' @keywords internal
.needleman_wunsch <- function(seq1, seq2, match = 2, mismatch = -1,
                               gap_open = -2, gap_extend = -1) {

    n <- nchar(seq1) + 1
    m <- nchar(seq2) + 1

    s1 <- strsplit(seq1, "")[[1]]
    s2 <- strsplit(seq2, "")[[1]]

    # Score matrix
    score <- matrix(0, nrow = n, ncol = m)

    # Initialize gaps
    for (i in 2:n) score[i, 1] <- score[i - 1, 1] + gap_extend
    for (j in 2:m) score[1, j] <- score[1, j - 1] + gap_extend

    # Fill matrix
    for (i in 2:n) {
        for (j in 2:m) {
            diag <- score[i - 1, j - 1] +
                ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
            up   <- score[i - 1, j] + gap_extend
            left <- score[i, j - 1] + gap_extend
            score[i, j] <- max(diag, up, left)
        }
    }

    # Traceback
    aln1 <- character(0)
    aln2 <- character(0)
    i <- n
    j <- m

    while (i > 1 || j > 1) {
        if (i > 1 && j > 1 &&
            score[i, j] == score[i - 1, j - 1] +
            ifelse(s1[i - 1] == s2[j - 1], match, mismatch)) {
            aln1 <- c(s1[i - 1], aln1)
            aln2 <- c(s2[j - 1], aln2)
            i <- i - 1
            j <- j - 1
        } else if (i > 1 && score[i, j] == score[i - 1, j] + gap_extend) {
            aln1 <- c(s1[i - 1], aln1)
            aln2 <- c("-", aln2)
            i <- i - 1
        } else {
            aln1 <- c("-", aln1)
            aln2 <- c(s2[j - 1], aln2)
            j <- j - 1
        }
    }

    list(
        ref_aligned = paste(aln1, collapse = ""),
        obs_aligned = paste(aln2, collapse = ""),
        score = score[n, m]
    )
}


#' Smith-Waterman local alignment
#'
#' @keywords internal
.smith_waterman <- function(seq1, seq2, match = 2, mismatch = -1,
                              gap_open = -2, gap_extend = -1) {

    n <- nchar(seq1) + 1
    m <- nchar(seq2) + 1

    s1 <- strsplit(seq1, "")[[1]]
    s2 <- strsplit(seq2, "")[[1]]

    score <- matrix(0, nrow = n, ncol = m)
    max_score <- 0
    max_i <- 1
    max_j <- 1

    for (i in 2:n) {
        for (j in 2:m) {
            diag <- score[i - 1, j - 1] +
                ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
            up   <- score[i - 1, j] + gap_extend
            left <- score[i, j - 1] + gap_extend
            score[i, j] <- max(0, diag, up, left)

            if (score[i, j] > max_score) {
                max_score <- score[i, j]
                max_i <- i
                max_j <- j
            }
        }
    }

    # Traceback from max
    aln1 <- character(0)
    aln2 <- character(0)
    i <- max_i
    j <- max_j

    while (i > 1 && j > 1 && score[i, j] > 0) {
        if (score[i, j] == score[i - 1, j - 1] +
            ifelse(s1[i - 1] == s2[j - 1], match, mismatch)) {
            aln1 <- c(s1[i - 1], aln1)
            aln2 <- c(s2[j - 1], aln2)
            i <- i - 1
            j <- j - 1
        } else if (score[i, j] == score[i - 1, j] + gap_extend) {
            aln1 <- c(s1[i - 1], aln1)
            aln2 <- c("-", aln2)
            i <- i - 1
        } else {
            aln1 <- c("-", aln1)
            aln2 <- c(s2[j - 1], aln2)
            j <- j - 1
        }
    }

    list(
        ref_aligned = paste(aln1, collapse = ""),
        obs_aligned = paste(aln2, collapse = ""),
        score = max_score,
        start_pos = i  # 0-based position where alignment starts
    )
}


#' Extract HGVS variants from a pairwise alignment
#'
#' @keywords internal
.extract_variants_from_alignment <- function(aln, accession,
                                               notation, start_pos) {

    ref_aln <- aln$ref_aligned
    obs_aln <- aln$obs_aligned

    variants <- character(0)

    # Walk through alignment and identify differences
    n <- nchar(ref_aln)
    i <- 1L
    ref_pos <- start_pos  # position in reference coordinates

    while (i <= n) {
        ref_char <- substr(ref_aln, i, i)
        obs_char <- substr(obs_aln, i, i)

        if (ref_char == obs_char) {
            # Match: advance both
            ref_pos <- ref_pos + 1L
            i <- i + 1L
        } else if (ref_char == "-") {
            # Insertion in observed sequence
            ins_start <- i
            while (i <= n && substr(ref_aln, i, i) == "-") {
                i <- i + 1L
            }
            ins_seq <- gsub("-", "",
                            substr(obs_aln, ins_start, i - 1))

            # In HGVS c. notation, insertions are between positions
            if (notation == "c") {
                variants <- c(variants,
                              sprintf("%s:%s.%d_%dins%s",
                                      accession,
                                      ref_pos - 1, ref_pos, ins_seq))
            } else {
                variants <- c(variants,
                              sprintf("%s:g.%d_%dins%s",
                                      accession,
                                      ref_pos - 1, ref_pos, ins_seq))
            }
        } else if (obs_char == "-") {
            # Deletion in observed sequence
            del_start <- i
            while (i <= n && substr(obs_aln, i, i) == "-") {
                i <- i + 1L
            }
            del_seq <- gsub("-", "",
                            substr(ref_aln, del_start, i - 1))
            del_end_pos <- ref_pos + nchar(del_seq) - 1L

            if (nchar(del_seq) == 1) {
                variants <- c(variants,
                              sprintf("%s:%s.%ddel%s", accession,
                                      notation, ref_pos, del_seq))
            } else {
                variants <- c(variants,
                              sprintf("%s:%s.%d_%ddel%s", accession,
                                      notation, ref_pos, del_end_pos,
                                      del_seq))
            }
            ref_pos <- ref_pos + nchar(del_seq)
        } else {
            # Substitution (mismatch)
            variants <- c(variants,
                          sprintf("%s:%s.%d%s>%s", accession,
                                  notation, ref_pos,
                                  ref_char, obs_char))
            ref_pos <- ref_pos + 1L
            i <- i + 1L
        }
    }

    variants
}
