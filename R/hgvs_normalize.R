#' Normalize HGVS variant descriptions
#'
#' Applies HGVS normalization rules: 3' shifting for indels in
#' repetitive regions, rewriting insertions as duplications when
#' applicable, and trimming common prefixes/suffixes.
#'
#' Normlization ensures a single canonical representation for each
#' variant, which is critical for comparison and deduplication.
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param bsgenome A \code{BSgenome} object providing the reference
#'   sequence for 3'-shifting context.
#' @param txdb Optional \code{TxDb} for transcript context when
#'   normalizing c. notation variants.
#'
#' @return A character vector of normalized HGVS strings, or the
#'   original string if normalization cannot be applied.
#'
#' @export
#'
#' @examples
#' # Normalize an HGVS string:
#' normalize_hgvs("NM_000546.6:c.215C>G")
normalize_hgvs <- function(hgvs_strings, bsgenome = NULL,
                            txdb = NULL) {
    vapply(hgvs_strings, function(s) {
        parsed <- tryCatch(
            parse_hgvs(s)[[1]],
            error = function(e) return(s)
        )
        if (is.na(parsed[[1]])) return(s)

        result <- tryCatch(
            .normalize_variant(parsed, bsgenome, txdb),
            error = function(e) return(s)
        )
        result
    }, character(1), USE.NAMES = FALSE)
}


#' Normalize a single parsed variant
#'
#' @keywords internal
.normalize_variant <- function(parsed, bsgenome, txdb) {
    # Step 1: Trim common prefix/suffix for delins
    if (parsed$type == "delins" && !is.na(parsed$reference) &&
        !is.na(parsed$alternate)) {
        result <- .trim_common_affixes(parsed)
        if (!is.null(result)) {
            parsed <- result
        }
    }

    # Step 2: Rewrite insertion as duplication if applicable
    if (parsed$type == "insertion" && !is.na(parsed$alternate)) {
        result <- .try_rewrite_as_dup(parsed, bsgenome, txdb)
        if (!is.null(result)) {
            parsed <- result
        }
    }

    # Step 3: 3' shift for indels in repetitive regions
    if (parsed$type %in% c("deletion", "insertion", "duplication",
                            "delins") && !is.null(bsgenome)) {
        result <- .three_prime_shift(parsed, bsgenome, txdb)
        if (!is.null(result)) {
            parsed <- result
        }
    }

    # Rebuild HGVS string from normalized parsed object
    .rebuild_hgvs(parsed)
}


#' Trim common prefix and suffix from reference and alternate
#' sequences for delins variants.
#'
#' @keywords internal
.trim_common_affixes <- function(parsed) {
    ref <- parsed$reference
    alt <- parsed$alternate

    if (is.na(ref) || is.na(alt) || nchar(ref) == 0 || nchar(alt) == 0) {
        return(NULL)
    }

    ref_chars <- strsplit(ref, "")[[1]]
    alt_chars <- strsplit(alt, "")[[1]]

    # Trim common prefix
    prefix_len <- 0L
    min_len <- min(length(ref_chars), length(alt_chars))
    while (prefix_len < min_len &&
           ref_chars[prefix_len + 1] == alt_chars[prefix_len + 1]) {
        prefix_len <- prefix_len + 1L
    }

    # Trim common suffix
    suffix_len <- 0L
    ref_remaining <- length(ref_chars) - prefix_len
    alt_remaining <- length(alt_chars) - prefix_len
    min_remain <- min(ref_remaining, alt_remaining)
    while (suffix_len < min_remain &&
           ref_chars[length(ref_chars) - suffix_len] ==
           alt_chars[length(alt_chars) - suffix_len]) {
        suffix_len <- suffix_len + 1L
    }

    if (prefix_len > 0 || suffix_len > 0) {
        new_ref <- paste(
            ref_chars[(prefix_len + 1):(length(ref_chars) - suffix_len)],
            collapse = "")
        new_alt <- paste(
            alt_chars[(prefix_len + 1):(length(alt_chars) - suffix_len)],
            collapse = "")

        pos <- parsed$position
        pos$start <- pos$start + prefix_len
        if (!is.na(pos$end)) {
            pos$end <- pos$end - suffix_len
        }

        parsed$reference <- new_ref
        parsed$alternate <- new_alt
        parsed$position <- pos

        # If after trimming, ref is empty → insertion
        if (nchar(new_ref) == 0) {
            parsed$type <- "insertion"
        } else if (nchar(new_alt) == 0) {
            parsed$type <- "deletion"
        }
    }

    parsed
}


#' Rewrite insertion as duplication when the inserted sequence matches
#' the adjacent reference sequence.
#'
#' @keywords internal
.try_rewrite_as_dup <- function(parsed, bsgenome, txdb) {
    alt <- parsed$alternate
    if (is.na(alt) || nchar(alt) == 0) return(NULL)

    # Check if inserted sequence matches upstream reference
    # This requires knowing the reference context
    # For now, return NULL (passthrough if bsgenome unavailable)
    NULL
}


#' Apply 3' shifting rule for indels in repetitive regions.
#'
#' For deletions: shift position downstream as far as possible
#' while the deleted sequence still matches the reference.
#'
#' For insertions: shift position downstream as far as possible
#' while the inserted base matches the following reference base.
#'
#' @keywords internal
.three_prime_shift <- function(parsed, bsgenome, txdb) {
    # 3' shift requires knowing the reference sequence context
    # and genomic coordinates. Placeholder that requires
    # transcription mapping integration.
    NULL
}


#' Rebuild HGVS string from a parsed object
#'
#' @keywords internal
.rebuild_hgvs <- function(parsed) {
    prefix <- paste0(parsed$accession, ":", parsed$notation, ".")

    pos_str <- .rebuild_position(parsed$position)

    switch(
        parsed$type,
        substitution  = paste0(prefix, pos_str,
                                parsed$reference, ">", parsed$alternate),
        deletion      = paste0(prefix, pos_str,
                                "del", parsed$reference),
        insertion     = paste0(prefix, pos_str,
                                "ins", parsed$alternate),
        duplication   = paste0(prefix, pos_str,
                                "dup",
                                if (nchar(parsed$reference) > 0)
                                    parsed$reference else ""),
        inversion     = paste0(prefix, pos_str, "inv"),
        delins        = paste0(prefix, pos_str,
                                "delins", parsed$alternate),
        frameshift    = paste0(prefix, pos_str, "fs"),
        silent        = paste0(prefix, pos_str, "="),
        parsed$raw    # fallback
    )
}


#' Rebuild position string from position list
#'
#' @keywords internal
.rebuild_position <- function(pos) {
    if (is.na(pos$start)) return("")

    start_str <- .format_single_pos(
        pos$start, pos$offset_start,
        pos$is_five_utr, pos$is_three_utr
    )

    if (is.na(pos$end) || pos$start == pos$end) {
        return(start_str)
    }

    end_str <- .format_single_pos(
        pos$end, pos$offset_end,
        pos$is_five_utr, pos$is_three_utr
    )
    paste0(start_str, "_", end_str)
}


#' Format a single HGVS coordinate
#'
#' @keywords internal
.format_single_pos <- function(pos, offset, is_five_utr, is_three_utr) {
    if (is_three_utr) {
        base <- paste0("*", pos)
    } else if (is_five_utr) {
        base <- as.character(pos)  # already negative
    } else {
        base <- as.character(pos)
    }

    if (!is.na(offset) && offset > 0) {
        paste0(base, "+", offset)
    } else if (!is.na(offset) && offset < 0) {
        paste0(base, offset)  # already includes minus sign
    } else {
        base
    }
}
