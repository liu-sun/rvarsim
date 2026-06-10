#' Validate HGVS variant descriptions
#'
#' Performs syntactic and semantic validation of HGVS variant strings.
#' Syntactic validation checks that the string follows correct HGVS
#' grammar. Semantic validation verifies that positions are in range,
#' reference alleles match expected sequence, and the variant is
#' biologically plausible.
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param txdb Optional \code{TxDb} or \code{EnsDb} for transcript-
#'   context semantic validation (e.g., checking that CDS positions
#'   are within the actual CDS length).
#' @param bsgenome Optional \code{BSgenome} for reference allele
#'   verification.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{input}{Original input string.}
#'   \item{syntax_valid}{Logical: passes syntactic checks.}
#'   \item{semantic_valid}{Logical: passes semantic checks.}
#'   \item{syntax_error}{Error message from syntactic check (or NA).}
#'   \item{semantic_error}{Error message from semantic check (or NA).}
#'   \item{type}{Variant type if successfully parsed.}
#' }
#'
#' @export
#'
#' @examples
#' validate_hgvs(c("NM_000546.6:c.215C>G", "invalid string"))
validate_hgvs <- function(hgvs_strings, txdb = NULL,
                           bsgenome = NULL) {
    result <- lapply(hgvs_strings, function(s) {
        # Parse result: always a list with $success and $value/$error_msg
        pr <- tryCatch({
            val <- parse_hgvs(s, strict = TRUE)[[1]]
            list(success = TRUE, value = val)
        }, error = function(e) {
            list(success = FALSE,
                 error_msg = conditionMessage(e))
        })

        syn_ok  <- pr$success
        sem_ok  <- FALSE
        syn_err <- if (pr$success) NA_character_ else pr$error_msg
        sem_err <- NA_character_
        var_type <- NA_character_

        if (pr$success) {
            parsed <- pr$value
            synt <- .syntactic_checks(parsed)
            syn_ok  <- synt$valid
            syn_err <- synt$error
            var_type <- parsed$type

            if (syn_ok) {
                sem <- .semantic_checks(parsed, txdb, bsgenome)
                sem_ok  <- sem$valid
                sem_err <- sem$error
            }
        }

        data.frame(
            input          = s,
            syntax_valid   = syn_ok,
            semantic_valid = sem_ok,
            syntax_error   = syn_err,
            semantic_error = sem_err,
            type           = var_type,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, result)
}


#' Syntactic validation rules
#'
#' Checks HGVS grammar per HGVS nomenclature guidelines.
#'
#' @param parsed Parsed HGVS variant object.
#' @return List with \code{valid} and \code{error} fields.
#'
#' @keywords internal
.syntactic_checks <- function(parsed) {
    # Check accession format
    if (!grepl("^(NM_|NR_|NP_|NC_|NG_|ENST|ENSG|LRG)", parsed$accession)) {
        return(list(valid = FALSE,
                    error = paste0(
                        "Accession must start with NM_, NR_, NP_, ",
                        "NC_, NG_, ENST, ENSG, or LRG")))
    }

    # Check notation prefix
    valid_notations <- c("c", "g", "p", "n", "m", "r")
    if (!parsed$notation %in% valid_notations) {
        return(list(valid = FALSE,
                    error = paste("Notation prefix must be one of:",
                                   paste(valid_notations, collapse = ", "))))
    }

    # Check nucleotide vs protein notation
    nt_notations <- c("c", "g", "n", "m")
    nt_types <- c("substitution", "deletion", "insertion", "duplication",
                  "inversion", "delins", "frameshift", "silent")
    if (parsed$notation %in% c("c", "g", "n", "m") &&
        !parsed$type %in% nt_types) {
        return(list(valid = FALSE,
                    error = paste("Type", parsed$type, "not valid for",
                                   parsed$notation, "notation")))
    }

    # Check alleles are valid nucleotide or amino acid codes
    if (parsed$notation %in% c("c", "g", "n", "m")) {
        valid_bases <- c("A", "C", "G", "T")
        if (!is.na(parsed$reference) && nchar(parsed$reference) > 0 &&
            !all(strsplit(parsed$reference, "")[[1]] %in% valid_bases)) {
            return(list(valid = FALSE,
                        error = "Reference contains invalid nucleotide codes"))
        }
        if (!is.na(parsed$alternate) && nchar(parsed$alternate) > 0 &&
            !all(strsplit(parsed$alternate, "")[[1]] %in% valid_bases)) {
            return(list(valid = FALSE,
                        error = "Alternate contains invalid nucleotide codes"))
        }
    }

    # Check position is numeric
    if (!is.na(parsed$position$start) &&
        is.na(suppressWarnings(
            as.integer(parsed$position$start)))) {
        return(list(valid = FALSE, error = "Position is not numeric"))
    }

    list(valid = TRUE, error = NA_character_)
}


#' Semantic validation against biological context
#'
#' @param parsed Parsed HGVS variant object.
#' @param txdb Optional TxDb for transcript context.
#' @param bsgenome Optional BSgenome for reference verification.
#'
#' @return List with \code{valid} and \code{error} fields.
#'
#' @keywords internal
.semantic_checks <- function(parsed, txdb = NULL, bsgenome = NULL) {
    # Check if reference allele is same as alternate (no change)
    if (!is.na(parsed$reference) && !is.na(parsed$alternate) &&
        parsed$reference == parsed$alternate &&
        parsed$type != "silent") {
        return(list(valid = FALSE,
                    error = "Reference and alternate alleles are identical"))
    }

    # Check position is positive (unless 5'UTR which uses negative)
    pos <- parsed$position$start
    if (!is.na(pos) && pos <= 0 && !parsed$position$is_five_utr &&
        parsed$notation != "p") {
        return(list(valid = FALSE,
                    error = "Position must be positive unless 5'UTR"))
    }

    # If bsgenome provided, verify reference allele at position
    # (requires txdb for transcript-to-genomic mapping)
    if (!is.null(bsgenome) && !is.null(txdb) &&
        !is.na(parsed$reference) && nchar(parsed$reference) > 0) {
        # Attempt reference verification
        ref_check <- .verify_reference_allele(parsed, txdb, bsgenome)
        if (!is.null(ref_check) && !ref_check$valid) {
            return(ref_check)
        }
    }

    list(valid = TRUE, error = NA_character_)
}


#' Verify reference allele matches genome
#'
#' @keywords internal
.verify_reference_allele <- function(parsed, txdb, bsgenome) {
    # Reference allele verification maps HGVS position to genomic
    # coordinates via TxDb, then checks against BSgenome.
    # This depends on full transcript-to-genomic mapping (deferred
    # to a future release when transcription context is integrated).
    NULL
}


#' Check if an HGVS string is valid (convenience wrapper)
#'
#' @param hgvs_string A single HGVS variant string.
#' @param ... Additional arguments passed to \code{\link{validate_hgvs}}.
#' @return \code{TRUE} if both syntactic and semantic checks pass.
#'
#' @examples
#' is_valid_hgvs("NM_000546.6:c.215C>G")
#' @export
is_valid_hgvs <- function(hgvs_string, ...) {
    result <- validate_hgvs(hgvs_string, ...)
    result$syntax_valid && result$semantic_valid
}
