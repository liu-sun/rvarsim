#' Parse HGVS variant descriptions
#'
#' Parses HGVS (Human Genome Variation Society)-formatted variant
#' strings into structured R list objects.
#' Supports all major HGVS notation types: substitution (>), deletion
#' (del), insertion (ins), duplication (dup), inversion (inv), deletion-
#' insertion (delins), repeat, and frameshift variants across c., g., p.,
#' n., m., and r. notation prefixes.
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param strict Logical. If \code{TRUE}, raises errors on parse
#'   failures. If \code{FALSE}, returns \code{NA} for unparseable
#'   strings with a warning. Default: \code{FALSE}.
#'
#' @return A list of parsed HGVS objects, each with components:
#' \describe{
#'   \item{type}{Variant type: "substitution", "deletion", "insertion",
#'     "duplication", "inversion", "delins", "frameshift", "repeat",
#'     or "unknown".}
#'   \item{accession}{Accession number (e.g., "NM_000546.6").}
#'   \item{notation}{Notation prefix: "c", "g", "p", "n", "m", "r".}
#'   \item{reference}{Reference allele/sequence.}
#'   \item{alternate}{Alternate allele/sequence.}
#'   \item{position}{List with \code{start}, \code{end},
#'     \code{offset_start}, \code{offset_end}, \code{is_utr},
#'     \code{utr_offset}.}
#'   \item{raw}{Original input string.}
#' }
#'
#' @export
#'
#' @examples
#' parse_hgvs(c("NM_000546.6:c.215C>G", "NC_000001.11:g.123456delA"))
parse_hgvs <- function(hgvs_strings, strict = FALSE) {
    lapply(hgvs_strings, function(s) {
        if (is.na(s)) return(NA)
        result <- tryCatch(
            .parse_single_hgvs(s),
            error = function(e) {
                if (strict) stop(e) else {
                    warning("Parse failed for '", s, "': ", e$message)
                    NA
                }
            }
        )
        result
    })
}


#' @keywords internal
.parse_single_hgvs <- function(hgvs_str) {
    # Split on first colon to separate accession:variant
    parts <- strsplit(hgvs_str, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Missing colon separator")

    accession <- parts[1]
    variant_part <- parts[2]

    # Extract notation prefix (c., g., p., n., m., r.)
    notation <- substr(variant_part, 1, 1)

    # Remove prefix for further parsing
    variant_body <- substring(variant_part, 3)  # after "c."

    # Determine variant type and parse accordingly
    type <- .detect_hgvs_type(variant_body)
    pos_ref_alt <- .parse_variant_body(variant_body, type)

    result <- list(
        type      = type,
        accession = accession,
        notation  = notation,
        reference = pos_ref_alt$reference,
        alternate = pos_ref_alt$alternate,
        position  = pos_ref_alt$position,
        raw       = hgvs_str
    )
    class(result) <- c("hgvs_variant", "list")
    result
}


#' Detect HGVS variant type from the variant body string
#'
#' @keywords internal
.detect_hgvs_type <- function(body) {
    if (grepl("delins", body, fixed = TRUE)) return("delins")
    if (grepl("fs", body, fixed = TRUE))    return("frameshift")
    if (grepl("del", body, fixed = TRUE))   return("deletion")
    if (grepl("ins", body, fixed = TRUE))   return("insertion")
    if (grepl("dup", body, fixed = TRUE))   return("duplication")
    if (grepl("inv", body, fixed = TRUE))   return("inversion")
    if (grepl(">", body, fixed = TRUE))     return("substitution")
    if (grepl("=", body, fixed = TRUE))     return("silent")
    if (grepl("fs", body, fixed = TRUE))    return("frameshift")
    if (grepl("[", body, fixed = TRUE))     return("repeat")
    "unknown"
}


#' Parse position, reference, and alternate from variant body
#'
#' @keywords internal
.parse_variant_body <- function(body, type) {
    switch(
        type,
        substitution  = .parse_substitution(body),
        deletion      = .parse_deletion(body),
        insertion     = .parse_insertion(body),
        duplication   = .parse_duplication(body),
        inversion     = .parse_inversion(body),
        delins        = .parse_delins(body),
        frameshift    = .parse_frameshift(body),
        `repeat`      = .parse_repeat(body),
        silent        = .parse_silent(body),
        list(reference = NA_character_, alternate = NA_character_,
             position = .parse_position(""))
    )
}


#' Parse substitution: e.g., 215C>G, -14G>A, *32T>C, 453+1G>T
#' @keywords internal
.parse_substitution <- function(body) {
    # Split on '>'
    parts <- strsplit(body, ">", fixed = TRUE)[[1]]
    ref_allele <- substr(parts[1], nchar(parts[1]), nchar(parts[1]))
    alt_allele <- parts[2]
    pos_str <- substr(parts[1], 1, nchar(parts[1]) - 1)

    list(
        reference = ref_allele,
        alternate = alt_allele,
        position  = .parse_position(pos_str)
    )
}


#' Parse deletion: e.g., 215delC, 215_217delTCA, 215+1delG
#' @keywords internal
.parse_deletion <- function(body) {
    parts <- strsplit(body, "del", fixed = TRUE)[[1]]
    deleted_seq <- parts[2]
    pos_str <- parts[1]

    list(
        reference = deleted_seq,
        alternate = "",
        position  = .parse_position(pos_str)
    )
}


#' Parse insertion: e.g., 215_216insTAG, 215+1insG
#' @keywords internal
.parse_insertion <- function(body) {
    parts <- strsplit(body, "ins", fixed = TRUE)[[1]]
    inserted_seq <- parts[2]
    pos_str <- parts[1]

    list(
        reference = "",
        alternate = inserted_seq,
        position  = .parse_position(pos_str)
    )
}


#' Parse duplication: e.g., 215dupC, 215_217dup
#' @keywords internal
.parse_duplication <- function(body) {
    parts <- strsplit(body, "dup", fixed = TRUE)[[1]]
    dup_seq <- if (nchar(parts[2]) > 0) parts[2] else ""
    pos_str <- parts[1]

    list(
        reference = dup_seq,
        alternate = paste0(dup_seq, dup_seq),
        position  = .parse_position(pos_str)
    )
}


#' Parse inversion: e.g., 215_217inv
#' @keywords internal
.parse_inversion <- function(body) {
    pos_str <- sub("inv", "", body, fixed = TRUE)
    list(
        reference = NA_character_,
        alternate = NA_character_,
        position  = .parse_position(pos_str)
    )
}


#' Parse delins: e.g., 215_217delinsTAG
#' @keywords internal
.parse_delins <- function(body) {
    parts <- strsplit(body, "delins", fixed = TRUE)[[1]]
    pos_str <- parts[1]
    inserted_seq <- parts[2]

    list(
        reference = NA_character_,  # requires reference context
        alternate = inserted_seq,
        position  = .parse_position(pos_str)
    )
}


#' Parse frameshift: e.g., 215fs, 215delGfs*5
#' @keywords internal
.parse_frameshift <- function(body) {
    # Extract the numeric position: everything after c. and before the
    # first non-numeric character. e.g., "215delGfs*5" → "215"
    m <- regexec("^([0-9]+)", body)
    pos_match <- regmatches(body, m)[[1]]
    pos_str <- if (length(pos_match) >= 2) pos_match[2] else ""
    list(
        reference = NA_character_,
        alternate = NA_character_,
        position  = .parse_position(pos_str)
    )
}


#' Parse repeat variant: e.g., \code{215[3]}, 215_217(2)
#' @keywords internal
.parse_repeat <- function(body) {
    parts <- strsplit(body, "[", fixed = TRUE)[[1]]
    pos_str <- parts[1]
    count <- as.integer(gsub("]", "", parts[2], fixed = TRUE))
    list(
        reference = NA_character_,
        alternate = NA_character_,
        position  = .parse_position(pos_str)
    )
}


#' Parse silent (no change): e.g., 215=
#' @keywords internal
.parse_silent <- function(body) {
    pos_str <- sub("=", "", body, fixed = TRUE)
    list(
        reference = "=",
        alternate = "=",
        position  = .parse_position(pos_str)
    )
}


#' Parse HGVS position notation
#'
#' Parses position expressions like:
#' \itemize{
#'   \item 215 (simple CDS position)
#'   \item -14 (5' UTR position)
#'   \item *32 (3' UTR position)
#'   \item 453+1 (intronic donor offset)
#'   \item 612-2 (intronic acceptor offset)
#'   \item 215_217 (range)
#'   \item -14_-12 (UTR range)
#'   \item 453+1_453+2 (intronic range)
#' }
#'
#' @param pos_str Position string from HGVS notation.
#'
#' @return A list with \code{start}, \code{end}, \code{offset_start},
#'   \code{offset_end}, \code{is_utr}, \code{utr_offset}.
#'
#' @keywords internal
.parse_position <- function(pos_str) {
    if (nchar(pos_str) == 0) {
        return(list(start = NA_integer_, end = NA_integer_,
                    offset_start = 0L, offset_end = 0L,
                    is_utr = FALSE, is_five_utr = FALSE,
                    is_three_utr = FALSE))
    }

    # Split on range separator "_"
    if (grepl("_", pos_str, fixed = TRUE)) {
        range_parts <- strsplit(pos_str, "_", fixed = TRUE)[[1]]
        start_parsed <- .parse_single_pos(range_parts[1])
        end_parsed   <- .parse_single_pos(range_parts[2])

        # If one is UTR, both should be treated as UTR positions
        if (start_parsed$is_utr) end_parsed$is_utr <- TRUE
        if (end_parsed$is_utr)   start_parsed$is_utr <- TRUE

        return(list(
            start          = start_parsed$pos,
            end            = end_parsed$pos,
            offset_start   = start_parsed$offset,
            offset_end     = end_parsed$offset,
            is_utr         = start_parsed$is_utr || end_parsed$is_utr,
            is_five_utr    = start_parsed$is_five_utr,
            is_three_utr   = start_parsed$is_three_utr
        ))
    }

    # Single position
    parsed <- .parse_single_pos(pos_str)
    list(
        start        = parsed$pos,
        end          = parsed$pos,
        offset_start = parsed$offset,
        offset_end   = parsed$offset,
        is_utr       = parsed$is_utr,
        is_five_utr  = parsed$is_five_utr,
        is_three_utr = parsed$is_three_utr
    )
}


#' Parse a single HGVS coordinate (not a range)
#'
#' @keywords internal
.parse_single_pos <- function(coord) {
    is_five_utr <- FALSE
    is_three_utr <- FALSE
    offset <- 0L

    # 3' UTR: *32
    if (grepl("^\\*", coord)) {
        is_three_utr <- TRUE
        pos <- as.integer(sub("^\\*", "", coord))
        return(list(pos = pos, offset = 0L, is_utr = TRUE,
                    is_five_utr = FALSE, is_three_utr = TRUE))
    }

    # Check for intronic offsets: +N or -N
    offset <- 0L
    pos_with_offset <- coord

    if (grepl("+", coord, fixed = TRUE)) {
        parts <- strsplit(coord, "+", fixed = TRUE)[[1]]
        pos_with_offset <- parts[1]
        offset <- as.integer(parts[2])
    } else if (grepl("-", coord, fixed = TRUE) &&
               !grepl("^-", coord)) {
        # Handle case like 612-2 but not -14
        # Split on last '-' not at start
        m <- regexec("^(-?[0-9]+|\\*[0-9]+)-([0-9]+)$", coord)
        m_parts <- regmatches(coord, m)[[1]]
        if (length(m_parts) >= 3) {
            pos_with_offset <- m_parts[2]
            offset <- -as.integer(m_parts[3])
        }
    }

    if (grepl("^-", pos_with_offset)) {
        is_five_utr <- TRUE
        pos <- as.integer(pos_with_offset)  # negative
    } else {
        pos <- as.integer(pos_with_offset)
    }

    list(
        pos          = pos,
        offset       = offset,
        is_utr       = is_five_utr || is_three_utr,
        is_five_utr  = is_five_utr,
        is_three_utr = is_three_utr
    )
}


#' @export
print.hgvs_variant <- function(x, ...) {
    cat("HGVS Variant: ", x$raw, "\n", sep = "")
    cat("  Type:      ", x$type, "\n", sep = "")
    cat("  Accession: ", x$accession, "\n", sep = "")
    cat("  Notation:  ", x$notation, ".\n", sep = "")
    cat("  Reference: ", x$reference, "\n", sep = "")
    cat("  Alternate: ", x$alternate, "\n", sep = "")
    if (!is.na(x$position$start)) {
        pos_str <- as.character(x$position$start)
        if (!is.na(x$position$end) &&
            x$position$start != x$position$end) {
            pos_str <- paste0(pos_str, "_", x$position$end)
        }
        cat("  Position:  ", pos_str, "\n", sep = "")
    }
    invisible(x)
}


#' @export
as.character.hgvs_variant <- function(x, ...) {
    x$raw
}
