#' Map variants between aligned genomic sequences
#'
#' Converts HGVS variants between genome assemblies (e.g., hg19 to
#' hg38) using chain files via the \code{rtracklayer} package.
#' Supports both g. and c. notation (c. requires TxDb for the target
#' assembly).
#'
#' @param hgvs_strings Character vector of HGVS variant descriptions.
#' @param chain_file Path to a UCSC chain file (e.g.,
#'   \code{"hg19ToHg38.over.chain"}).
#' @param target_bsgenome A \code{BSgenome} object for the target
#'   assembly (used to resolve reference alleles in the target).
#'
#' @return Character vector of HGVS strings mapped to the target
#'   assembly, or \code{NA} for variants that cannot be mapped.
#'
#' @export
#'
#' @examples
#' # Liftover requires rtracklayer and a chain file:
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'     cat("rtracklayer available for liftover")
#' }
liftover_hgvs <- function(hgvs_strings, chain_file,
                            target_bsgenome = NULL) {

    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("rtracklayer is required for liftover. ",
             "Install with BiocManager::install('rtracklayer').")
    }

    # Validate that all inputs are g. notation before importing chain
    parsed_list <- lapply(hgvs_strings, function(s) {
        tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
    })

    has_non_g <- vapply(parsed_list, function(p) {
        !is.na(p[[1]]) && !identical(p$notation, "g")
    }, logical(1))

    if (any(has_non_g)) {
        warning("liftover currently supports g. notation only")
    }

    chain <- rtracklayer::import.chain(chain_file)

    vapply(hgvs_strings, function(s) {
        parsed <- tryCatch(parse_hgvs(s)[[1]], error = function(e) NA)
        if (is.na(parsed[[1]])) return(NA_character_)
        if (parsed$notation != "g") {
            warning("liftover_g_hgvs supports g. notation only")
            return(NA_character_)
        }

        .liftover_single(parsed, chain, target_bsgenome)
    }, character(1), USE.NAMES = FALSE)
}


#' Lift over a single g. variant
#'
#' @keywords internal
.liftover_single <- function(parsed, chain, target_bsgenome) {

    pos <- parsed$position$start
    end_pos <- if (!is.na(parsed$position$end)) parsed$position$end
        else pos
    chr <- parsed$accession

    # Create a GRanges for the variant position
    gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(start = pos, end = end_pos),
        strand   = "+"
    )

    # Perform liftover
    lifted <- tryCatch(
        rtracklayer::liftOver(gr, chain),
        error = function(e) {
            warning("Liftover failed: ", e$message)
            GenomicRanges::GRangesList()
        }
    )

    if (length(lifted) == 0 || length(lifted[[1]]) == 0) {
        warning("Position ", pos, " could not be lifted over")
        return(NA_character_)
    }

    target_gr <- lifted[[1]]

    if (length(target_gr) > 1) {
        warning("Position split into multiple mappings: using first")
    }

    target_chr <- as.character(GenomeInfoDb::seqnames(target_gr)[1])
    target_pos <- BiocGenerics::start(target_gr)[1]
    target_end <- BiocGenerics::end(target_gr)[1]

    # Check if lengths match (could be inversion or indel in chain)
    orig_width <- end_pos - pos + 1
    target_width <- target_end - target_pos + 1

    if (orig_width != target_width) {
        warning("Target interval width differs (", target_width,
                " vs ", orig_width, "). Variant may be ambiguous.")
    }

    # Verify reference allele in target assembly if BSgenome provided
    ref <- parsed$reference
    alt <- parsed$alternate

    if (!is.null(target_bsgenome) && !is.na(ref) &&
        nchar(ref) == 1 && nchar(alt) == 1) {

        target_ref <- tryCatch({
            s <- BSgenome::getSeq(target_bsgenome, target_chr,
                                   start = target_pos, end = target_pos)
            as.character(s)
        }, error = function(e) NA_character_)

        if (!is.na(target_ref) && target_ref != ref) {
            warning("Reference allele differs in target assembly: ",
                    "source=", ref, " target=", target_ref,
                    ". Using target allele.")
            ref <- target_ref
        }
    }

    # Build target HGVS string
    if (parsed$type == "substitution") {
        sprintf("%s:g.%d%s>%s", target_chr, target_pos, ref, alt)
    } else if (parsed$type == "deletion") {
        if (target_pos == target_end) {
            sprintf("%s:g.%ddel%s", target_chr, target_pos, ref)
        } else {
            sprintf("%s:g.%d_%ddel%s", target_chr, target_pos,
                    target_end, ref)
        }
    } else if (parsed$type == "insertion") {
        sprintf("%s:g.%d_%dins%s", target_chr,
                target_pos - 1L, target_pos, alt)
    } else {
        # Generic conversion
        sprintf("%s:g.%d%s>%s", target_chr, target_pos, ref, alt)
    }
}


#' Lift over a c. variant by converting to genomic, lifting,
#' then converting back
#'
#' @param hgvs_c A HGVS c. notation string.
#' @param source_txdb TxDb for the source assembly.
#' @param source_bsgenome BSgenome for the source assembly.
#' @param chain_file Path to chain file.
#' @param target_txdb TxDb for the target assembly.
#' @param target_bsgenome BSgenome for the target assembly.
#'
#' @return HGVS c. notation string in the target assembly, or NA.
#'
#' @examples
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'     cat("rtracklayer available for c. liftover")
#' }
#'
#' @export
liftover_c_hgvs <- function(hgvs_c, source_txdb, source_bsgenome,
                              chain_file, target_txdb,
                              target_bsgenome) {

    # Step 1: c. → g. in source assembly
    g_source <- c_to_g(hgvs_c, source_txdb, source_bsgenome)

    if (is.na(g_source[[1]])) return(NA_character_)

    # Step 2: Lift over genomic position
    g_target <- liftover_hgvs(g_source[[1]], chain_file,
                                target_bsgenome)

    if (is.na(g_target)) return(NA_character_)

    # Step 3: g. → c. in target assembly
    c_target <- g_to_c(g_target, target_txdb, target_bsgenome)

    if (is.na(c_target[[1]])) return(NA_character_)

    c_target[[1]]
}
