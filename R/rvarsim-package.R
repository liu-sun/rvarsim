#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom BiocGenerics start end strand
#' @importFrom Biostrings getSeq
#' @importFrom ensembldb EnsDb transcripts TxidFilter filter listColumns
#' @importFrom GenomeInfoDb seqnames seqlevelsStyle seqlevels
#' @importFrom GenomicFeatures cdsBy fiveUTRsByTranscript threeUTRsByTranscript
#' @importFrom GenomicFeatures exonsBy transcripts isActiveSeq
#' @importFrom GenomicRanges GRanges GRangesList reduce start end strand
#' @importFrom GenomicRanges width resize flank
#' @importFrom IRanges IRanges subsetByOverlaps overlapsAny width
#' @importFrom S4Vectors mcols metadata
## usethis namespace: end
NULL

utils::globalVariables("tx_id")

#' rvarsim: R/Bioconductor Variant Simulator with HGVS Notation
#'
#' Simulates all possible single nucleotide variants across MANE Select
#' transcripts, generates HGVS-compliant variant descriptions, and provides
#' a comprehensive toolkit for HGVS parsing, validation, normalization,
#' format conversion, transcription mapping, translation, and liftover.
#'
#' @section Variant simulation:
#' \code{\link{simulate_variants}} orchestrates the full pipeline:
#' fetch MANE Select transcripts, extract CDS/UTR/splice structure,
#' generate all possible SNVs, and format in HGVS notation.
#'
#' @section HGVS manipulation:
#' \itemize{
#'   \item \code{\link{parse_hgvs}} — parse HGVS strings into
#'     structured objects
#'   \item \code{\link{validate_hgvs}} — syntactic and semantic
#'     validation
#'   \item \code{\link{normalize_hgvs}} — 3' shift and canonicalize
#'   \item \code{\link{hgvs_to_vcf}}, \code{\link{vcf_to_hgvs}} —
#'     VCF conversion
#'   \item \code{\link{transcribe_hgvs}}, \code{\link{c_to_g}},
#'     \code{\link{g_to_c}} — genome↔transcript mapping
#'   \item \code{\link{translate_hgvs}} — protein consequence
#'     prediction
#'   \item \code{\link{backtranslate_hgvs}} — protein→coding
#'     enumeration
#'   \item \code{\link{extract_hgvs}} — sequence→HGVS extraction
#'   \item \code{\link{liftover_hgvs}} — assembly liftover
#' }
#'
#' @name rvarsim-package
NULL
