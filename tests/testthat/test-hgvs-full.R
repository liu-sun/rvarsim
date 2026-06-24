test_that("parse_hgvs handles substitution", {
    result <- parse_hgvs("NM_000546.6:c.215C>G")[[1]]
    expect_equal(result$type, "substitution")
    expect_equal(result$accession, "NM_000546.6")
    expect_equal(result$notation, "c")
    expect_equal(result$reference, "C")
    expect_equal(result$alternate, "G")
    expect_equal(result$position$start, 215)
})

test_that("parse_hgvs handles deletion", {
    result <- parse_hgvs("NM_000546.6:c.215delC")[[1]]
    expect_equal(result$type, "deletion")
    expect_equal(result$reference, "C")
    expect_equal(result$alternate, "")
})

test_that("parse_hgvs handles insertion", {
    result <- parse_hgvs("NM_000546.6:c.215_216insTAG")[[1]]
    expect_equal(result$type, "insertion")
    expect_equal(result$alternate, "TAG")
})

test_that("parse_hgvs handles duplication", {
    result <- parse_hgvs("NM_000546.6:c.215dupC")[[1]]
    expect_equal(result$type, "duplication")
})

test_that("parse_hgvs handles delins", {
    result <- parse_hgvs("NM_000546.6:c.215_217delinsTAG")[[1]]
    expect_equal(result$type, "delins")
    expect_equal(result$alternate, "TAG")
})

test_that("parse_hgvs handles 5'UTR", {
    result <- parse_hgvs("NM_000546.6:c.-14G>A")[[1]]
    expect_true(result$position$is_five_utr)
    expect_equal(result$position$start, -14)
})

test_that("parse_hgvs handles 3'UTR", {
    result <- parse_hgvs("NM_000546.6:c.*32T>C")[[1]]
    expect_true(result$position$is_three_utr)
    expect_equal(result$position$start, 32)
})

test_that("parse_hgvs handles splice donor", {
    result <- parse_hgvs("NM_000546.6:c.453+1G>T")[[1]]
    expect_equal(result$position$offset_start, 1)
    expect_equal(result$position$start, 453)
})

test_that("parse_hgvs handles splice acceptor", {
    result <- parse_hgvs("NM_000546.6:c.612-2A>G")[[1]]
    expect_equal(result$position$offset_start, -2)
})

test_that("parse_hgvs handles genomic notation", {
    result <- parse_hgvs("NC_000001.11:g.123456A>G")[[1]]
    expect_equal(result$notation, "g")
    expect_equal(result$accession, "NC_000001.11")
})

test_that("parse_hgvs handles silent variant", {
    result <- parse_hgvs("NM_000546.6:c.215=")[[1]]
    expect_equal(result$type, "silent")
})

test_that("parse_hgvs handles frameshift", {
    result <- parse_hgvs("NM_000546.6:c.215delGfs*5")[[1]]
    expect_equal(result$type, "frameshift")
})

test_that("parse_hgvs handles multiple inputs", {
    results <- parse_hgvs(c("NM_000546.6:c.215C>G",
                             "NC_000001.11:g.123456A>G"))
    expect_length(results, 2)
    expect_equal(results[[1]]$type, "substitution")
    expect_equal(results[[2]]$type, "substitution")
})

test_that("parse_hgvs returns NA for invalid input in non-strict mode", {
    expect_warning(result <- parse_hgvs("garbage")[[1]])
    expect_true(is.na(result))
})

test_that("parse_hgvs raises error in strict mode", {
    expect_error(parse_hgvs("garbage", strict = TRUE))
})

# --- Validation tests ---

test_that("validate_hgvs accepts valid variants", {
    result <- validate_hgvs("NM_000546.6:c.215C>G")
    expect_true(result$syntax_valid)
})

test_that("validate_hgvs rejects invalid strings", {
    result <- validate_hgvs("not valid hgvs")
    expect_false(result$syntax_valid)
})

test_that("is_valid_hgvs is a convenience wrapper", {
    expect_true(is_valid_hgvs("NM_000546.6:c.215C>G"))
    expect_false(is_valid_hgvs("garbage"))
})

# --- Conversion tests ---

test_that("vcf_to_hgvs converts SNV correctly", {
    vcf <- data.frame(CHROM = "1", POS = 123456L, REF = "A", ALT = "G")
    result <- vcf_to_hgvs(vcf, "GRCh38")
    expect_true(grepl(":g\\.123456A>G", result))
})

test_that("hgvs_to_spdi converts correctly", {
    result <- hgvs_to_spdi("NC_000001.11:g.123456A>G")
    expect_equal(result, "NC_000001.11:123455:A:G")
})

test_that("spdi_to_hgvs converts correctly", {
    result <- spdi_to_hgvs("NC_000001.11:123455:A:G")
    expect_equal(result, "NC_000001.11:g.123456A>G")
})

# --- Normalization tests ---

test_that("normalize_hgvs handles valid input", {
    result <- normalize_hgvs("NM_000546.6:c.215C>G")
    expect_equal(result, "NM_000546.6:c.215C>G")
})

test_that(".trim_common_affixes removes shared prefix/suffix", {
    # Simulate a parsed delins with common prefix
    # This tests the internal helper
    parsed <- list(
        type = "delins", reference = "ATG", alternate = "ACG",
        accession = "TEST", notation = "c",
        position = list(start = 1, end = 3, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    result <- .trim_common_affixes(parsed)
    expect_equal(result$reference, "T")
    expect_equal(result$alternate, "C")
})

# --- Extraction tests ---

test_that("extract_hgvs detects substitution", {
    result <- extract_hgvs("ATGCGT", "ATGCAT", "NM_000546.6", "c", 1)
    expect_true(any(grepl("4G>A|4C>", result)) ||
                any(grepl("5G>A|5C>", result)))
})

test_that("extract_hgvs detects deletion", {
    result <- extract_hgvs("ATGCGT", "ATGGT", "NM_000546.6", "c", 1)
    expect_true(any(grepl("del", result)))
})

test_that("extract_hgvs returns = for identical sequences", {
    result <- extract_hgvs("ATGCGT", "ATGCGT", "NM_000546.6", "c", 1)
    expect_true(grepl("=$", result))
})

# --- Translation tests ---

test_that(".translate_sequence translates correctly", {
    gc <- .standard_genetic_code()
    # ATG = M (Met), GCA = A (Ala), TGA = * (stop)
    result <- .translate_sequence("ATGGCA", gc)
    expect_equal(result, "MA")
    result2 <- .translate_sequence("ATGTGA", gc)
    expect_equal(result2, "M*")
})

test_that(".aa_three_letter converts correctly", {
    expect_equal(.aa_three_letter("A"), "Ala")
    expect_equal(.aa_three_letter("R"), "Arg")
    expect_equal(.aa_three_letter("*"), "Ter")
})

# --- Backtranslation tests ---

test_that(".enumerate_missense finds codon changes", {
    gc <- .standard_genetic_code()
    ref_cds <- "CGTGGC"  # Arg-Gly
    ref_prot <- "RG"
    # Arg(CGT) → Gly(Gly): GGT, GGC, GGA, GGG
    pv <- list(type = "missense", pos = 1L, ref_aa = "R", alt_aa = "G")
    result <- .enumerate_missense(pv, ref_cds, ref_prot, gc, "TEST")

    # CGT→GGT: C>G at pos 1
    # CGT→GGC: C>G at pos 1, T>C at pos 2
    # CGT→GGA: C>G at pos 1, T>A at pos 2
    # CGT→GGG: C>G at pos 1
    expect_true(any(grepl("c\\.1C>G", result)))
})

test_that(".aa_one_letter converts correctly", {
    expect_equal(.aa_one_letter("Arg"), "R")
    expect_equal(.aa_one_letter("Gly"), "G")
    expect_equal(.aa_one_letter("Ter"), "*")
})

# --- Liftover tests ---

test_that("liftover_hgvs requires rtracklayer", {
    skip("rtracklayer not available in test environment")
})

test_that("liftover_hgvs warns for non-g notation", {
    skip("rtracklayer not available in test environment")
})

# --- Normalization tests ---

test_that("normalize_hgvs returns input on unparseable strings", {
    result <- normalize_hgvs("not an hgvs string")
    expect_equal(result, "not an hgvs string")
})

test_that(".trim_common_affixes handles no common affix", {
    parsed <- list(
        type = "delins", reference = "ATG", alternate = "TAA",
        position = list(start = 1, end = 3, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    result <- .trim_common_affixes(parsed)
    expect_equal(result$reference, "ATG")
    expect_equal(result$alternate, "TAA")
})

test_that(".rebuild_hgvs reconstructs substitution", {
    parsed <- list(
        type = "substitution", accession = "NM_TEST", notation = "c",
        reference = "C", alternate = "G",
        position = list(start = 215, end = 215, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.215C>G")
})

test_that(".rebuild_hgvs reconstructs deletion", {
    parsed <- list(
        type = "deletion", accession = "NM_TEST", notation = "c",
        reference = "A", alternate = "",
        position = list(start = 215, end = 215, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.215delA")
})

test_that(".rebuild_hgvs reconstructs 5'UTR variant", {
    parsed <- list(
        type = "substitution", accession = "NM_TEST", notation = "c",
        reference = "G", alternate = "A",
        position = list(start = -14, end = -14, offset_start = 0,
                        offset_end = 0, is_utr = TRUE,
                        is_five_utr = TRUE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.-14G>A")
})

test_that(".rebuild_hgvs reconstructs 3'UTR variant", {
    parsed <- list(
        type = "substitution", accession = "NM_TEST", notation = "c",
        reference = "T", alternate = "C",
        position = list(start = 32, end = 32, offset_start = 0,
                        offset_end = 0, is_utr = TRUE,
                        is_five_utr = FALSE, is_three_utr = TRUE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.*32T>C")
})

test_that(".rebuild_hgvs reconstructs splice donor", {
    parsed <- list(
        type = "substitution", accession = "NM_TEST", notation = "c",
        reference = "G", alternate = "T",
        position = list(start = 453, end = 453, offset_start = 1,
                        offset_end = 1, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.453+1G>T")
})

test_that(".rebuild_hgvs reconstructs insertion", {
    parsed <- list(
        type = "insertion", accession = "NM_TEST", notation = "c",
        reference = "", alternate = "TAG",
        position = list(start = 215, end = 216, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.215_216insTAG")
})

test_that(".rebuild_hgvs reconstructs silent variant", {
    parsed <- list(
        type = "silent", accession = "NM_TEST", notation = "c",
        reference = "=", alternate = "=",
        position = list(start = 215, end = 215, offset_start = 0,
                        offset_end = 0, is_utr = FALSE,
                        is_five_utr = FALSE, is_three_utr = FALSE)
    )
    expect_equal(.rebuild_hgvs(parsed), "NM_TEST:c.215=")
})

test_that(".rebuild_hgvs handles NA position gracefully", {
    parsed <- list(
        type = "substitution", accession = "NM_TEST", notation = "c",
        reference = "C", alternate = "G",
        position = list(start = NA_integer_, end = NA_integer_,
                        offset_start = 0, offset_end = 0,
                        is_utr = FALSE, is_five_utr = FALSE,
                        is_three_utr = FALSE)
    )
    expect_true(grepl(">G", .rebuild_hgvs(parsed)))
})

# --- Needleman-Wunsch alignment tests ---

test_that(".needleman_wunsch aligns identical sequences", {
    aln <- .needleman_wunsch("ATGC", "ATGC")
    expect_equal(aln$ref_aligned, "ATGC")
    expect_equal(aln$obs_aligned, "ATGC")
})

test_that(".needleman_wunsch handles substitution", {
    aln <- .needleman_wunsch("ATGC", "ATAC")
    expect_true(grepl("T", aln$ref_aligned))
    expect_true(grepl("A", aln$obs_aligned))
    expect_equal(nchar(aln$ref_aligned), nchar(aln$obs_aligned))
})

test_that(".needleman_wunsch handles deletion", {
    aln <- .needleman_wunsch("ATGC", "ATC")
    expect_true(grepl("-", aln$ref_aligned) ||
                grepl("-", aln$obs_aligned))
})

test_that(".needleman_wunsch handles insertion", {
    aln <- .needleman_wunsch("ATC", "ATGC")
    expect_true(grepl("-", aln$ref_aligned) ||
                grepl("-", aln$obs_aligned))
})

# --- Smith-Waterman alignment tests ---

test_that(".smith_waterman finds local alignment", {
    aln <- .smith_waterman("GGGATCGGG", "ATC")
    expect_gt(aln$score, 0)
})

# --- Conversion edge case tests ---

test_that("vcf_to_hgvs handles deletion", {
    vcf <- data.frame(CHROM = "1", POS = 123456L, REF = "AG", ALT = "A")
    result <- vcf_to_hgvs(vcf, "GRCh38")
    expect_true(grepl("delAG", result))
})

test_that("vcf_to_hgvs handles insertion", {
    vcf <- data.frame(CHROM = "1", POS = 123456L, REF = "A", ALT = "ATC")
    result <- vcf_to_hgvs(vcf, "GRCh38")
    expect_true(grepl("delinsATC", result))
})

test_that("spdi_to_hgvs handles deletion SPDI", {
    result <- spdi_to_hgvs("NC_000001.11:123455:AG:A")
    expect_true(grepl("delAG", result) || grepl(">", result))
})

test_that("hgvs_to_vcf returns NA on unparseable input", {
    result <- hgvs_to_vcf("not valid")
    expect_true(is.na(result$POS[1]))
})

# --- Format HGVS integration tests ---

test_that("format_hgvs handles mixed region types", {
    df <- data.frame(
        transcript_id = rep("NM_TEST", 3),
        region = c("cds", "five_utr", "three_utr"),
        genomic_pos = c(215L, 30L, 500L),
        genomic_ref = c("C", "G", "T"),
        genomic_alt = c("G", "A", "C"),
        cds_pos = c(215L, -14L, 32L),
        exon_boundary = NA_integer_,
        offset = NA_integer_,
        stringsAsFactors = FALSE
    )
    result <- format_hgvs(df)
    expect_equal(result$hgvs_c[1], "NM_TEST:c.215C>G")
    expect_equal(result$hgvs_c[2], "NM_TEST:c.-14G>A")
    expect_equal(result$hgvs_c[3], "NM_TEST:c.*32T>C")
})

test_that("format_hgvs handles splice sites", {
    df <- data.frame(
        transcript_id = rep("NM_TEST", 2),
        region = c("splice_donor", "splice_acceptor"),
        genomic_pos = c(454L, 610L),
        genomic_ref = c("G", "A"),
        genomic_alt = c("T", "G"),
        cds_pos = NA_integer_,
        exon_boundary = c(453L, 612L),
        offset = c(1L, -2L),
        stringsAsFactors = FALSE
    )
    result <- format_hgvs(df)
    expect_equal(result$hgvs_c[1], "NM_TEST:c.453+1G>T")
    expect_equal(result$hgvs_c[2], "NM_TEST:c.612-2A>G")
})

# --- Parse edge cases ---

test_that("parse_hgvs handles range positions", {
    result <- parse_hgvs("NM_000546.6:c.215_217delTCA")[[1]]
    expect_equal(result$position$start, 215)
    expect_equal(result$position$end, 217)
    expect_equal(result$reference, "TCA")
})

test_that("parse_hgvs handles duplication with no sequence", {
    result <- parse_hgvs("NM_000546.6:c.215_217dup")[[1]]
    expect_equal(result$type, "duplication")
})

test_that("parse_hgvs handles inversion", {
    result <- parse_hgvs("NM_000546.6:c.215_217inv")[[1]]
    expect_equal(result$type, "inversion")
})

test_that("as.character.hgvs_variant returns raw string", {
    v <- parse_hgvs("NM_000546.6:c.215C>G")[[1]]
    expect_equal(as.character(v), "NM_000546.6:c.215C>G")
})

# --- Protein variant parser tests ---

test_that(".parse_protein_variant handles missense", {
    result <- .parse_protein_variant("p.Arg215Gly")
    expect_equal(result$type, "missense")
    expect_equal(result$ref_aa, "R")
    expect_equal(result$pos, 215)
    expect_equal(result$alt_aa, "G")
})

test_that(".parse_protein_variant handles nonsense", {
    result <- .parse_protein_variant("p.Arg215Ter")
    expect_equal(result$type, "nonsense")
})

test_that(".parse_protein_variant handles silent", {
    result <- .parse_protein_variant("p.(=)")
    expect_equal(result$type, "silent")
})

test_that(".parse_protein_variant handles deletion", {
    result <- .parse_protein_variant("p.Lys215del")
    expect_equal(result$type, "deletion")
    expect_equal(result$pos, 215)
})

test_that(".parse_protein_variant handles range deletion", {
    result <- .parse_protein_variant("p.Lys215_Gly217del")
    expect_equal(result$type, "deletion")
    expect_equal(result$pos, 215)
    expect_equal(result$end_pos, 217)
})

# --- Validation edge cases ---

test_that("validate_hgvs rejects identical ref and alt", {
    result <- validate_hgvs("NM_000546.6:c.215A>A")
    expect_false(result$semantic_valid)
})

test_that("validate_hgvs rejects non-standard accession prefix", {
    result <- validate_hgvs("XM_000546.6:c.215C>G")
    expect_false(result$syntax_valid)
})
