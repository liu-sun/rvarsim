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
