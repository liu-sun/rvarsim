test_that("HGVS CDS notation is correct", {
    # Test CDS variant: c.215C>G
    hgvs <- .build_hgvs_c(
        region = "cds", cds_pos = 215L,
        exon_boundary = NA_integer_, offset = NA_integer_,
        ref = "C", alt = "G",
        transcript = "NM_000546.6", use_version = FALSE
    )
    expect_equal(hgvs, "NM_000546.6:c.215C>G")
})

test_that("HGVS 5'UTR notation is correct", {
    hgvs <- .build_hgvs_c(
        region = "five_utr", cds_pos = -14L,
        exon_boundary = NA_integer_, offset = NA_integer_,
        ref = "G", alt = "A",
        transcript = "NM_000314.8", use_version = FALSE
    )
    expect_equal(hgvs, "NM_000314.8:c.-14G>A")
})

test_that("HGVS 3'UTR notation uses asterisk", {
    hgvs <- .build_hgvs_c(
        region = "three_utr", cds_pos = 32L,
        exon_boundary = NA_integer_, offset = NA_integer_,
        ref = "T", alt = "C",
        transcript = "NM_001126114.2", use_version = FALSE
    )
    expect_equal(hgvs, "NM_001126114.2:c.*32T>C")
})

test_that("HGVS splice donor notation is correct", {
    hgvs <- .build_hgvs_c(
        region = "splice_donor", cds_pos = NA_integer_,
        exon_boundary = 453L, offset = 1L,
        ref = "G", alt = "T",
        transcript = "NM_000546.6", use_version = FALSE
    )
    expect_equal(hgvs, "NM_000546.6:c.453+1G>T")
})

test_that("HGVS splice acceptor notation is correct", {
    hgvs <- .build_hgvs_c(
        region = "splice_acceptor", cds_pos = NA_integer_,
        exon_boundary = 612L, offset = -2L,
        ref = "A", alt = "G",
        transcript = "NM_000546.6", use_version = FALSE
    )
    expect_equal(hgvs, "NM_000546.6:c.612-2A>G")
})

test_that("HGVS genomic notation is correct", {
    hgvs <- .build_hgvs_g("NC_000001.11", 123456L, "A", "G")
    expect_equal(hgvs, "NC_000001.11:g.123456A>G")
})

test_that("NA ref/alt returns NA HGVS", {
    expect_true(is.na(.build_hgvs_c(
        "cds", 10L, NA_integer_, NA_integer_, NA_character_, "G",
        "TX", FALSE
    )))
})

test_that("format_hgvs handles empty input", {
    empty <- data.frame(
        transcript_id   = character(0),
        region          = character(0),
        genomic_pos     = integer(0),
        genomic_ref     = character(0),
        genomic_alt     = character(0),
        cds_pos         = integer(0),
        exon_boundary   = integer(0),
        offset          = integer(0),
        stringsAsFactors = FALSE
    )
    result <- format_hgvs(empty)
    expect_equal(nrow(result), 0)
    expect_true("hgvs_c" %in% colnames(result))
    expect_true("hgvs_g" %in% colnames(result))
})
