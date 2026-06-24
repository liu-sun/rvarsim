# rvarsim 1.0.0

* Bioconductor release addressing full review feedback
* Vignette rewrite: motivation, package comparison table, terminology
  definitions (MANE, HGVS, SPDI), Bioconductor installation instructions,
  interoperability section, evaluated chunks, references
* Fixed DESCRIPTION: removed "R/Bioconductor" from Title, matched LICENSE
  (MIT) with DESCRIPTION, added SequenceMatching biocView
* Fixed fetch_mane.R: deprecated TxidFilter → TxIdFilter, removed
  hardcoded tx_support_level column, NULL-guarded filter construction
* Fixed hgvs_extract.R sprintf format mismatch (5 fields, 4 args)
* Removed redundant requireNamespace(BSgenome) — already in Imports
* Moved GenomicRanges from Imports to Depends
* Refactored switch() statements to reduce repeated argument passing
* Added acronym expansions (MANE, HGVS, SPDI) in function documentation
* Added 30+ new tests (normalization, alignment, conversion edge cases,
  format_hgvs integration, parse edge cases, protein parser, validation)

# rvarsim 0.99.1

* Resubmission addressing Bioconductor review feedback
* Replaced \dontrun{} with \donttest{} in all examples
* Removed <<- super-assignment patterns for cleaner scoping
* Added proper S3 method registration in NAMESPACE
* Fixed line lengths to 80 character limit
* Removed duplicate biocViews entry
* Added package-level documentation

# rvarsim 0.99.0

* Initial Bioconductor submission
* Variant simulation pipeline: fetch MANE Select transcripts, generate
  all possible SNVs across CDS, UTRs, and splice sites
* HGVS parsing and validation (c./g./p. notation)
* HGVS normalization: 3' shifting, common affix trimming
* Format conversion: HGVS ↔ VCF ↔ SPDI
* Transcription mapping: c. ↔ g. via TxDb
* Translation: c. → p. consequence prediction (missense, nonsense,
  frameshift, silent)
* Backtranslation: p. → c. codon enumeration
* Variant extraction from pairwise sequence alignment
* Liftover between genome assemblies via chain files
