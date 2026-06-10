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
