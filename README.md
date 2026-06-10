# rvarsim

R/Bioconductor variant simulator with HGVS notation and MANE Select transcripts.

## Features

- **Variant simulation**: Generate all possible SNVs across CDS, UTR, and
  canonical splice sites for MANE Select transcripts
- **HGVS parsing**: Parse c./g./p. notation into structured R objects
- **Validation**: Syntactic and semantic HGVS validation
- **Normalization**: 3' shifting, canonical representation
- **Format conversion**: HGVS ↔ VCF ↔ SPDI
- **Transcription mapping**: g. ↔ c. coordinate mapping via TxDb
- **Translation**: c. → p. protein consequence prediction
- **Backtranslation**: p. → c. codon enumeration
- **Variant extraction**: Sequence alignment → HGVS descriptions
- **Liftover**: Assembly mapping via chain files

## Installation

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rvarsim")
```

## Quick Start

```r
library(rvarsim)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# Parse an HGVS variant
parse_hgvs("NM_000546.6:c.215C>G")[[1]]

# Validate
is_valid_hgvs("NM_000546.6:c.215C>G")

# Convert to VCF
hgvs_to_vcf("NC_000001.11:g.123456A>G")

# Simulate all SNVs for a transcript
result <- simulate_variants(
    txdb     = EnsDb.Hsapiens.v86,
    bsgenome = BSgenome.Hsapiens.UCSC.hg38,
    transcript_ids = "ENST00000357654"
)
```

## License

MIT. See [LICENSE](LICENSE).
