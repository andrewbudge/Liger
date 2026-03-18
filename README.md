<div align="center">
  <img src="liger_logo_test.svg" width="400">
</div>

# Liger - a fast, smart supermatrix concatenation tool

## Overview

Liger is a fast, composable concatenation tool that creates supermatrices from multiple FASTA alignments. It runs in two modes:

- **Exact match (default):** headers must match exactly across files, like FASconCAT and AMAS.
- **Smart match (`-a alias.txt`):** pass an alias list of clean output names that get matched to messy input headers via case-insensitive substring search. Underscores in aliases match spaces in headers, so `Mus_musculus` finds `AB123.1 Mus musculus COX1 gene, partial cds`. Longer aliases match first to prevent partial collisions. The alias list doubles as a rename map — input headers stay messy, output gets clean names. Requires `-l` for a provenance TSV that records exactly which original header matched each alias.

Liger auto-detects DNA vs amino acid data per gene and adjusts missing characters and partition labels accordingly. FASTA output goes to stdout, partition boundaries to stderr in RAxML/IQ-TREE format by default. NEXUS bundles everything into one file.

Liger is also available as the `concat` subcommand in [phylo](https://github.com/andrewbudge/phylo).

## Performance

| Scale | Taxa x Genes | Liger | FASconCAT-G | Speedup |
|---|---|---|---|---|
| Small | 100 x 20 | 19ms | 12s | **637x** |
| Medium | 300 x 50 | 146ms | 4 min | **1,646x** |
| Large | 1,000 x 30 | 225ms | 5.4 min | **1,438x** |
| Mega | 4,000 x 30 | 968ms | crash | - |

**Test System:**
- CPU: Intel Ultra 5 125U (14 cores) @ 4.3GHz
- RAM: 16GB
- OS: Ubuntu 24.04.3 LTS

## Install

Requires [Rust](https://www.rust-lang.org/tools/install).

```bash
cargo install --git https://github.com/andrewbudge/Liger
```

To update to the latest version:

```bash
cargo install --force --git https://github.com/andrewbudge/Liger
```

The previous Go version (v1) is preserved in the `v1-go/` directory.

## Usage

```bash
liger [FLAGS] [INPUT FASTA FILES]
```

### Flags
- `-a, --alias` — alias list for smart matching (clean output names that map to messy input headers)
- `-l, --log` — provenance TSV output file (required with `-a`)
- `-f, --format` — output format: fasta (default), nexus (also accepts `n` or `nex`)
- `-m, --missing` — override missing data character (default: auto per data type — N for DNA, X for amino acid, ? for mixed)
- `-p, --partitions` — partition format: raxml (default, also used by IQ-TREE) or nexus

### Exact match — clean headers

```bash
$ liger gene1.fasta gene2.fasta > supermatrix.fasta
DNA, gene1.fasta = 1-500
DNA, gene2.fasta = 501-1000
```

### Smart match — messy headers with an alias list

```bash
$ cat alias.txt
Mus_musculus
Rattus_rattus
Xenopus_laevis

$ liger -a alias.txt -l prov.tsv gene1.fasta gene2.fasta > supermatrix.fasta
DNA, gene1.fasta = 1-4
DNA, gene2.fasta = 5-8

$ cat supermatrix.fasta
>Mus_musculus
ATCGATCG
>Rattus_rattus
ATCGNNNN
>Xenopus_laevis
NNNNATCG

$ cat prov.tsv
alias.txt	gene1.fasta	gene2.fasta
Mus_musculus	AB123.1 Mus musculus gene1 cds	XM456.1 Mus musculus gene2 cds
Rattus_rattus	AB124.1 Rattus rattus gene1 cds	MISSING
Xenopus_laevis	MISSING	XM789.1 Xenopus laevis gene2 cds
```

### NEXUS output

```bash
$ liger -a alias.txt -l prov.tsv -f nexus gene1.fasta gene2.fasta
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=N GAP=-;
  MATRIX
  Mus_musculus    ATCGATCG
  Rattus_rattus   ATCGNNNN
  Xenopus_laevis  NNNNATCG
;
END;
BEGIN SETS;
  CHARSET gene1.fasta = 1-4;
  CHARSET gene2.fasta = 5-8;
END;
```

### Partition formats

Partitions always go to stderr so they can be piped or redirected.

```bash
# Default: RAxML/IQ-TREE format
$ liger dna_gene.fasta protein_gene.fasta > matrix.fasta
DNA, dna_gene.fasta = 1-500
WAG, protein_gene.fasta = 501-700

# NEXUS charset format
$ liger -p nexus dna_gene.fasta protein_gene.fasta > matrix.fasta
CHARSET dna_gene.fasta = 1-500;
CHARSET protein_gene.fasta = 501-700;
```

## Development Note

Liger v2 is a rewrite in Rust of the original Go tool (v1). Development is assisted by Claude (Anthropic), which serves as a teaching aid and coding partner. The design, domain knowledge, and direction are the author's own.

## Author

Andrew Budge
