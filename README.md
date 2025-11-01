# Liger - a smart concatenation tool to create supermatrices 
## Overview
Liger is a fast, composable, smart concatenation tool that can process an unlimited number of FASTA files and create a supermatrix with accompanying partitions. The key feature of Liger is smart matching of FASTA headers. Users provide a list of taxa to include with canonical headers (i.e., the headers you wish to use for the final supermatrix). Liger takes this list and matches headers across all input files. This allows users to preserve metadata without having to manually alter headers in the input files.

## Performance

Liger is built for speed. Benchmarked against FASconCAT v1.11 on a dataset of 13 genes, 61 taxa, 43,363 bp:

| Tool | Time |
|------|------|
| Liger | 95ms |
| FASconCAT (Perl) | 967ms |

**Liger is ~10x faster** while maintaining the same functionality plus fuzzy matching.

**Test System:**
- CPU: Intel Ultra 5 125U (14 cores) @ 4.3GHz
- RAM: 16GB
- OS: Ubuntu 24.04.3 LTS

## Installation and Dependencies

Liger written is written in go. Go 1.21 or higher is reqiured. No other dependencies are required. It can be downloaded and complied as follows:
```
# Download the file
wget https://raw.githubusercontent.com/andrewbudge/Liger/refs/heads/main/liger.go

# Compile it
go build -o liger liger.go

# Move to PATH (optional)
sudo mv liger /usr/local/bin/
```

## Usage
```bash
liger [INPUT FASTA FILES] [TAXA LIST] > [SUPERMATRIX] 2> [PARTITONS]
```
- **Input**: Pre-aligned gene files + taxa list (one name per line)
- **Output**: Supermatrix in FASTA format (stdout) + NEXUS partitions (stderr)
- **Missing taxa**: Automatically filled with N's

## How Matching Works

Liger searches for your taxon name anywhere in the FASTA header.

### Basic Example
```bash
$ cat COX1.fasta
>AB123.1 Mus mus COX1 gene, partial cds
ATCGATCGATCG
>AB124.1 Rattus rat COX1 gene, partial cds
GCTAGCTAGCTA
>AB125.1 Ovis sheep COX1 gene, partial cds
CGATCGATCGAT

$ cat ND2.fasta
>XM456.1 Mus mus ND2 gene, complete cds
TACGTACGTACG
>XM457.1 Rattus rat ND2 gene, complete cds
ATATATATATAT
>XM458.1 Ovis sheep ND2 gene, complete cds
GCGCGCGCGCGC

$ cat taxa.txt
Mus mus
Rattus rat
Ovis sheep

$ liger COX1.fasta ND2.fasta taxa.txt > matrix.fasta 2> parts.nex

$ cat matrix.fasta
>Mus mus
ATCGATCGATCGTACGTACGTACG
>Rattus rat
GCTAGCTAGCTAATATATATATAT
>Ovis sheep
CGATCGATCGATGCGCGCGCGCGC

$ cat parts.nex
#NEXUS
begin sets;
  charset COX1 = 1-12;
  charset ND2 = 13-24;
end;
```

### Missing Data
```bash
$ cat COX1.fasta
>AB123.1 Mus mus COX1
ATCGATCGATCG
>AB124.1 Rattus rat COX1
GCTAGCTAGCTA

$ cat ND2.fasta
>XM456.1 Mus mus ND2
TACGTACGTACG
# Rattus rat missing!

$ cat taxa.txt
Mus mus
Rattus rat

$ liger COX1.fasta ND2.fasta taxa.txt > matrix.fasta 2>/dev/null

$ cat matrix.fasta
>Mus mus
ATCGATCGATCGTACGTACGTACG
>Rattus rat
GCTAGCTAGCTANNNNNNNNNNNN
```

