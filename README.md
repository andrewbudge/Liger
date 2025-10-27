# Liger - a smart concatenation tool to create supermatrices 
## Overview
Liger is a smart concatenation tool that can process an unlimited number of FASTA files and create a supermatrix with accompanying partitions. The key feature of Liger is smart matching of FASTA headers. Users provide a list of taxa to include with canonical headers (i.e., the headers you wish to use for the final supermatrix). Liger takes this list and matches headers across all input files. This allows users to preserve metadata without having to manually alter headers in the input files.

## Installation and Dependencies

The bash version of liger is depend on Seqkit. Be sure to have it installed before using Liger. Seqkit is open source and available [here](https://github.com/shenwei356/seqkit). 
```bash
# Download liger
wget https://raw.githubusercontent.com/andrewbudge/Liger/main/liger
chmod +x liger
sudo mv liger /usr/local/bin/

# Verify
liger
```
There is also a faster version of liger written in go. It can be downloaded an complied as follows:
```
# Download the file
wget https://raw.githubusercontent.com/andrewbudge/Liger/main/liger.go

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

**WARNING: Liger will not check if input files are aligned!!! Be sure to check files before putting them into Liger.** 

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

### Test Matching
```bash
$ seqkit grep -r -n -p "Mus mus" COX1.fasta
>AB123.1 Mus mus COX1 gene, partial cds
ATCGATCGATCG
```    
Should return exactly one sequence per file.
