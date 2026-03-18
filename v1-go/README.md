<div align="center">
  <img src="liger_logo_test.svg" width="400">
</div>

# Liger - a smart concatenation tool to create supermatrices 
## Overview
Liger is a fast, composable, smart concatenation tool that can process an unlimited number of FASTA files and create a supermatrix with accompanying partitions. The key feature of Liger is smart matching of FASTA headers. Users provide a list of taxa to include with canonical headers (i.e., the headers you wish to use for the final supermatrix). Liger takes this list and matches headers across all input files. This allows users to preserve metadata without having to alter headers in the input files.

## Performance

Liger is built for speed. Benchmarked against FASconCAT v1.11 on a dataset of 13 genes, 61 taxa, 43,363 bp:

| Tool | Time |
|------|------|
| Liger | ~45ms |
| FASconCAT (Perl) | ~965ms |

**Liger is ~20x faster** while maintaining the same functionality plus fuzzy matching.

**Test System:**
- CPU: Intel Ultra 5 125U (14 cores) @ 4.3GHz
- RAM: 16GB
- OS: Ubuntu 24.04.3 LTS

## Installation and Dependencies

Liger written is written in go. Go 1.21 or higher is required. No other dependencies are required. It can be downloaded and complied as follows:
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
liger [FLAGS] [TAXA LIST] [INPUT FASTA FILES] > [OUTPUT]
```

### Flags
- `-f` Output format: `fasta` (default) or `nexus`
- `-m` Character for missing data (default: `N`)

### Input/Output
- **Input**: Taxa list (one name per line) + pre-aligned gene files
- **Output (FASTA mode)**: Supermatrix to stdout, NEXUS partitions to stderr
- **Output (NEXUS mode)**: Complete NEXUS file to stdout (includes supermatrix + partitions)
- **Missing taxa**: Automatically filled with specified missing character

## How Matching Works

Liger searches for your taxon name anywhere in the FASTA header.

### Basic Example (FASTA output)
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

$ liger taxa.txt COX1.fasta ND2.fasta > matrix.fasta 2> parts.nex

$ cat matrix.fasta
>Mus mus
ATCGATCGATCGTACGTACGTACG
>Ovis sheep
CGATCGATCGATGCGCGCGCGCGC
>Rattus rat
GCTAGCTAGCTAATATATATATAT

$ cat parts.nex
#NEXUS
begin sets;
  charset COX1 = 1-12;
  charset ND2 = 13-24;
end;
```

### NEXUS Output Example
```bash
$ liger -f nexus taxa.txt COX1.fasta ND2.fasta > output.nex

$ cat output.nex
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=3;
    TAXLABELS
        Mus mus
        Ovis sheep
        Rattus rat
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=24;
    FORMAT DATATYPE=DNA MISSING=N GAP=-;
    MATRIX
        Mus mus     ATCGATCGATCGTACGTACGTACG
        Ovis sheep  CGATCGATCGATGCGCGCGCGCGC
        Rattus rat  GCTAGCTAGCTAATATATATATAT
    ;
END;

BEGIN SETS;
    CHARSET COX1 = 1-12;
    CHARSET ND2 = 13-24;
END;
```

### Missing Data Example
```bash
$ cat COX1.fasta
>AB123.1 Mus mus COX1
ATCGATCGATCG
>AB124.1 Rattus rat COX1
GCTAGCTAGCTA

$ cat ND2.fasta
>XM456.1 Mus mus ND2
TACGTACGTACG
# Ovis sheep missing from ND2!

$ cat taxa.txt
Mus mus
Rattus rat
Ovis sheep

$ liger taxa.txt COX1.fasta ND2.fasta 2>/dev/null
>Mus mus
ATCGATCGATCGTACGTACGTACG
>Ovis sheep
CGATCGATCGATNNNNNNNNNNNN
>Rattus rat
GCTAGCTAGCTAATATATATATAT
```

### Custom Missing Character
```bash
$ liger -m ? taxa.txt COX1.fasta ND2.fasta 2>/dev/null
>Mus mus
ATCGATCGATCGTACGTACGTACG
>Ovis sheep
CGATCGATCGAT????????????
>Rattus rat
GCTAGCTAGCTAATATATATATAT
```

