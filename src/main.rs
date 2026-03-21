use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// Exit silently on broken pipe (e.g., piping to head/tail)
#[cfg(unix)]
fn reset_sigpipe() {
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

#[derive(Parser)]
#[command(name = "liger", version, about = "Fast, smart supermatrix concatenation for phylogenetics")]
struct Args {
    /// FASTA alignment files
    pub files: Vec<String>,

    /// Alias list for smart matching — output names that map to messy input headers
    #[arg(short, long, requires = "log")]
    pub alias: Option<String>,

    /// Output format (FASTA or Nexus)
    #[arg(short, long, default_value = "FASTA")]
    pub format: String,

    /// Override missing data character (default: N for DNA, X for amino acid)
    #[arg(short, long)]
    pub missing: Option<String>,

    /// Partition format: raxml (default) or nexus
    #[arg(short = 'p', long = "partitions", default_value = "raxml")]
    pub partitions: String,

    /// Provenance TSV output file (required with -a)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,

    /// Number of threads (default: all available cores)
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,
}

enum DataType {
    Dna,
    AminoAcid,
}

fn detect_data_type(sequence: &str) -> DataType {
    if sequence
        .chars()
        .any(|c| matches!(c, 'E' | 'F' | 'I' | 'L' | 'P' | 'Q'))
    {
        DataType::AminoAcid
    } else {
        DataType::Dna
    }
}

fn parse_fasta(
    filename: &str,
    validate_equal: bool,
) -> Result<(HashMap<String, String>, usize), String> {
    let file = File::open(filename).map_err(|e| format!("Could not open {}: {}", filename, e))?;
    let reader = BufReader::new(file);

    let mut sequences = HashMap::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();
    let mut expected_length: Option<usize> = None;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading file: {}", e))?;
        let line = line.trim();

        if line.starts_with('>') {
            if !current_header.is_empty() {
                current_seq.make_ascii_uppercase();
                if validate_equal {
                    match expected_length {
                        None => expected_length = Some(current_seq.len()),
                        Some(len) => {
                            if current_seq.len() != len {
                                return Err(format!(
                                    "Error: Sequence length mismatch in {} : {}",
                                    filename, current_header
                                ));
                            }
                        }
                    }
                }
                sequences.insert(
                    std::mem::take(&mut current_header),
                    std::mem::take(&mut current_seq),
                );
            }
            current_header = line[1..].to_string();
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }

    if !current_header.is_empty() {
        current_seq.make_ascii_uppercase();
        if validate_equal {
            if let Some(len) = expected_length {
                if current_seq.len() != len {
                    return Err(format!(
                        "Error: Sequence length mismatch in {} : {}",
                        filename, current_header
                    ));
                }
            }
        }
        sequences.insert(current_header, current_seq);
    }

    let length = sequences.values().next().map_or(0, |s| s.len());
    Ok((sequences, length))
}

fn load_taxa_list(filename: &str) -> Result<Vec<String>, String> {
    let file = File::open(filename).map_err(|e| format!("Could not open {}: {}", filename, e))?;
    let reader = BufReader::new(file);

    let mut taxa = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading file: {}", e))?;
        let line = line.trim().to_string();
        if !line.is_empty() {
            taxa.push(line);
        }
    }

    Ok(taxa)
}

/// Match taxa names to FASTA headers using case-insensitive substring search.
/// Longer taxa names match first to prevent partial match collisions.
/// Once a header is claimed, no other taxon can match it.
fn match_taxa(
    taxa: &[String],
    sequences: &HashMap<String, String>,
) -> HashMap<String, (String, String)> {
    let mut sorted_taxa = taxa.to_vec();
    sorted_taxa.sort_by(|a, b| b.len().cmp(&a.len()));

    let mut claimed_headers = HashSet::new();
    let mut results = HashMap::new();

    for taxon in &sorted_taxa {
        for (header, sequence) in sequences {
            if claimed_headers.contains(header) {
                continue;
            }
            if header
                .to_lowercase()
                .contains(&taxon.to_lowercase().replace("_", " "))
            {
                claimed_headers.insert(header.clone());
                results.insert(taxon.clone(), (header.clone(), sequence.clone()));
                break;
            }
        }
    }
    results
}

fn main() {
    #[cfg(unix)]
    reset_sigpipe();

    let args = Args::parse();

    // Configure thread pool
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to set thread count");
    }

    // Parse all gene files in parallel
    let gene_data: Vec<_> = args
        .files
        .par_iter()
        .map(|file| {
            let (sequences, length) = parse_fasta(file, true).expect("Failed to parse fasta file");
            (file, sequences, length)
        })
        .collect();

    let smart_matching = args.alias.is_some();

    // Determine taxa list
    let taxa: Vec<String>;

    // Smart matching needs an intermediate matched_genes structure;
    // exact match can look up directly from gene_data.
    let matched_genes: Option<Vec<(&String, HashMap<String, (String, String)>, &usize)>>;

    if let Some(ref taxa_file) = args.alias {
        // Smart matching mode: load taxa list, match by substring
        taxa = load_taxa_list(taxa_file).expect("Failed to load taxa list");
        matched_genes = Some(
            gene_data
                .par_iter()
                .map(|(file, sequences, length)| {
                    let matched = match_taxa(&taxa, sequences);
                    (*file, matched, length)
                })
                .collect(),
        );
    } else {
        // Exact match mode: union of all headers across files becomes the taxa list
        let mut seen = HashSet::new();
        let mut taxa_order = Vec::new();
        for (_file, sequences, _length) in &gene_data {
            for header in sequences.keys() {
                if seen.insert(header.clone()) {
                    taxa_order.push(header.clone());
                }
            }
        }
        taxa = taxa_order;
        matched_genes = None;
    }

    // Detect data type per gene using first available sequence
    let gene_types: Vec<DataType> = gene_data
        .iter()
        .map(|(_, sequences, _)| {
            sequences
                .values()
                .next()
                .map(|seq| detect_data_type(seq))
                .unwrap_or(DataType::Dna)
        })
        .collect();

    // Determine overall data type mix
    let has_dna = gene_types.iter().any(|t| matches!(t, DataType::Dna));
    let has_aa = gene_types.iter().any(|t| matches!(t, DataType::AminoAcid));
    let is_mixed = has_dna && has_aa;

    // Pre-compute missing fill strings per gene
    let missing_chars: Vec<String> = gene_types
        .iter()
        .zip(gene_data.iter())
        .map(|(dtype, (_, _, length))| {
            let ch = match &args.missing {
                Some(m) => m.as_str(),
                None => {
                    if is_mixed {
                        "?"
                    } else {
                        match dtype {
                            DataType::Dna => "N",
                            DataType::AminoAcid => "X",
                        }
                    }
                }
            };
            ch.repeat(*length)
        })
        .collect();

    // Build supermatrix in parallel: one entry per taxon, indexed by taxa order
    let supermatrix: Vec<String> = if let Some(ref mg) = matched_genes {
        // Smart matching: look up via matched_genes
        taxa.par_iter()
            .map(|taxon| {
                let mut seq = String::new();
                for (i, (_file, matched, _length)) in mg.iter().enumerate() {
                    if let Some((_, s)) = matched.get(taxon) {
                        seq.push_str(s);
                    } else {
                        seq.push_str(&missing_chars[i]);
                    }
                }
                seq
            })
            .collect()
    } else {
        // Exact match: look up directly from gene_data — no clones needed
        taxa.par_iter()
            .map(|taxon| {
                let mut seq = String::new();
                for (i, (_file, sequences, _length)) in gene_data.iter().enumerate() {
                    if let Some(s) = sequences.get(taxon) {
                        seq.push_str(s);
                    } else {
                        seq.push_str(&missing_chars[i]);
                    }
                }
                seq
            })
            .collect()
    };

    // Build partition boundaries with data type
    let mut partitions = Vec::new();
    let mut position = 1;
    for (i, (file, _sequences, length)) in gene_data.iter().enumerate() {
        let name = Path::new(file).file_name().unwrap().to_str().unwrap();
        partitions.push((name.to_string(), position, position + length - 1, &gene_types[i]));
        position = position + length;
    }
    let total_length = position - 1;

    // Determine overall data type label for NEXUS output
    let nexus_datatype = match (has_dna, has_aa) {
        (true, true) => "MIXED",
        (false, true) => "PROTEIN",
        _ => "DNA",
    };

    // Default missing char for NEXUS format line
    let nexus_missing = match &args.missing {
        Some(m) => m.clone(),
        None => match (has_dna, has_aa) {
            (true, true) => "?".to_string(),
            (false, true) => "X".to_string(),
            _ => "N".to_string(),
        },
    };

    // Buffered output
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let fmt = args.format.to_lowercase();
    if fmt == "nexus" || fmt == "n" || fmt == "nex" {
        // NEXUS: complete file to stdout (data + partitions in one)
        writeln!(out, "#NEXUS").unwrap();
        writeln!(out, "BEGIN DATA;").unwrap();
        writeln!(out, "  DIMENSIONS NTAX={} NCHAR={};", taxa.len(), total_length).unwrap();
        writeln!(
            out,
            "  FORMAT DATATYPE={} MISSING={} GAP=-;",
            nexus_datatype, nexus_missing
        )
        .unwrap();
        writeln!(out, "  MATRIX").unwrap();
        for (i, taxon) in taxa.iter().enumerate() {
            writeln!(out, "  {}    {}", taxon, supermatrix[i]).unwrap();
        }
        writeln!(out, ";").unwrap();
        writeln!(out, "END;").unwrap();
        writeln!(out, "BEGIN SETS;").unwrap();
        for (name, start, end, _) in &partitions {
            writeln!(out, "  CHARSET {} = {}-{};", name, start, end).unwrap();
        }
        writeln!(out, "END;").unwrap();
    } else {
        // FASTA: supermatrix to stdout
        for (i, taxon) in taxa.iter().enumerate() {
            writeln!(out, ">{}", taxon).unwrap();
            writeln!(out, "{}", supermatrix[i]).unwrap();
        }
    }

    // Partitions to stderr
    let stderr = std::io::stderr();
    let mut err = BufWriter::new(stderr.lock());
    let part_fmt = args.partitions.to_lowercase();
    if part_fmt == "nexus" || part_fmt == "n" || part_fmt == "nex" {
        for (name, start, end, _) in &partitions {
            writeln!(err, "CHARSET {} = {}-{};", name, start, end).unwrap();
        }
    } else {
        // RAxML/IQ-TREE format (default)
        for (name, start, end, dtype) in &partitions {
            let model = match dtype {
                DataType::Dna => "DNA",
                DataType::AminoAcid => "WAG",
            };
            writeln!(err, "{}, {} = {}-{}", model, name, start, end).unwrap();
        }
    }

    // Write provenance TSV (only in smart matching mode, -l is required with -a)
    if smart_matching {
        let mg = matched_genes.as_ref().unwrap();
        let log_path = args.log.as_ref().unwrap();
        let taxa_file = args.alias.as_ref().unwrap();
        let mut log_file = File::create(log_path).expect("Failed to create provenance log file");
        // Header row: taxa list filename, then each gene filename
        let gene_names: Vec<String> = mg
            .iter()
            .map(|(file, _, _)| {
                Path::new(file)
                    .file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .to_string()
            })
            .collect();
        writeln!(log_file, "{}\t{}", taxa_file, gene_names.join("\t"))
            .expect("Failed to write to log file");
        // One row per taxon: taxon name, then matched header or MISSING
        for taxon in &taxa {
            let mut row = vec![taxon.clone()];
            for (_file, matched, _length) in mg {
                if matched.contains_key(taxon) {
                    row.push(matched[taxon].0.clone());
                } else {
                    row.push("MISSING".to_string());
                }
            }
            writeln!(log_file, "{}", row.join("\t")).expect("Failed to write to log file");
        }
    }
}
