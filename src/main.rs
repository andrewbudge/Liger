use clap::Parser;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
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
        let line = line.trim().to_string();

        if line.starts_with('>') {
            if !current_header.is_empty() {
                current_seq = current_seq.to_uppercase();
                sequences.insert(current_header.clone(), current_seq.clone());
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
            }
            current_header = line[1..].to_string();
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.push_str(&line);
        }
    }

    if !current_header.is_empty() {
        current_seq = current_seq.to_uppercase();
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

    // Parse all gene files
    let mut gene_data = Vec::new();
    for file in &args.files {
        let (sequences, length) = parse_fasta(file, true).expect("Failed to parse fasta file");
        gene_data.push((file, sequences, length));
    }

    let smart_matching = args.alias.is_some();

    // Determine taxa list and build matched_genes
    let taxa: Vec<String>;
    let mut matched_genes: Vec<(&String, HashMap<String, (String, String)>, &usize)> = Vec::new();

    if let Some(ref taxa_file) = args.alias {
        // Smart matching mode: load taxa list, match by substring
        taxa = load_taxa_list(taxa_file).expect("Failed to load taxa list");
        for (file, sequences, length) in &gene_data {
            let matched = match_taxa(&taxa, sequences);
            matched_genes.push((file, matched, length));
        }
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
        // Direct lookup — header IS the taxon name, no matching needed
        for (file, sequences, length) in &gene_data {
            let mut matched = HashMap::new();
            for (header, sequence) in sequences {
                matched.insert(header.clone(), (header.clone(), sequence.clone()));
            }
            matched_genes.push((file, matched, length));
        }
    }

    // Detect data type per gene using first available sequence
    let gene_types: Vec<DataType> = matched_genes
        .iter()
        .map(|(_, matched, _)| {
            matched
                .values()
                .next()
                .map(|(_, seq)| detect_data_type(seq))
                .unwrap_or(DataType::Dna)
        })
        .collect();

    // Determine overall data type mix
    let has_dna = gene_types.iter().any(|t| matches!(t, DataType::Dna));
    let has_aa = gene_types.iter().any(|t| matches!(t, DataType::AminoAcid));
    let is_mixed = has_dna && has_aa;

    // Build supermatrix: concatenate matched sequences per taxon, fill gaps with missing char
    let mut supermatrix: HashMap<String, String> = HashMap::new();
    for taxon in &taxa {
        for (i, (_file, matched, length)) in matched_genes.iter().enumerate() {
            let entry = supermatrix.entry(taxon.clone()).or_insert(String::new());
            if matched.contains_key(taxon) {
                entry.push_str(&matched[taxon].1);
            } else {
                let missing_char = match &args.missing {
                    Some(m) => m.clone(),
                    None => if is_mixed {
                        "?".to_string()
                    } else {
                        match gene_types[i] {
                            DataType::Dna => "N".to_string(),
                            DataType::AminoAcid => "X".to_string(),
                        }
                    },
                };
                entry.push_str(&missing_char.repeat(**length));
            }
        }
    }

    // Build partition boundaries with data type
    let mut partitions = Vec::new();
    let mut position = 1;
    for (i, (file, _matched, length)) in matched_genes.iter().enumerate() {
        let name = Path::new(file).file_name().unwrap().to_str().unwrap();
        partitions.push((name.to_string(), position, position + *length - 1, &gene_types[i]));
        position = position + *length;
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

    let fmt = args.format.to_lowercase();
    if fmt == "nexus" || fmt == "n" || fmt == "nex" {
        // NEXUS: complete file to stdout (data + partitions in one)
        println!("#NEXUS");
        println!("BEGIN DATA;");
        println!("  DIMENSIONS NTAX={} NCHAR={};", taxa.len(), total_length);
        println!("  FORMAT DATATYPE={} MISSING={} GAP=-;", nexus_datatype, nexus_missing);
        println!("  MATRIX");
        for taxon in &taxa {
            println!("  {}    {}", taxon, supermatrix[taxon]);
        }
        println!(";");
        println!("END;");
        println!("BEGIN SETS;");
        for (name, start, end, _) in &partitions {
            println!("  CHARSET {} = {}-{};", name, start, end);
        }
        println!("END;");
    } else {
        // FASTA: supermatrix to stdout
        for taxon in &taxa {
            println!(">{}", taxon);
            println!("{}", supermatrix[taxon]);
        }
    }

    // Partitions to stderr
    let part_fmt = args.partitions.to_lowercase();
    if part_fmt == "nexus" || part_fmt == "n" || part_fmt == "nex" {
        for (name, start, end, _) in &partitions {
            eprintln!("CHARSET {} = {}-{};", name, start, end);
        }
    } else {
        // RAxML/IQ-TREE format (default)
        for (name, start, end, dtype) in &partitions {
            let model = match dtype {
                DataType::Dna => "DNA",
                DataType::AminoAcid => "WAG",
            };
            eprintln!("{}, {} = {}-{}", model, name, start, end);
        }
    }

    // Write provenance TSV (only in smart matching mode, -l is required with -a)
    if smart_matching {
        let log_path = args.log.as_ref().unwrap();
        let taxa_file = args.alias.as_ref().unwrap();
        let mut log_file = File::create(log_path).expect("Failed to create provenance log file");
        // Header row: taxa list filename, then each gene filename
        let gene_names: Vec<String> = matched_genes
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
            for (_file, matched, _length) in &matched_genes {
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
