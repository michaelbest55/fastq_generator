use clap::{Command, Arg};
use rand::seq::SliceRandom;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;


/// Takes a DNA sequence and returns its reverse complement
/// This function is used both for the reverse_complement subcommand
/// and within the FASTQ generation functions
fn reverse_complement(sequence: &str) -> String {
    let sequence = sequence.to_uppercase();
    
    let complement: String = sequence
        .chars()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect();
    
    complement.chars().rev().collect()
}

/// Processes stdin for the reverse_complement subcommand
fn process_reverse_complement() -> io::Result<()> {
    let stdin = io::stdin();
    let handle = stdin.lock();

    for line in handle.lines() {
        let line = line?;
        if line.starts_with('>') {
            println!("{}", line);
        } else {
            let rc = reverse_complement(line.trim());
            println!("{}", rc);
        }
    }

    Ok(())
}

/// Generates a random DNA sequence of specified length
fn generate_dna(sequence_size: usize) -> String {
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut rng = rand::thread_rng();
    
    (0..sequence_size)
        .map(|_| *nucleotides.choose(&mut rng).unwrap())
        .collect()
}


/// Generates and prints a FASTA format sequence
fn generate_fasta(sequence_size: usize, nb_seq: usize) {
    for i in 1..=nb_seq {
        println!(">seq_{}\n{}", i, generate_dna(sequence_size));
    }
}

/// Generates and prints a FASTQ format sequence
fn generate_fastq(sequence: &str, header: &str) {
    println!("{}", header);
    println!("{}", sequence);
    println!("+");
    // Generate quality string of same length as sequence, using 'I' for best quality
    println!("{}", "I".repeat(sequence.len()));
}

/// Generates random single-end FASTQ sequences
fn generate_random_fastq_se(sequence_size: usize, nb_seq: usize) {
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{} 1:N:0:TGACCAAT",
            index
        );
        generate_fastq(&sequence, &header);
    }
}

/// Generates random paired-end FASTQ sequences
fn generate_random_fastq_pe(sequence_size: usize, nb_seq: usize) {
    // Generate first reads
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/1",
            index
        );
        generate_fastq(&sequence, &header);
    }
    
    // Generate second reads
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/2",
            index
        );
        generate_fastq(&sequence, &header);
    }
}

/// Generates mapped single-end FASTQ sequences from a reference
fn generate_mapped_fastq_se(ref_fasta: &str, sequence_size: usize, coverage: usize) -> io::Result<()> {
    let file = File::open(ref_fasta)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    
    // Read sequences from FASTA file
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            sequences.push(line.trim().to_string());
        }
    }
    
    let mut index = 1;
    
    // Generate reads for each coverage level
    for coverage_level in 1..=coverage {
        for sequence in &sequences {
            let mut start = coverage_level - 1;
            while start + sequence_size < sequence.len() {
                let kmer = &sequence[start..start + sequence_size];
                let header = format!(
                    "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{} 1:N:0:ATGT",
                    index
                );
                generate_fastq(kmer, &header);
                start += sequence_size;
                index += 1;
            }
        }
    }
    
    Ok(())
}

/// Generates mapped paired-end FASTQ sequences from a reference
fn generate_mapped_fastq_pe(
    ref_fasta: &str,
    sequence_size: usize,
    insertion_size: usize,
    coverage: usize,
) -> io::Result<()> {
    let file = File::open(ref_fasta)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    
    // Read sequences from FASTA file
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            sequences.push(line.trim().to_string());
        }
    }
    
    let mut kmers1 = Vec::new();
    let mut kmers2 = Vec::new();
    
    // Generate paired reads for each coverage level
    for coverage_level in 1..=coverage {
        for sequence in &sequences {
            let mut start = coverage_level - 1;
            while start + sequence_size * 2 + insertion_size < sequence.len() {
                let kmer1 = &sequence[start..start + sequence_size];
                kmers1.push(kmer1.to_string());
                
                let kmer2_start = start + insertion_size + sequence_size;
                let kmer2 = &sequence[kmer2_start..kmer2_start + sequence_size];
                kmers2.push(reverse_complement(kmer2));
                
                start += (sequence_size * 2) + insertion_size;
            }
        }
    }
    
    // Print first reads
    for (index, kmer1) in kmers1.iter().enumerate() {
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/1",
            index + 1
        );
        generate_fastq(kmer1, &header);
    }
    
    // Print second reads
    for (index, kmer2) in kmers2.iter().enumerate() {
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/2",
            index + 1
        );
        generate_fastq(kmer2, &header);
    }
    
    Ok(())
}

fn exit_with_error(msg: &str) -> ! {
    eprintln!("Error: {}", msg);
    process::exit(1);
}

// Helper function to add common arguments
fn add_sequence_size_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("sequence_size")
            .long("sequence-size")
            .short('s')
            .required(true)
            .value_parser(clap::value_parser!(usize))
            .help("Length of each generated sequence")
    );
    cmd
}

fn add_nb_seq_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("nb_seq")
            .required(true)
            .long("nb_seq")
            .short('n')
            .value_parser(clap::value_parser!(usize))
            .help("Number of sequences to generate")
    );
    cmd
}

fn add_coverage_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("coverage")
            .required(true)
            .value_parser(clap::value_parser!(usize))
            .help("Desired coverage")
    );
    cmd
}

fn add_ref_fasta_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("ref_fasta")
            .required(true)
            .help("Reference FASTA file")
    );
    cmd
}

fn main() {
    let matches = Command::new("dna_tools")
        .version("1.0")
        .author("Michael Best")
        .about("DNA sequence manipulation and FASTQ generation tools")
        .subcommand(
            Command::new("reverse_complement")
                .about("Generates reverse complement of DNA sequences from stdin")
        )
        .subcommand(
            add_sequence_size_arg(
                add_nb_seq_arg(
                    Command::new("generate_fasta")
                        .about("Creates a de novo fasta file")
                )
            )
        )
        .subcommand(
            add_nb_seq_arg(
                add_sequence_size_arg(
                    Command::new("generate_random_fastq_se")
                        .about("Generates random single-end reads")
                )
            )
        )
        .subcommand(
            add_nb_seq_arg(
                add_sequence_size_arg(
                    Command::new("generate_random_fastq_pe")
                        .about("Generates random paired-end reads")
                )
            )
        )
        .subcommand(
            add_coverage_arg(
                add_sequence_size_arg(
                    add_ref_fasta_arg(
                        Command::new("generate_mapped_fastq_se")
                            .about("Generates mapped single-end reads")
                    )
                )
            )
        )
        .subcommand(
            add_coverage_arg(
                add_sequence_size_arg(
                    add_ref_fasta_arg(
                        Command::new("generate_mapped_fastq_pe")
                            .about("Generates mapped paired-end reads")
                            .arg(
                                Arg::new("insertion_size")
                                    .long("insertion-size")
                                    .short('i')
                                    .required(true)
                                    .value_parser(clap::value_parser!(usize))
                                    .help("Size of insertion")
                            )
                    )
                )
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("reverse_complement", _)) => {
            if let Err(e) = process_reverse_complement() {
                exit_with_error(&format!("Error processing input: {}", e));
            }
        }
        Some(("generate_fasta", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            generate_fasta(sequence_size, nb_seq);
        }
        Some(("generate_random_fastq_se", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            generate_random_fastq_se(sequence_size, nb_seq);
        }
        Some(("generate_random_fastq_pe", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            generate_random_fastq_pe(sequence_size, nb_seq);
        }
        Some(("generate_mapped_fastq_se", sub_m)) => {
            let ref_fasta = sub_m.get_one::<String>("ref_fasta").unwrap();
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let coverage = *sub_m.get_one::<usize>("coverage").unwrap();
            
            if let Err(e) = generate_mapped_fastq_se(ref_fasta, sequence_size, coverage) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_mapped_fastq_pe", sub_m)) => {
            let ref_fasta = sub_m.get_one::<String>("ref_fasta").unwrap();
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let insertion_size = *sub_m.get_one::<usize>("insertion_size").unwrap();
            let coverage = *sub_m.get_one::<usize>("coverage").unwrap();
            
            if let Err(e) = generate_mapped_fastq_pe(ref_fasta, sequence_size, insertion_size, coverage) {
                exit_with_error(&e.to_string());
            }
        }
        _ => exit_with_error("Please specify a valid subcommand. Use --help for usage information."),
    }
}

#[cfg(test)]
mod tests {
    use std::{fs, path::Path, io::Write};

    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACCTTGG"), "CCAAGGT");
        assert_eq!(reverse_complement("accttgg"), "CCAAGGT");
        assert_eq!(reverse_complement(""), "");
        assert_eq!(reverse_complement("A"), "T");
        assert_eq!(reverse_complement("T"), "A");
        assert_eq!(reverse_complement("C"), "G");
        assert_eq!(reverse_complement("G"), "C");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("GGCC"), "GGCC");
    }

    #[test]
    fn test_generate_dna() {
        let dna = generate_dna(10);
        assert_eq!(dna.len(), 10);
        assert!(dna.chars().all(|c| ['A', 'C', 'G', 'T'].contains(&c)));

        let dna = generate_dna(0);
        assert_eq!(dna.len(), 0);

        let dna = generate_dna(100);
        assert_eq!(dna.len(), 100);
        assert!(dna.chars().all(|c| ['A', 'C', 'G', 'T'].contains(&c)));
    }

    #[test]
    fn test_generate_fasta() {
        let sequence_size = 10;
        let nb_seq = 10;
        generate_fasta(sequence_size, nb_seq);
        // This test is more about checking the output manually
    }

    #[test]
    fn test_generate_fastq() {
        let sequence = "ACGTACGTAC";
        let header = "@test_header";
        generate_fastq(sequence, header);
        // This test is more about checking the output manually
    }

    #[test]
    fn test_generate_random_fastq_se() {
        generate_random_fastq_se(10, 5);
        // This test is more about checking the output manually
    }

    #[test]
    fn test_generate_random_fastq_pe() {
        generate_random_fastq_pe(10, 5);
        // This test is more about checking the output manually
    }


fn create_test_fasta(temp_file: &mut NamedTempFile) -> io::Result<()> {
    let fasta_content = ">seq1\nACGTACGTACGT\n\
                        >seq2\nTGCATGCATGCA\n\
                        >seq3\nGCTAGCTAGCTA\n\
                        >seq4\nCGATCGATCGAT\n\
                        >seq5\nTACGTACGTACG\n\
                        >seq6\nGTACGTACGTAC\n\
                        >seq7\nCTAGCTAGCTAG\n\
                        >seq8\nATCGATCGATCG\n\
                        >seq9\nCGTACGTACGTA\n\
                        >seq10\nGATCGATCGATC\n";

    // Write content to provided temporary file
    temp_file.as_file_mut().write_all(fasta_content.as_bytes())?;

    Ok(())
}

#[test]
fn test_generate_mapped_fastq_se() {
    // Create directory for test data
    let test_dir = Path::new("test_data");
    fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
    // Create the temporary file
    let mut temp_file = NamedTempFile::new_in(test_dir)
        .expect("Failed to create temp file");
    
    // Populate the file with test data
    create_test_fasta(&mut temp_file)
        .expect("Failed to write test data");

    let sequence_size = 10;
    let coverage = 2;

    let result = generate_mapped_fastq_se(temp_file.path().to_str().unwrap(), sequence_size, coverage);
    assert!(result.is_ok());
}

    #[test]
    fn test_generate_mapped_fastq_pe() {
        // Create directory for test data
        let test_dir = Path::new("test_data");
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
        
        // Create the temporary file
        let mut temp_file = NamedTempFile::new_in(test_dir)
            .expect("Failed to create temp file");
        
        // Populate the file with test data
        create_test_fasta(&mut temp_file)
            .expect("Failed to write test data");
        let sequence_size = 10;
        let insertion_size = 5;
        let coverage = 2;
        let result = generate_mapped_fastq_pe(temp_file.path().to_str().unwrap(), sequence_size, insertion_size, coverage);
        assert!(result.is_ok());
    }

}

