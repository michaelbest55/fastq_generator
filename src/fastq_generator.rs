use clap::{Arg, Command};
use rand::seq::SliceRandom;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process;

/// Represents where output should be written. Supports stdout or writing to a file.
enum OutputDest {
    /// Wraps a buffered writer for stdout.
    Stdout(BufWriter<io::Stdout>),
    /// Wraps a buffered writer for a file.
    File(BufWriter<File>),
}

impl Write for OutputDest {
    /// Writes the given buffer into the underlying writer.
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputDest::Stdout(writer) => writer.write(buf),
            OutputDest::File(writer) => writer.write(buf),
        }
    }

    /// Flushes the underlying writer.
    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputDest::Stdout(writer) => writer.flush(),
            OutputDest::File(writer) => writer.flush(),
        }
    }
}

impl OutputDest {
    /// Creates a new `OutputDest` from an optional file path.
    ///
    /// # Arguments
    ///
    /// * `path` - An optional `PathBuf`. If `Some`, the output will be written to that file;
    ///            otherwise, output is written to stdout.
    ///
    /// # Returns
    ///
    /// Returns an `io::Result<Self>` representing either the created `OutputDest` or an error.
    fn new(path: Option<PathBuf>) -> io::Result<Self> {
        match path {
            Some(path) => Ok(OutputDest::File(BufWriter::new(File::create(path)?))),
            None => Ok(OutputDest::Stdout(BufWriter::new(io::stdout()))),
        }
    }

    /// Writes an entire buffer into the underlying writer.
    ///
    /// # Arguments
    ///
    /// * `buf` - A byte slice representing the data to be written.
    ///
    /// # Returns
    ///
    /// Returns an `io::Result<()>` indicating success or an error.
    #[allow(dead_code)]
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            OutputDest::Stdout(writer) => writer.write_all(buf),
            OutputDest::File(writer) => writer.write_all(buf),
        }
    }
}

/// Returns the reverse complement of a DNA sequence.
///
/// This function converts the input sequence to uppercase, maps each nucleotide to its complement,
/// and then reverses the result.
///
/// # Arguments
///
/// * `sequence` - A string slice containing the DNA sequence.
///
/// # Examples
///
/// ```
/// let rc = reverse_complement("ATCG");
/// assert_eq!(rc, "CGAT");
/// ```
///
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

/// Processes input to generate reverse complements for each DNA sequence.
///
/// Reads from either a file (if provided) or stdin, computes the reverse complement for
/// lines that do not start with '>', and writes the results to output (file or stdout).
///
/// # Arguments
///
/// * `input` - An `Option<PathBuf>` specifying the input file path; if `None`, reads from stdin.
/// * `output` - An `Option<PathBuf>` specifying the output file path; if `None`, writes to stdout.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or an error.
///
fn process_reverse_complement(input: Option<PathBuf>, output: Option<PathBuf>) -> io::Result<()> {
    let input_reader: Box<dyn BufRead> = match input {
        Some(path) => Box::new(BufReader::new(File::open(path)?)),
        None => Box::new(BufReader::new(io::stdin())),
    };

    let mut output_dest = OutputDest::new(output)?;

    for line in input_reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(output_dest, "{}", line)?;
        } else {
            let rc = reverse_complement(line.trim());
            writeln!(output_dest, "{}", rc)?;
        }
    }

    output_dest.flush()?;
    Ok(())
}

/// Generates a random DNA sequence of a specified length.
///
/// # Arguments
///
/// * `sequence_size` - The desired length of the DNA sequence.
///
/// # Returns
///
/// Returns a `String` containing the random DNA sequence.
///
fn generate_dna(sequence_size: usize) -> String {
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut rng = rand::thread_rng();

    (0..sequence_size)
        .map(|_| *nucleotides.choose(&mut rng).unwrap())
        .collect()
}

/// Generates and writes a FASTA formatted file with random DNA sequences.
///
/// For each sequence, writes a header (`>seq_N`) and a generated DNA sequence.
///
/// # Arguments
///
/// * `sequence_size` - Length of each DNA sequence.
/// * `nb_seq` - Number of sequences to generate.
/// * `output` - An `OutputDest` where the FASTA content is written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_fasta(sequence_size: usize, nb_seq: usize, mut output: OutputDest) -> io::Result<()> {
    for i in 1..=nb_seq {
        writeln!(output, ">seq_{}", i)?;
        writeln!(output, "{}", generate_dna(sequence_size))?;
    }
    output.flush()
}

/// Generates and writes a single FASTQ formatted sequence.
///
/// Writes a header, the DNA sequence, a plus sign, and a quality string (defaulting to 'I').
///
/// # Arguments
///
/// * `sequence` - The DNA sequence as a string slice.
/// * `header` - The header for the sequence (includes '@').
/// * `output` - A mutable reference to an `OutputDest` for writing the sequence.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_fastq_seq(sequence: &str, header: &str, output: &mut OutputDest) -> io::Result<()> {
    writeln!(output, "{}", header)?;
    writeln!(output, "{}", sequence)?;
    writeln!(output, "+")?;
    writeln!(output, "{}", "I".repeat(sequence.len()))?;
    Ok(())
}

/// Generates and writes multiple FASTQ formatted sequences with random DNA sequences.
///
/// # Arguments
///
/// * `sequence_size` - Length of each DNA sequence.
/// * `nb_seq` - Number of sequences to generate.
/// * `output` - An `OutputDest` where the FASTQ data will be written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_fastq(sequence_size: usize, nb_seq: usize, mut output: OutputDest) -> io::Result<()> {
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!("@SEQ_ID_{}", index);
        generate_fastq_seq(&sequence, &header, &mut output)?;
    }
    output.flush()
}

/// Generates and writes random single-end FASTQ sequences.
///
/// Each sequence gets a formatted header and a randomly generated DNA sequence.
///
/// # Arguments
///
/// * `sequence_size` - Length of each DNA sequence.
/// * `nb_seq` - Number of sequences to generate.
/// * `output` - An `OutputDest` where the FASTQ data will be written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_random_fastq_se(
    sequence_size: usize,
    nb_seq: usize,
    mut output: OutputDest,
) -> io::Result<()> {
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!(
            "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{} 1:N:0:TGACCAAT",
            index
        );
        generate_fastq_seq(&sequence, &header, &mut output)?;
    }
    output.flush()
}

/// Generates and writes random paired-end FASTQ sequences.
///
/// First writes the first read for each sequence, then the corresponding second read.
///
/// # Arguments
///
/// * `sequence_size` - Length of each DNA sequence.
/// * `nb_seq` - Number of sequence pairs to generate.
/// * `output` - An `OutputDest` where the FASTQ data will be written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_random_fastq_pe(
    sequence_size: usize,
    nb_seq: usize,
    mut output: OutputDest,
) -> io::Result<()> {
    // Generate first reads
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!("@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/1", index);
        generate_fastq_seq(&sequence, &header, &mut output)?;
    }

    // Generate second reads
    for index in 1..=nb_seq {
        let sequence = generate_dna(sequence_size);
        let header = format!("@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/2", index);
        generate_fastq_seq(&sequence, &header, &mut output)?;
    }
    output.flush()
}

/// Generates and writes mapped single-end FASTQ sequences from a reference FASTA file.
///
/// Reads sequences from the reference file and prints sequence fragments with formatted headers.
/// The sequences are broken into kmers of fixed length based on the provided coverage.
///
/// # Arguments
///
/// * `ref_fasta` - Path to the reference FASTA file as a string slice.
/// * `sequence_size` - Length of each generated kmer.
/// * `coverage` - The desired coverage level (number of passes over each sequence).
/// * `output` - An `OutputDest` where the FASTQ data will be written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_mapped_fastq_se(
    ref_fasta: &str,
    sequence_size: usize,
    coverage: usize,
    mut output: OutputDest,
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

    let mut index = 1;

    // Generate reads for each coverage level
    for coverage_level in 1..=coverage {
        for sequence in &sequences {
            let mut start = coverage_level - 1;

            while start + sequence_size < sequence.len() {
                let kmer = &sequence[start..start + sequence_size];

                writeln!(
                    output,
                    "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{} 1:N:0:ATGT",
                    index
                )?;
                writeln!(output, "{}", kmer)?;
                writeln!(output, "+")?;
                writeln!(output, "{}", "I".repeat(sequence_size))?;

                start += sequence_size;
                index += 1;
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Generates and writes mapped paired-end FASTQ sequences from a reference FASTA file.
///
/// Processes the reference and collects kmers to simulate paired-end reads.
/// For each valid region, computes a forward read and the reverse complement of a corresponding region,
/// then writes both to the output.
///
/// # Arguments
///
/// * `ref_fasta` - Path to the reference FASTA file as a string slice.
/// * `sequence_size` - Length of each kmer.
/// * `insertion_size` - Size of the insertion between paired-end reads.
/// * `coverage` - The desired sequencing coverage.
/// * `output` - An `OutputDest` where the FASTQ data will be written.
///
/// # Returns
///
/// Returns an `io::Result<()>` indicating success or error.
///
fn generate_mapped_fastq_pe(
    ref_fasta: &str,
    sequence_size: usize,
    insertion_size: usize,
    coverage: usize,
    mut output: OutputDest,
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

    // First pass: generate and store all kmers
    let mut kmers = Vec::new();

    for coverage_level in 1..=coverage {
        for sequence in &sequences {
            let mut start = coverage_level - 1;

            while start + sequence_size * 2 + insertion_size < sequence.len() {
                // Forward read
                let kmer1 = &sequence[start..start + sequence_size];

                // Reverse read
                let kmer2_start = start + sequence_size + insertion_size;
                let kmer2 = &sequence[kmer2_start..kmer2_start + sequence_size];
                let kmer2_rc = reverse_complement(kmer2);

                kmers.push((kmer1.to_string(), kmer2_rc));

                start += (sequence_size * 2) + insertion_size;
            }
        }
    }

    // Write forward reads
    for (i, (kmer1, _)) in kmers.iter().enumerate() {
        writeln!(output, "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/1", i + 1)?;
        writeln!(output, "{}", kmer1)?;
        writeln!(output, "+")?;
        writeln!(output, "{}", "I".repeat(sequence_size))?;
    }

    // Write reverse reads
    for (i, (_, kmer2)) in kmers.iter().enumerate() {
        writeln!(output, "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{}#0/2", i + 1)?;
        writeln!(output, "{}", kmer2)?;
        writeln!(output, "+")?;
        writeln!(output, "{}", "I".repeat(sequence_size))?;
    }

    output.flush()?;
    Ok(())
}

/// Prints an error message to stderr and exits the process.
///
/// # Arguments
///
/// * `msg` - The error message to display.
///
/// # Note
///
/// This function does not return.
///
fn exit_with_error(msg: &str) -> ! {
    eprintln!("Error: {}", msg);
    process::exit(1);
}

/// Adds the sequence size argument to a Clap command.
///
/// # Arguments
///
/// * `cmd` - A `clap::Command` to which the argument is added.
///
/// # Returns
///
/// Returns the modified `Command`.
///
fn add_sequence_size_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("sequence_size")
            .long("sequence-size")
            .short('s')
            .required(true)
            .value_parser(clap::value_parser!(usize))
            .help("Length of each generated sequence"),
    );
    cmd
}

/// Adds the number of sequences argument to a Clap command.
///
/// # Arguments
///
/// * `cmd` - A `clap::Command` to which the argument is added.
///
/// # Returns
///
/// Returns the modified `Command`.
///
fn add_nb_seq_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("nb_seq")
            .required(true)
            .long("nb_seq")
            .short('n')
            .value_parser(clap::value_parser!(usize))
            .help("Number of sequences to generate"),
    );
    cmd
}

/// Adds the coverage argument to a Clap command.
///
/// # Arguments
///
/// * `cmd` - A `clap::Command` to which the argument is added.
///
/// # Returns
///
/// Returns the modified `Command`.
///
fn add_coverage_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("coverage")
            .required(true)
            .long("coverage")
            .short('c')
            .value_parser(clap::value_parser!(usize))
            .help("Desired coverage"),
    );
    cmd
}

/// Adds the reference FASTA argument to a Clap command.
///
/// # Arguments
///
/// * `cmd` - A `clap::Command` to which the argument is added.
///
/// # Returns
///
/// Returns the modified `Command`.
///
fn add_ref_fasta_arg(mut cmd: Command) -> Command {
    cmd = cmd.arg(
        Arg::new("ref_fasta")
            .required(true)
            .long("ref_fasta")
            .short('r')
            .help("Reference FASTA file"),
    );
    cmd
}

/// Adds the output file argument to a Clap command.
///
/// # Arguments
///
/// * `cmd` - A `clap::Command` to which the argument is added.
///
/// # Returns
///
/// Returns the modified `Command`.
///
fn add_output_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("output")
            .long("output")
            .short('o')
            .value_name("FILE")
            .help("Output file path (defaults to stdout)")
            .value_parser(clap::value_parser!(PathBuf)),
    )
}

/// Main entry point for the application.
///
/// Sets up the CLI with subcommands and arguments, then executes the chosen functionality.
///
/// # Returns
///
/// Returns a `Result<(), Box<dyn Error>>` indicating success or error.
///
fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("dna_tools")
        .version("1.0")
        .author("Michael Best")
        .about("DNA sequence manipulation and FASTQ generation tools")
        .subcommand(add_output_arg(
            Command::new("reverse_complement")
                .arg(
                    Arg::new("input")
                        .long("input")
                        .short('i')
                        .value_name("FILE")
                        .help("Input file path (defaults to stdin)")
                        .value_parser(clap::value_parser!(PathBuf)),
                )
                .about("Generates reverse complement of DNA sequences from stdin"),
        ))
        .subcommand(add_output_arg(add_sequence_size_arg(add_nb_seq_arg(
            Command::new("generate_fasta").about("Creates a de novo fasta file"),
        ))))
        .subcommand(add_output_arg(add_nb_seq_arg(add_sequence_size_arg(
            Command::new("generate_fastq").about("Creates a de novo fastq file"),
        ))))
        .subcommand(add_output_arg(add_nb_seq_arg(add_sequence_size_arg(
            Command::new("generate_random_fastq_se").about("Generates random single-end reads"),
        ))))
        .subcommand(add_output_arg(add_nb_seq_arg(add_sequence_size_arg(
            Command::new("generate_random_fastq_pe").about("Generates random paired-end reads"),
        ))))
        .subcommand(add_output_arg(add_coverage_arg(add_sequence_size_arg(
            add_ref_fasta_arg(
                Command::new("generate_mapped_fastq_se").about("Generates mapped single-end reads"),
            ),
        ))))
        .subcommand(add_output_arg(add_coverage_arg(add_sequence_size_arg(
            add_ref_fasta_arg(
                Command::new("generate_mapped_fastq_pe")
                    .about("Generates mapped paired-end reads")
                    .arg(
                        Arg::new("insertion_size")
                            .long("insertion-size")
                            .short('i')
                            .required(true)
                            .value_parser(clap::value_parser!(usize))
                            .help("Size of insertion"),
                    ),
            ),
        ))))
        .get_matches();

    match matches.subcommand() {
        Some(("reverse_complement", sub_m)) => {
            let input_path = sub_m.get_one::<PathBuf>("input").cloned();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            if let Err(e) = process_reverse_complement(input_path, output_path) {
                exit_with_error(&format!("Error processing input: {}", e));
            }
        }
        Some(("generate_fasta", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) = generate_fasta(sequence_size, nb_seq, output) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_fastq", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) = generate_fastq(sequence_size, nb_seq, output) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_random_fastq_se", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) = generate_random_fastq_se(sequence_size, nb_seq, output) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_random_fastq_pe", sub_m)) => {
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let nb_seq = *sub_m.get_one::<usize>("nb_seq").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) = generate_random_fastq_pe(sequence_size, nb_seq, output) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_mapped_fastq_se", sub_m)) => {
            let ref_fasta = sub_m.get_one::<String>("ref_fasta").unwrap();
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let coverage = *sub_m.get_one::<usize>("coverage").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) = generate_mapped_fastq_se(ref_fasta, sequence_size, coverage, output) {
                exit_with_error(&e.to_string());
            }
        }
        Some(("generate_mapped_fastq_pe", sub_m)) => {
            let ref_fasta = sub_m.get_one::<String>("ref_fasta").unwrap();
            let sequence_size = *sub_m.get_one::<usize>("sequence_size").unwrap();
            let insertion_size = *sub_m.get_one::<usize>("insertion_size").unwrap();
            let coverage = *sub_m.get_one::<usize>("coverage").unwrap();
            let output_path = sub_m.get_one::<PathBuf>("output").cloned();
            let output = OutputDest::new(output_path)?;
            if let Err(e) =
                generate_mapped_fastq_pe(ref_fasta, sequence_size, insertion_size, coverage, output)
            {
                exit_with_error(&e.to_string());
            }
        }
        _ => {
            exit_with_error("Please specify a valid subcommand. Use --help for usage information.")
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::{fs, io::Write};

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
    fn test_generate_fasta_stdout() {
        let sequence_size = 10;
        let nb_seq = 10;
        let output = OutputDest::new(None).unwrap();
        generate_fasta(sequence_size, nb_seq, output).unwrap();
    }

    #[test]
    fn test_generate_fasta_to_file() {
        let sequence_size = 10;
        let nb_seq = 10;

        // Create a named temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_path.clone())).unwrap();
        generate_fasta(sequence_size, nb_seq, output).unwrap();

        assert!(temp_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_path).unwrap();

        // Add assertions to check the content
        assert!(content.starts_with(">"), "FASTA file should start with >");
        assert_eq!(
            content.lines().count(),
            nb_seq * 2,
            "Incorrect number of lines in FASTA file"
        );
    }

    #[test]
    fn test_generate_fastq_stdout() {
        let sequence_size = 10;
        let nb_seq = 10;
        let output = OutputDest::new(None).unwrap();
        generate_fastq(sequence_size, nb_seq, output).unwrap();
    }

    #[test]
    fn test_generate_fastq_to_file() {
        let sequence_size = 10;
        let nb_seq = 10;

        // Create a named temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_path.clone())).unwrap();
        generate_fastq(sequence_size, nb_seq, output).unwrap();
        assert!(temp_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_path).unwrap();

        // Add assertions to check the content
        assert!(content.starts_with("@"), "FASTQ file should start with @");
        assert_eq!(
            content.lines().count(),
            nb_seq * 4,
            "Incorrect number of lines in FASTQ file"
        );
    }

    #[test]
    fn test_generate_random_fastq_se_stdout() {
        let sequence_size = 10;
        let nb_seq = 10;
        let output = OutputDest::new(None).unwrap();
        generate_random_fastq_se(sequence_size, nb_seq, output).unwrap();
    }

    #[test]
    fn test_generate_random_fastq_se_to_file() {
        let sequence_size = 10;
        let nb_seq = 10;

        // Create a named temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_path.clone())).unwrap();
        generate_random_fastq_se(sequence_size, nb_seq, output).unwrap();
        assert!(temp_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_path).unwrap();

        // Add assertions to check the content
        assert!(content.starts_with("@"), "FASTQ file should start with @");
        assert_eq!(
            content.lines().count(),
            nb_seq * 4,
            "Incorrect number of lines in FASTQ file"
        );
    }

    #[test]
    fn test_generate_random_fastq_pe_stdout() {
        let sequence_size = 10;
        let nb_seq = 10;
        let output = OutputDest::new(None).unwrap();
        generate_random_fastq_pe(sequence_size, nb_seq, output).unwrap();
    }

    #[test]
    fn test_generate_random_fastq_pe_to_file() {
        let sequence_size = 10;
        let nb_seq = 10;

        // Create a named temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_path.clone())).unwrap();
        generate_random_fastq_pe(sequence_size, nb_seq, output).unwrap();
        assert!(temp_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_path).unwrap();

        // Add assertions to check the content
        assert!(content.starts_with("@"), "FASTQ file should start with @");
        assert_eq!(
            content.lines().count(),
            nb_seq * 8,
            "Incorrect number of lines in FASTQ file"
        );
    }

    // fn create_test_fasta(temp_file: &mut NamedTempFile) -> io::Result<()> {
    //     let fasta_content = ">seq1\nACGTACGTACGT\n\

    //                         >seq2\nTGCATGCATGCA\n\
    //                         >seq3\nGCTAGCTAGCTA\n\
    //                         >seq4\nCGATCGATCGAT\n\
    //                         >seq5\nTACGTACGTACG\n\
    //                         >seq6\nGTACGTACGTAC\n\
    //                         >seq7\nCTAGCTAGCTAG\n\
    //                         >seq8\nATCGATCGATCG\n\
    //                         >seq9\nCGTACGTACGTA\n\
    //                         >seq10\nGATCGATCGATC\n";

    //     // Write content to provided temporary file
    //     temp_file
    //         .as_file_mut()
    //         .write_all(fasta_content.as_bytes())?;

    //     Ok(())
    // }
    fn create_test_fasta(
        temp_file: &mut NamedTempFile,
        num_seq: usize,
        seq_length: usize,
    ) -> io::Result<()> {
        let nucleotides = ['A', 'C', 'G', 'T'];
        let mut rng = rand::thread_rng();
        let mut fasta_content = String::new();

        for i in 1..=num_seq {
            let sequence: String = (0..seq_length)
                .map(|_| *nucleotides.choose(&mut rng).unwrap())
                .collect();
            fasta_content.push_str(&format!(">seq{}\n{}\n", i, sequence));
        }

        // Write content to provided temporary file
        temp_file
            .as_file_mut()
            .write_all(fasta_content.as_bytes())?;

        Ok(())
    }

    #[test]
    fn test_generate_mapped_fastq_se_stdout() {
        // Create the temporary file
        let mut temp_file = NamedTempFile::new().unwrap();

        // Populate the file with test data
        create_test_fasta(&mut temp_file, 10, 12).unwrap();

        let sequence_size = 10;
        let coverage = 2;

        let output = OutputDest::new(None).unwrap();
        generate_mapped_fastq_se(
            temp_file.path().to_str().unwrap(),
            sequence_size,
            coverage,
            output,
        )
        .unwrap();
    }

    #[test]
    fn test_generate_mapped_fastq_se_to_file() {
        // Create the temporary file
        let mut temp_file = NamedTempFile::new().unwrap();

        // Populate the file with test data
        create_test_fasta(&mut temp_file, 10, 12).unwrap();

        let sequence_size = 10;
        let coverage = 2;

        // Create a named temporary file for output
        let temp_output_file = NamedTempFile::new().unwrap();
        let temp_output_path = temp_output_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_output_path.clone())).unwrap();
        generate_mapped_fastq_se(
            temp_file.path().to_str().unwrap(),
            sequence_size,
            coverage,
            output,
        )
        .unwrap();
        assert!(temp_output_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_output_path).unwrap();

        // Add assertions to check the content
        assert!(content.starts_with("@"), "FASTQ file should start with @");
    }

    #[test]
    fn test_generate_mapped_fastq_pe_stdout() {
        // Create the temporary file
        let mut temp_file = NamedTempFile::new().unwrap();

        // Populate the file with test data
        create_test_fasta(&mut temp_file, 10, 12).unwrap();

        let sequence_size = 10;
        let insertion_size = 5;
        let coverage = 2;

        let output = OutputDest::new(None).unwrap();
        generate_mapped_fastq_pe(
            temp_file.path().to_str().unwrap(),
            sequence_size,
            insertion_size,
            coverage,
            output,
        )
        .unwrap();
    }

    #[test]
    fn test_generate_mapped_fastq_pe_to_file() {
        // Create the temporary file
        let mut temp_file = NamedTempFile::new().unwrap();

        // Populate the file with test data
        create_test_fasta(&mut temp_file, 10, 36).unwrap();

        let sequence_size = 10;
        let insertion_size = 5;
        let coverage = 2;

        // Create a named temporary file for output
        let temp_output_file = NamedTempFile::new().unwrap();
        let temp_output_path = temp_output_file.path().to_owned();

        // Create OutputDest with the temporary file path
        let output = OutputDest::new(Some(temp_output_path.clone())).unwrap();
        generate_mapped_fastq_pe(
            temp_file.path().to_str().unwrap(),
            sequence_size,
            insertion_size,
            coverage,
            output,
        )
        .unwrap();
        assert!(temp_output_path.exists(), "Temporary file was not created");

        // Read the content of the file
        let content = fs::read_to_string(&temp_output_path).unwrap();
        // Add assertions to check the content
        assert!(content.starts_with("@"), "FASTQ file should start with @");
    }
}
