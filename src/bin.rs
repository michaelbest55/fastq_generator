use clap::{Arg, Command};
use fastq_generator::*;
use std::error::Error;
use std::path::PathBuf;
use std::process; // Import everything from your library

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
pub fn main() -> Result<(), Box<dyn Error>> {
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
