# `fastq-generator`

![Rust Version][rustc-image]
[![crates.io][crate-image]][crate-link]
[![Documentation][docs-image]][docs-link]
[![Dependency Status][deps-image]][deps-link]

A simple cli tool for generating fastq files.

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [`fastq-generator`](#project-name)
- [Installation](#installation)
    - [Compile from Source](#compile-from-source)
- [Usage](#usage)
    - [Command Line Interface](#command-line-interface)
- [License](#license)
    - [Contribution](#contribution)

<!-- markdown-toc end -->

This project started at a time when I needed to generate a large number of fastq files quickly and with a large number of reads. I wanted a simple tool that could generate fastq files with a large number of reads quickly. I also wanted to learn
more about Rust and how to write a simple CLI tool in Rust. This project is the result of that. As this is my first time writing in Rust, I am sure there are many things that could be done better. I am open to suggestions and PRs. 

# Installation

```
cargo install fastq-generator
```

`fastq-generator` is a single binary that must be placed somewhere in your
`$PATH`.

## Compile from Source

Ensure you have a [Rust toolchain installed](https://rustup.rs). Some of the
dependencies also require `gcc` to be installed.

```
$ git clone https://github.com/michaelbest55/fastq-generator
$ cd fastq-generator
$ cargo build --release
$ sudo cp target/release/fastq-generator /usr/local/bin/
```

# Usage

```
# Reverse complement
fastq_generator reverse_complement --input path/to/input.fasta --output path/to/output.fasta

# Generate a FASTA file
fastq_generator generate_fasta --sequence-size 100 --nb-seq 10 --output path/to/output.fasta

# Generate random single-end FASTQ
fastq_generator generate_random_fastq_se --sequence-size 50 --nb-seq 5 --output path/to/output.fastq

# Generate random paired-end FASTQ
fastq_generator generate_random_fastq_pe --sequence-size 50 --nb-seq 5 --output path/to/output.fastq

# Generate mapped single-end FASTQ
fastq_generator generate_mapped_fastq_se --ref-fasta path/to/reference.fasta --sequence-size 75 --coverage 10 --output path/to/output.fastq

# Generate mapped paired-end FASTQ
fastq_generator generate_mapped_fastq_pe --ref-fasta path/to/reference.fasta --sequence-size 75 --coverage 10 --insertion-size 300 --output path/to/output.fastq
```

# License

This crate is licensed under [MIT license](http://opensource.org/licenses/MIT).

[//]: # (badges)

[rustc-image]: https://img.shields.io/badge/rustc-1.53+-blue.svg
[crate-image]: https://img.shields.io/crates/v/fastq-generator.svg
[crate-link]: https://crates.io/crates/fastq-generator
[docs-image]: https://docs.rs/fastq-generator/badge.svg
[docs-link]: https://docs.rs/fastq-generator
[deps-image]: https://deps.rs/repo/github/michaelbest55/fastq-generator/status.svg
[deps-link]: https://deps.rs/repo/github/michaelbest55/fastq-generator
