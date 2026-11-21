# Sequence Handling – FASTA Sequence Analysis

## Project Overview

This project is a Python-based bioinformatics tool designed to read and analyze nucleotide sequences stored in FASTA format.

The script automatically parses FASTA files and calculates basic sequence properties, including sequence length, GC content, and nucleotide composition.

## Objectives

- Parse FASTA files using Python.
- Calculate the length of each sequence.
- Determine GC content as a percentage.
- Count individual nucleotide bases.
- Identify ambiguous nucleotide symbols.
- Generate CSV and text-based analysis reports.

## Dataset

- Organism: *Escherichia coli* K-12 MG1655
- Database: NCBI
- Accession: NC_000913.3
- File format: FASTA

## Technologies Used

- Python 3
- CSV module
- Collections module
- Argparse module
- pathlib module

## Project Structure

```text
SequenceHandling/
│
├── data/
│   └── ecoli_k12.fasta
│
├── results/
│   ├── sequence_statistics.csv
│   └── summary_report.txt
│
├── analyze_fasta.py
└── README.md
```

## Features

### 1. FASTA Parsing
Reads sequence identifiers and nucleotide sequences automatically from a FASTA file.

### 2. Sequence Length
Calculates the total number of characters in each sequence.

### 3. GC Content
Calculates the percentage of guanine and cytosine bases among canonical nucleotide bases.

### 4. Nucleotide Composition
Reports the counts of A, T, G, C and other supported nucleotide symbols.

### 5. Automated Reports
Generates:
- `sequence_statistics.csv`
- `summary_report.txt`

## How to Run

Run the following command from the project directory:

```bash
python analyze_fasta.py data/ecoli_k12.fasta --outdir results
```

To analyze another FASTA file:

```bash
python analyze_fasta.py data/another_file.fasta --outdir results
```

## GC Content Formula

GC content is calculated using:

GC Content (%) = ((G + C) / (A + T + G + C + U)) × 100

Ambiguous nucleotide symbols and gaps are excluded from the denominator.

## Expected Results

The program produces individual sequence statistics and an overall summary containing:

- Number of sequences
- Total sequence length
- Overall GC content
- Nucleotide composition

## Applications

This tool demonstrates basic sequence handling and nucleotide analysis used in bioinformatics workflows. It can be used as a starting point for larger sequence-analysis projects.
