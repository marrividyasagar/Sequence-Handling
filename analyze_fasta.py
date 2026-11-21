
import argparse
import csv
from pathlib import Path
from collections import Counter

# Valid DNA/RNA nucleotide symbols, including ambiguity codes
VALID_SYMBOLS = set("ACGTUNRYKMSWBDHV-")


def parse_fasta(file_path):
    """Read a FASTA file and return sequence IDs and sequences."""
    sequences = []
    current_header = None
    current_sequence = []

    with open(file_path, "r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequences.append(
                        (current_header, "".join(current_sequence))
                    )

                current_header = line[1:].strip()

                if not current_header:
                    raise ValueError(
                        f"Empty FASTA header at line {line_number}"
                    )

                current_sequence = []

            else:
                if current_header is None:
                    raise ValueError(
                        f"Sequence found before FASTA header "
                        f"at line {line_number}"
                    )

                current_sequence.append(
                    "".join(line.split()).upper()
                )

    if current_header is not None:
        sequences.append(
            (current_header, "".join(current_sequence))
        )

    if not sequences:
        raise ValueError("No sequences found in the FASTA file.")

    return sequences


def analyze_sequence(sequence):
    """Calculate length, composition and GC content."""

    counts = Counter(sequence)
    length = len(sequence)

    canonical_bases = (
        counts["A"] + counts["C"] + counts["G"]
        + counts["T"] + counts["U"]
    )

    gc_count = counts["G"] + counts["C"]

    gc_content = (
        (gc_count / canonical_bases) * 100
        if canonical_bases else 0.0
    )

    composition = {}

    for base in "ACGTUNRYKMSWBDHV-":
        composition[base] = counts.get(base, 0)

    return {
        "length": length,
        "gc_content_percent": round(gc_content, 2),
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"],
        "U": counts["U"],
        "ambiguous_bases": sum(
            counts[b] for b in "NRYKMSWBDHV"
        ),
        "gaps": counts["-"],
        "composition": composition
    }


def main():
    parser = argparse.ArgumentParser(
        description="Analyze nucleotide sequences in FASTA format."
    )

    parser.add_argument(
        "input",
        help="Path to the input FASTA file"
    )

    parser.add_argument(
        "--outdir",
        default="results",
        help="Folder for output files"
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.outdir)

    if not input_file.is_file():
        raise FileNotFoundError(
            f"Input FASTA file not found: {input_file}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    sequences = parse_fasta(input_file)

    rows = []
    total_counts = Counter()
    total_length = 0
    total_canonical = 0
    total_gc = 0

    for header, sequence in sequences:

        invalid = set(sequence) - VALID_SYMBOLS

        if invalid:
            raise ValueError(
                f"Invalid nucleotide symbols in {header}: "
                f"{', '.join(sorted(invalid))}"
            )

        if not sequence:
            raise ValueError(
                f"Empty sequence found: {header}"
            )

        result = analyze_sequence(sequence)

        row = {
            "sequence_id": header,
            "length": result["length"],
            "GC_content_percent": result["gc_content_percent"],
            "A": result["A"],
            "T": result["T"],
            "G": result["G"],
            "C": result["C"],
            "U": result["U"],
            "ambiguous_bases": result["ambiguous_bases"],
            "gaps": result["gaps"]
        }

        rows.append(row)

        total_length += result["length"]
        total_counts.update(result["composition"])
        total_gc += result["G"] + result["C"]
        total_canonical += (
            result["A"] + result["T"]
            + result["G"] + result["C"]
            + result["U"]
        )

    # Save per-sequence results
    csv_path = output_dir / "sequence_statistics.csv"

    with open(csv_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=rows[0].keys()
        )
        writer.writeheader()
        writer.writerows(rows)

    # Calculate combined GC content
    overall_gc = (
        (total_gc / total_canonical) * 100
        if total_canonical else 0.0
    )

    # Save overall summary
    summary_path = output_dir / "summary_report.txt"

    with open(summary_path, "w", encoding="utf-8") as file:
        file.write("FASTA SEQUENCE ANALYSIS REPORT\n")
        file.write("=" * 35 + "\n\n")

        file.write(f"Input file: {input_file.name}\n")
        file.write(f"Number of sequences: {len(sequences)}\n")
        file.write(f"Total sequence length: {total_length}\n")
        file.write(f"Overall GC content: {overall_gc:.2f}%\n\n")

        file.write("Combined nucleotide composition:\n")

        for base in "ACGTUNRYKMSWBDHV-":
            count = total_counts[base]
            percentage = (
                count / total_length * 100
                if total_length else 0
            )
            file.write(
                f"{base}: {count} ({percentage:.2f}%)\n"
            )

    print("FASTA analysis completed successfully!")
    print(f"Sequences analyzed: {len(sequences)}")
    print(f"Total length: {total_length}")
    print(f"Overall GC content: {overall_gc:.2f}%")
    print(f"CSV output: {csv_path}")
    print(f"Summary output: {summary_path}")


if __name__ == "__main__":
    main()
