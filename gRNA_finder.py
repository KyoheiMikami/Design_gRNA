import argparse
import subprocess
import re
import tempfile
import os
from collections import defaultdict

# --- Constants ---
DNA_COMPLEMENT = str.maketrans("ATCG", "TAGC")
MIN_GRNA_LENGTH = 12
DEFAULT_GRNA_LENGTH = 20
DEFAULT_MISMATCHES = 3
DEFAULT_BULGE = 0
DEFAULT_OUTPUT_FILE = "output_GF.txt"
CASOFFINDER_DEFAULT_PATH = "cas-offinder"

# --- Utility Functions ---
def reverse_complement(seq: str) -> str:
    """Calculates the reverse complement of a DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]

def count_mismatches(seq: str) -> int:
    """Counts lowercase characters (mismatches) in a DNA sequence."""
    return sum(1 for nt in seq if nt.islower())

def split_regions(seq: str) -> tuple[str, str, str]:
    """Splits a sequence into distal, proximal, and PAM regions."""
    return seq[:-15], seq[-15:-3], seq[-3:]

# --- Core Functions ---
def find_gRNA_candidates(input_file: str, gRNA_length: int) -> list[str]:
    """Finds unique forward and reverse gRNA candidates ending with NGG."""
    if gRNA_length < MIN_GRNA_LENGTH:
        raise ValueError(f"Invalid gRNA length: must be ≥ {MIN_GRNA_LENGTH} nt")

    with open(input_file, "r") as f:
        dna_seq = f.read().upper()

    pattern = re.compile(f"(?=([ATGC]{{{gRNA_length}}}[ATGC]GG))")

    forward = pattern.findall(dna_seq)
    reverse = pattern.findall(reverse_complement(dna_seq))

    set_forward = set(forward)
    unique_reverse = [m for m in reverse if reverse_complement(m) not in set_forward]

    print(f"✓ Found {len(forward)} forward gRNA candidates.")
    print(f"✓ Found {len(reverse)} reverse gRNA candidates ({len(unique_reverse)} unique after filtering).")

    return forward + unique_reverse

def make_CO_input(gRNAs: list[str], genome_path: str, gRNA_length: int, mismatches: int, bulge: int, out_file: str):
    """Writes Cas-OFFinder input file."""
    with open(out_file, "w") as out:
        out.write(f"{genome_path}\n")
        out.write(f"{'N' * gRNA_length}NGG {bulge} {bulge}\n")
        for i, g in enumerate(gRNAs, 1):
            out.write(f"{g} {mismatches} g{i}\n")

def run_CO(input_file: str, output_file: str, path: str = CASOFFINDER_DEFAULT_PATH):
    """Runs Cas-OFFinder with error capture."""
    cmd = [path, input_file, "C", output_file]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✓ Cas-OFFinder completed successfully.")
    except subprocess.CalledProcessError as e:
        print("❌ Error running Cas-OFFinder:")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        raise

def passes_stringency(dna_seq: str, level: str) -> bool:
    """Returns True if off-target fails the stringency check (should be excluded)."""
    distal, proximal, _ = split_regions(dna_seq)
    dmm = count_mismatches(distal)
    pmm = count_mismatches(proximal)

    if level == "high":
        return (pmm == 0) or (pmm == 1 and dmm < 2)
    elif level == "maximum":
        return (pmm == 0) or (pmm == 1 and dmm < 5) or (pmm == 2 and dmm < 2)
    else:
        raise ValueError(f"Invalid stringency level: {level}")

def filter_offtargets(cas_output: str, output_file: str, stringency: str):
    """Filters Cas-OFFinder output and writes accepted gRNAs."""
    with open(cas_output, encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    headers = [l for l in lines if l.startswith("#")]
    entries = [l for l in lines if not l.startswith("#")]

    grouped = defaultdict(list)
    for line in entries:
        cols = line.split("\t")
        gRNA_id, aligned_seq = cols[0], cols[3]
        grouped[gRNA_id].append((line, aligned_seq))

    accepted_lines, accepted_ids = [], []

    for gRNA_id, records in grouped.items():
        has_on_target = False
        exclude = False

        for _, seq in records:
            mm = count_mismatches(seq)
            if mm == 0:
                has_on_target = True
                continue
            if passes_stringency(seq, stringency):
                exclude = True
                break

        if has_on_target and not exclude:
            accepted_ids.append(gRNA_id)
            for line, seq in records:
                mm = count_mismatches(seq)
                accepted_lines.append(line + (" # 0MM" if mm == 0 else ""))

    with open(output_file, "w", encoding="utf-8") as f:
        for h in headers:
            f.write(h + "\n")
        for l in accepted_lines:
            f.write(l + "\n")

    print(f"✓ Filtering complete. Valid gRNAs: {len(accepted_ids)}. Output: {output_file}")

# --- Main ---
def main():
    parser = argparse.ArgumentParser(description="gRNA candidate filter pipeline using Cas-OFFinder")
    parser.add_argument("-i", "--input_file", required=True, help="DNA sequence input (TXT)")
    parser.add_argument("-l", "--gRNA_length", type=int, default=DEFAULT_GRNA_LENGTH, help="gRNA length before PAM (default: 20)")
    parser.add_argument("-g", "--genome_path", required=True, help="Path to genome file (2bit or FASTA)")
    parser.add_argument("-m", "--mismatches", type=int, default=DEFAULT_MISMATCHES, help="Max allowed mismatches")
    parser.add_argument("-b", "--bulge", type=int, default=DEFAULT_BULGE, help="Max bulges allowed")
    parser.add_argument("--CO_path", default=CASOFFINDER_DEFAULT_PATH, help="Cas-OFFinder executable path")
    parser.add_argument("-s", "--stringency", choices=["high", "maximum"], default="high", help="Filtering stringency level")
    parser.add_argument("-t", "--save_temp", action="store_true", help="Save intermediate files")
    parser.add_argument("-o", "--output_file", default=DEFAULT_OUTPUT_FILE, help="Output file name")
    args = parser.parse_args()

    print("✓ Starting gRNA pipeline...\n")

    # Step 1: gRNA search
    gRNAs = find_gRNA_candidates(args.input_file, args.gRNA_length)
    if not gRNAs:
        print("❌ No gRNA candidates found. Exiting.")
        return

    # Step 2: Prepare temp or named files
    if args.save_temp:
        input_CO = "Cas-OFFinder_input.txt"
        output_CO = "Cas-OFFinder_output.txt"
        print(f"✓ Intermediate files will be saved: {input_CO}, {output_CO}")
    else:
        input_CO = tempfile.NamedTemporaryFile(delete=False).name
        output_CO = tempfile.NamedTemporaryFile(delete=False).name

    try:
        make_CO_input(gRNAs, args.genome_path, args.gRNA_length, args.mismatches, args.bulge, input_CO)
        run_CO(input_CO, output_CO, args.CO_path)
        filter_offtargets(output_CO, args.output_file, args.stringency)
    finally:
        if not args.save_temp:
            for f in [input_CO, output_CO]:
                if os.path.exists(f):
                    os.remove(f)
        else:
            print(f"✓ Intermediate files saved: {input_CO}, {output_CO}")

    print(f"✓ Pipeline completed. Output: {args.output_file}")

if __name__ == "__main__":
    main()
