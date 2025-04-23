import argparse
import subprocess
import re
import tempfile
import os
from collections import defaultdict

DNA_COMPLEMENT = str.maketrans("ATCG", "TAGC")

def reverse_complement(seq):
    return seq.translate(DNA_COMPLEMENT)[::-1]

def find_gRNA_candidates(input_file, gRNA_length):
    if gRNA_length < 12:
        raise ValueError(f"invalid gRNA length: set > 12 nt value")
        
    with open(input_file, "r") as f:
        dna_seq = f.read().upper()

    pattern = f"(?=([ATGC]{{{gRNA_length}}}[ATGC]GG))"

    forward_matches = re.findall(pattern, dna_seq)
    reverse_matches = re.findall(pattern, reverse_complement(dna_seq))

    set_f = set(forward_matches)
    set_r = set([reverse_complement(match) for match in reverse_matches])

    unique_reverse_matches = [match for match in reverse_matches if reverse_complement(match) not in set_f]

    return forward_matches + unique_reverse_matches

def make_CO_input(gRNA_candidates, genome_path, gRNA_length, mismatches, bulge, input_file):    
    with open(input_file, "w") as out:
        out.write(genome_path + "\n")
        out.write(f"{'N'*gRNA_length}NGG {bulge} {bulge}\n")
        for i, g in enumerate(gRNA_candidates, 1):
            out.write(f"{g} {mismatches} g{i}\n")

def run_CO(input_CO, output_CO, CO_path="cas-offinder"):
    command = [CO_path, input_CO, "C", output_CO]
    subprocess.run(command, check=True)
    print(f"✓ Cas-OFFinder finished.")

def split_regions(dna_seq):
    pam = dna_seq[-3:]
    proximal = dna_seq[-15:-3]
    distal = dna_seq[:-15]
    return distal, proximal, pam

def count_mismatches(dna_seq):
    return sum(1 for nt in dna_seq if nt.islower())

def check_condition(dna_seq, stringency="high", is_exact_match=False):
    if is_exact_match:
        return True
    distal, proximal, _ = split_regions(dna_seq)
    distal_mm = count_mismatches(distal)
    proximal_mm = count_mismatches(proximal)

    if stringency == "high":
        return not ((proximal_mm == 0) or (proximal_mm == 1 and distal_mm < 2))
    elif stringency == "maximum":
        return not ((proximal_mm == 0) or 
                    (proximal_mm == 1 and distal_mm < 5) or
                    (proximal_mm == 2 and distal_mm < 2)
                   )
    else:
        raise ValueError(f"Invalid stringency level: {stringency}")

def filter_offtargets(output_CO, output_file, stringency="high"):
    with open(output_CO, encoding="utf-8") as f:
        lines = f.readlines()

    header = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#") and line.strip()]

    # gRNA IDごとにグループ化
    grouped = defaultdict(list)
    for line in data:
        cols = line.strip().split("\t")
        id_ = cols[0]
        dna_seq = cols[3]
        grouped[id_].append((line.strip(), dna_seq))

    accepted_lines = []
    accepted_ids = []

    for id_, entries in grouped.items():
        exclude_this_gRNA = False

        for line, dna in entries:
            mm_count = count_mismatches(dna)
            if mm_count == 0:
                continue  # ミスマッチ0（オンターゲット候補）は常に許容
            distal, proximal, _ = split_regions(dna)
            distal_mm = count_mismatches(distal)
            proximal_mm = count_mismatches(proximal)

            if stringency == "high":
                if (proximal_mm == 0) or (proximal_mm == 1 and distal_mm < 2):
                    exclude_this_gRNA = True
                    break
            elif stringency == "maximum":
                if (proximal_mm == 0) or \
                   (proximal_mm == 1 and distal_mm < 5) or \
                   (proximal_mm == 2 and distal_mm < 2):
                    exclude_this_gRNA = True
                    break
            else:
                raise ValueError(f"Invalid stringency level: {stringency}")

        if not exclude_this_gRNA:
            for line, dna in entries:
                mm = count_mismatches(dna)
                flag = " # 0MM" if mm == 0 else ""
                accepted_lines.append(line + flag)
            accepted_ids.append(id_)

    with open(output_file, "w", encoding="utf-8") as f:
        for line in header + accepted_lines:
            f.write(line + "\n")

    print(f"\n✓ Filtering complete. Valid gRNAs: {len(accepted_ids)}. Output: {output_file}\n")


def main():
    parser = argparse.ArgumentParser(description="Full gRNA candidate pipeline with Cas-OFFinder")
    parser.add_argument("-i", "--input_file", required=True, help="Input DNA sequence (as .txt)")
    parser.add_argument("-l", "--gRNA_length", type=int, default=20, help="Number of bases before PAM (default: 20)")
    parser.add_argument("-g", "--genome_path", required=True, help="Path to genome file (e.g. 2bit)")
    parser.add_argument("-m", "--mismatches", type=int, default=3, help="Max number of mismatches")
    parser.add_argument("-b", "--bulge", type=int, default=0, help="Max number of bulges")
    parser.add_argument("--CO_path", default="cas-offinder", help="Cas-OFFinder executable path")
    parser.add_argument("-s", "--stringency", choices=["high", "maximum"], default="high", help="Stringency level for off-target filtering (default: high)")
    parser.add_argument("-t", "--save_temp", action="store_true", help="Save intermediate input/output files")
    parser.add_argument("-o", "--output_file", default="output_GF.txt", help="Filtered output file")
    args = parser.parse_args()

    gRNA_candidates = find_gRNA_candidates(args.input_file, args.gRNA_length)
    print(f"\n✓ Found {len(gRNA_candidates)} gRNA candidates (forward + reverse).\n")

    # 一時ファイルまたは保存
    if args.save_temp:
        input_CO = "Cas-OFFinder_input.txt"
        output_CO = "Cas-OFFinder_output.txt"
    else:
        input_CO_fd, input_CO = tempfile.mkstemp()
        os.close(input_CO_fd)
        output_CO_fd, output_CO = tempfile.mkstemp()
        os.close(output_CO_fd)

    make_CO_input(gRNA_candidates, args.genome_path, args.gRNA_length, args.mismatches, args.bulge, input_CO)
    run_CO(input_CO, output_CO, args.CO_path)
    filter_offtargets(output_CO, args.output_file, stringency=args.stringency)

    if not args.save_temp:
        os.remove(input_CO)
        os.remove(output_CO)
    else:
        print(f"✓ Intermediate files saved as: {input_CO}, {output_CO}\nß")
    
    print(f"✓ Final filtered output: {args.output_file}\n")
        
if __name__ == "__main__":
    main()