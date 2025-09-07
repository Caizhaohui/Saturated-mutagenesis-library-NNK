import csv
import math
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description="Generate saturation mutagenesis library oligos with overlaps.")
parser.add_argument('--wt_DNA', type=str, required=True, help="The wild-type DNA sequence for mutagenesis.")
parser.add_argument('--up_primer', type=str, default="ACAATTCTGCCTAGGAGATCT", help="Upstream primer site sequence.")
parser.add_argument('--down_primer', type=str, default="TGACATCTGtagtgcaACAAG", help="Downstream primer site sequence.")
parser.add_argument('--oligo_len', type=int, default=300, help="Target synthesis length for each oligo (default: 300 bp).")
parser.add_argument('--safe_len', type=int, default=30, help="Safe zone length at each end (default: 30 bp).")
parser.add_argument('--output', type=str, default="saturation_library_with_overlaps.csv", help="Output CSV filename (default: saturation_library_with_overlaps.csv).")

args = parser.parse_args()

# 从参数中获取值
wt_dna = args.wt_DNA
UPSTREAM_PRIMER_SITE = args.up_primer
DOWNSTREAM_PRIMER_SITE = args.down_primer
TARGET_SYNTHESIS_LEN = args.oligo_len
SAFE_ZONE_LEN = args.safe_len
output_filename = args.output

# 优选密码子（不包括终止子）
ECO_PREFERRED_CODONS = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}

# 核心区步长
step_size = TARGET_SYNTHESIS_LEN - (2 * SAFE_ZONE_LEN)

num_fragments = math.ceil(len(wt_dna) / step_size)

fragments = []
for i in range(num_fragments):
    frag_num = i + 1
    slice_start = step_size * i
    final_slice_start = slice_start
    final_slice_end = slice_start + TARGET_SYNTHESIS_LEN
    if frag_num == 1:
        final_slice_start = 0
        final_slice_end = TARGET_SYNTHESIS_LEN - len(UPSTREAM_PRIMER_SITE) if len(UPSTREAM_PRIMER_SITE) > 0 else TARGET_SYNTHESIS_LEN
    elif frag_num == num_fragments:
        remaining_gene = len(wt_dna) - slice_start
        final_oligo_len = min(TARGET_SYNTHESIS_LEN, remaining_gene + len(DOWNSTREAM_PRIMER_SITE))
        final_slice_end = slice_start + (final_oligo_len - len(DOWNSTREAM_PRIMER_SITE)) if len(DOWNSTREAM_PRIMER_SITE) > 0 else slice_start + remaining_gene
    fragments.append({
        'num': frag_num,
        'slice_start': final_slice_start,
        'slice_end': final_slice_end,
    })

last_frag_num = fragments[-1]['num']
variants = []
for frag in fragments:
    fragment_dna = wt_dna[frag['slice_start']:frag['slice_end']]
    frag_len = len(fragment_dna)
    first_codon_local_start = (3 - (frag['slice_start'] % 3)) % 3
    for codon_start_local in range(first_codon_local_start, frag_len, 3):
        if codon_start_local + 3 > frag_len: continue
        is_in_safe_zone = False
        if frag['num'] == 1 and (codon_start_local + 3) > (frag_len - SAFE_ZONE_LEN):
            is_in_safe_zone = True
        elif frag['num'] == last_frag_num and codon_start_local < SAFE_ZONE_LEN:
            is_in_safe_zone = True
        elif frag['num'] not in [1, last_frag_num] and (codon_start_local < SAFE_ZONE_LEN or (codon_start_local + 3) > (frag_len - SAFE_ZONE_LEN)):
            is_in_safe_zone = True
        if is_in_safe_zone: continue
        codon_start_global = frag['slice_start'] + codon_start_local
        global_aa_pos = (codon_start_global // 3) + 1
        if global_aa_pos == 1 and frag['num'] == 1: continue
        for aa_code, preferred_codon in ECO_PREFERRED_CODONS.items():
            mutant_dna = fragment_dna[:codon_start_local] + preferred_codon + fragment_dna[codon_start_local + 3:]
            description = f"Frag{frag['num']}_Pos{global_aa_pos}_{aa_code}"
            final_sequence = mutant_dna
            if frag['num'] == 1:
                final_sequence = UPSTREAM_PRIMER_SITE + mutant_dna
            elif frag['num'] == last_frag_num:
                final_sequence = mutant_dna + DOWNSTREAM_PRIMER_SITE
            variants.append((description, final_sequence, len(final_sequence)))

with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Name', 'Sequence', 'Length'])
    for name, seq, length in variants:
        writer.writerow([name, seq, length])

print(f"Generated {len(variants)} oligos. Saved to {output_filename}")
