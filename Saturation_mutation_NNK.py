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

# 从参数中获取值，并转为大写
wt_dna = args.wt_DNA.upper()
UPSTREAM_PRIMER_SITE = args.up_primer.upper()
DOWNSTREAM_PRIMER_SITE = args.down_primer.upper()
TARGET_SYNTHESIS_LEN = args.oligo_len
SAFE_ZONE_LEN = args.safe_len
output_filename = args.output

# 验证wt_DNA和引物序列
valid_bases = set('ATCG')
if not all(c in valid_bases for c in wt_dna):
    raise ValueError("wt_DNA contains invalid characters. Only A, T, C, G are allowed.")
if len(wt_dna) % 3 != 0:
    raise ValueError("wt_DNA length must be divisible by 3 to maintain codon frame.")
if not all(c in valid_bases for c in UPSTREAM_PRIMER_SITE):
    raise ValueError("up_primer contains invalid characters. Only A, T, C, G are allowed.")
if not all(c in valid_bases for c in DOWNSTREAM_PRIMER_SITE):
    raise ValueError("down_primer contains invalid characters. Only A, T, G, C are allowed.")

# 优选密码子（不包括终止子）
ECO_PREFERRED_CODONS = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}

# 遗传代码表，用于翻译WT密码子到氨基酸
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

fragments = []
overlap_len = 2 * SAFE_ZONE_LEN

# First fragment: gene part = oligo_len - up_primer, total <= oligo_len
gene_part_len = TARGET_SYNTHESIS_LEN - len(UPSTREAM_PRIMER_SITE)
slice_start = 0
slice_end = min(gene_part_len, len(wt_dna))
fragments.append({'num': 1, 'slice_start': slice_start, 'slice_end': slice_end})

# Subsequent fragments: slide to maintain overlap
slice_start = slice_end - overlap_len
frag_num = 2
mutated_positions = set()  # 跟踪已突变位点
all_positions = set(range(2, len(wt_dna) // 3 + 1))  # aa2到末尾

while True:
    if slice_start >= len(wt_dna):
        break
    is_last = (slice_start + TARGET_SYNTHESIS_LEN >= len(wt_dna))
    if is_last:
        remaining = len(wt_dna) - slice_start
        gene_part_len = min(TARGET_SYNTHESIS_LEN - len(DOWNSTREAM_PRIMER_SITE), remaining)
        if gene_part_len <= 0:
            break
    else:
        gene_part_len = TARGET_SYNTHESIS_LEN
    slice_end = slice_start + gene_part_len
    if slice_end > len(wt_dna):
        slice_end = len(wt_dna)
    fragments.append({'num': frag_num, 'slice_start': slice_start, 'slice_end': slice_end})
    if is_last:
        break
    slice_start = slice_end - overlap_len
    frag_num += 1

last_frag_num = fragments[-1]['num']
variants = []

for frag in fragments:
    fragment_dna = wt_dna[frag['slice_start']:frag['slice_end']]
    frag_len = len(fragment_dna)
    first_codon_local_start = (3 - (frag['slice_start'] % 3)) % 3
    mutable_positions = []

    for codon_start_local in range(first_codon_local_start, frag_len, 3):
        if codon_start_local + 3 > frag_len:
            continue

        # 统一安全区检查：避免重叠区域突变
        is_in_safe_zone = False
        # 所有片段的起始安全区（左侧重叠）
        if frag['num'] != 1 and codon_start_local < SAFE_ZONE_LEN:
            is_in_safe_zone = True
        # 所有片段的结束安全区（右侧重叠）
        if frag['num'] != last_frag_num and (codon_start_local + 3) > (frag_len - SAFE_ZONE_LEN):
            is_in_safe_zone = True
        # 第一片段的结束安全区
        if frag['num'] == 1 and (codon_start_local + 3) > (frag_len - SAFE_ZONE_LEN):
            is_in_safe_zone = True
        # 最后片段的起始安全区
        if frag['num'] == last_frag_num and codon_start_local < SAFE_ZONE_LEN:
            is_in_safe_zone = True

        if is_in_safe_zone:
            continue

        codon_start_global = frag['slice_start'] + codon_start_local
        global_aa_pos = (codon_start_global // 3) + 1
        if global_aa_pos == 1 and frag['num'] == 1:
            continue  # 跳过起始密码子
        if global_aa_pos not in mutated_positions:
            mutable_positions.append(global_aa_pos)
            wt_codon = fragment_dna[codon_start_local:codon_start_local + 3]
            wt_aa = CODON_TABLE.get(wt_codon, '*')  # 如果未知，视为终止
            for aa_code, preferred_codon in ECO_PREFERRED_CODONS.items():
                if aa_code == wt_aa:
                    continue  # 跳过WT氨基酸
                mutant_dna = fragment_dna[:codon_start_local] + preferred_codon + fragment_dna[codon_start_local + 3:]
                description = f"Frag{frag['num']}_Pos{global_aa_pos}_{aa_code}"
                final_sequence = mutant_dna
                if frag['num'] == 1:
                    final_sequence = UPSTREAM_PRIMER_SITE + mutant_dna
                if frag['num'] == last_frag_num:
                    final_sequence += DOWNSTREAM_PRIMER_SITE
                seq_len = len(final_sequence)
                if seq_len > TARGET_SYNTHESIS_LEN:
                    print(f"Warning: Oligo {description} length {seq_len} exceeds target {TARGET_SYNTHESIS_LEN} bp.")
                variants.append((description, final_sequence, seq_len))
            mutated_positions.add(global_aa_pos)

    print(f"Fragment {frag['num']}: slice [{frag['slice_start']}:{frag['slice_end']}], length {frag_len}, mutated aa {min(mutable_positions, default='none')} to {max(mutable_positions, default='none')}")

# 检查覆盖完整性
missing_positions = all_positions - mutated_positions
if missing_positions:
    print(f"Warning: Missing mutations for amino acid positions: {sorted(missing_positions)}")
else:
    print("All amino acid positions (except start codon) covered.")

with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Name', 'Sequence', 'Length'])
    for name, seq, length in variants:
        writer.writerow([name, seq, length])

print(f"Generated {len(variants)} oligos. Saved to {output_filename}")
