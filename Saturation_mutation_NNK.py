import csv
import math
import argparse

def main():
    # --- 命令行参数解析 ---
    parser = argparse.ArgumentParser(description="Generate saturation mutagenesis library oligos with overlaps.")
    parser.add_argument('--wt_DNA', type=str, required=True, help="The wild-type DNA sequence for mutagenesis.")
    parser.add_argument('--up_primer', type=str, default="ACAATTCTGCCTAGGAGATCT", help="Upstream primer site sequence.")
    parser.add_argument('--down_primer', type=str, default="TGACATCTGtagtgcaACAAG", help="Downstream primer site sequence.")
    parser.add_argument('--oligo_len', type=int, default=300, help="Target synthesis length for each oligo (default: 300 bp).")
    parser.add_argument('--safe_len', type=int, default=30, help="Safe zone length at each end for overlaps (default: 30 bp).")
    parser.add_argument('--restriction_site', type=str, default="", help="Restriction site sequence to avoid (e.g., GAATTC for EcoRI).")
    parser.add_argument('--output', type=str, default="saturation_library_with_overlaps.csv", help="Output CSV filename.")
    args = parser.parse_args()

    # --- 从参数中获取值并进行初始化 ---
    wt_dna = args.wt_DNA.upper()
    UPSTREAM_PRIMER_SITE = args.up_primer.upper()
    DOWNSTREAM_PRIMER_SITE = args.down_primer.upper()
    RESTRICTION_SITE = args.restriction_site.upper() if args.restriction_site else None
    TARGET_SYNTHESIS_LEN = args.oligo_len
    SAFE_ZONE_LEN = args.safe_len
    output_filename = args.output

    # --- 输入序列合法性验证 ---
    valid_bases = set('ATCG')
    if not all(c in valid_bases for c in wt_dna):
        raise ValueError("wt_DNA contains invalid characters. Only A, T, C, G are allowed.")
    if len(wt_dna) % 3 != 0:
        raise ValueError("wt_DNA length must be divisible by 3 to maintain codon frame.")
    if not all(c in valid_bases for c in UPSTREAM_PRIMER_SITE):
        raise ValueError("up_primer contains invalid characters. Only A, T, C, G are allowed.")
    if not all(c in valid_bases for c in DOWNSTREAM_PRIMER_SITE):
        raise ValueError("down_primer contains invalid characters. Only A, T, G, C are allowed.")
    if RESTRICTION_SITE and not all(c in valid_bases for c in RESTRICTION_SITE):
        raise ValueError("restriction_site contains invalid characters. Only A, T, C, G are allowed.")

    # --- 生物学数据表 ---
    ECO_PREFERRED_CODONS = {
        'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC', 'Q': 'CAG', 'E': 'GAA',
        'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
        'P': 'CCG', 'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    }
    SYNONYM_CODONS = {
        'A': ['GCG', 'GCT', 'GCC', 'GCA'], 'R': ['CGC', 'CGT', 'CGA', 'CGG', 'AGA', 'AGG'],
        'N': ['AAC', 'AAT'], 'D': ['GAT', 'GAC'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA'],
        'E': ['GAA', 'GAG'], 'G': ['GGC', 'GGT', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'], 'L': ['CTG', 'CTT', 'CTC', 'CTA', 'TTA', 'TTG'],
        'K': ['AAA', 'AAG'], 'M': ['ATG'], 'F': ['TTT', 'TTC'], 'P': ['CCG', 'CCT', 'CCC', 'CCA'],
        'S': ['AGC', 'TCT', 'TCC', 'TCA', 'AGT', 'TCG'], 'T': ['ACC', 'ACT', 'ACA', 'ACG'],
        'W': ['TGG'], 'Y': ['TAT', 'TAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    }
    CODON_TABLE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    # --- 基因片段划分逻辑 ---
    fragments = []
    overlap_len = 2 * SAFE_ZONE_LEN

    # 第一个片段
    gene_part_len_frag1 = TARGET_SYNTHESIS_LEN - len(UPSTREAM_PRIMER_SITE)
    slice_start = 0
    slice_end = min(gene_part_len_frag1, len(wt_dna))
    fragments.append({'num': 1, 'slice_start': slice_start, 'slice_end': slice_end})

    # 后续片段
    slice_start = slice_end - overlap_len
    frag_num = 2
    while slice_start < len(wt_dna):
        is_last = (slice_start + TARGET_SYNTHESIS_LEN - len(DOWNSTREAM_PRIMER_SITE) >= len(wt_dna))
        
        if is_last:
            gene_part_len = len(wt_dna) - slice_start
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
    mutated_positions = set()
    all_positions = set(range(2, len(wt_dna) // 3 + 1))

    # --- 突变生成主循环 ---
    for frag in fragments:
        fragment_dna = wt_dna[frag['slice_start']:frag['slice_end']]
        frag_len = len(fragment_dna)
        first_codon_local_start = (3 - (frag['slice_start'] % 3)) % 3
        mutable_positions_in_frag = []

        for codon_start_local in range(first_codon_local_start, frag_len, 3):
            if codon_start_local + 3 > frag_len:
                continue

            # 安全区检查，避免在重叠区域引入突变
            is_in_safe_zone = (
                (frag['num'] != 1 and codon_start_local < SAFE_ZONE_LEN) or
                (frag['num'] != last_frag_num and (codon_start_local + 3) > (frag_len - SAFE_ZONE_LEN))
            )
            if is_in_safe_zone:
                continue

            codon_start_global = frag['slice_start'] + codon_start_local
            global_aa_pos = (codon_start_global // 3) + 1

            if global_aa_pos == 1:
                continue
            
            if global_aa_pos not in mutated_positions:
                mutable_positions_in_frag.append(global_aa_pos)
                wt_codon = fragment_dna[codon_start_local:codon_start_local + 3]
                wt_aa = CODON_TABLE.get(wt_codon)
                if not wt_aa:
                    raise ValueError(f"Invalid codon '{wt_codon}' at global DNA position {codon_start_global}.")

                # 遍历所有目标氨基酸
                for aa_code in ECO_PREFERRED_CODONS.keys():
                    if aa_code == wt_aa:
                        continue  # 跳过同义突变

                    best_codon_found = False
                    # 遍历目标氨基酸的所有同义密码子，寻找一个不产生酶切位点的
                    for codon_candidate in SYNONYM_CODONS.get(aa_code, []):
                        mutant_dna = fragment_dna[:codon_start_local] + codon_candidate + fragment_dna[codon_start_local + 3:]
                        
                        final_sequence = mutant_dna
                        if frag['num'] == 1:
                            final_sequence = UPSTREAM_PRIMER_SITE + mutant_dna
                        if frag['num'] == last_frag_num:
                            final_sequence = final_sequence[:len(mutant_dna)] + DOWNSTREAM_PRIMER_SITE
                        
                        # 如果提供了限制性酶切位点，则进行检查
                        if RESTRICTION_SITE and RESTRICTION_SITE in final_sequence:
                            continue  # 存在酶切位点，尝试下一个同义密码子
                        
                        # 找到了一个好的序列，记录下来并跳出
                        description = f"Frag{frag['num']}_Pos{global_aa_pos}_{aa_code}"
                        if codon_candidate != ECO_PREFERRED_CODONS[aa_code]:
                            description += "_AltCodon"
                        
                        seq_len = len(final_sequence)
                        if seq_len > TARGET_SYNTHESIS_LEN:
                             print(f"Warning: Oligo {description} length {seq_len} exceeds target {TARGET_SYNTHESIS_LEN} bp. Trimming.")
                             final_sequence = final_sequence[:TARGET_SYNTHESIS_LEN]
                             seq_len = len(final_sequence)

                        variants.append((description, final_sequence, seq_len))
                        best_codon_found = True
                        break
                    
                    # 如果所有同义密码子都无法避免酶切位点，则记录警告并使用最优密码子
                    if not best_codon_found:
                        print(f"Warning: Could not avoid restriction site for mutation Pos{global_aa_pos}_{aa_code}. Using preferred codon despite the site.")
                        preferred_codon = ECO_PREFERRED_CODONS[aa_code]
                        mutant_dna = fragment_dna[:codon_start_local] + preferred_codon + fragment_dna[codon_start_local + 3:]
                        
                        final_sequence = mutant_dna
                        if frag['num'] == 1:
                            final_sequence = UPSTREAM_PRIMER_SITE + mutant_dna
                        if frag['num'] == last_frag_num:
                            final_sequence += DOWNSTREAM_PRIMER_SITE

                        description = f"Frag{frag['num']}_Pos{global_aa_pos}_{aa_code}_RS_Unavoided"
                        seq_len = len(final_sequence)
                        if seq_len > TARGET_SYNTHESIS_LEN:
                             print(f"Warning: Oligo {description} length {seq_len} exceeds target {TARGET_SYNTHESIS_LEN} bp. Trimming.")
                             final_sequence = final_sequence[:TARGET_SYNTHESIS_LEN]
                             seq_len = len(final_sequence)
                        
                        variants.append((description, final_sequence, seq_len))

                mutated_positions.add(global_aa_pos)
        
        mutable_range = f"{min(mutable_positions_in_frag)} to {max(mutable_positions_in_frag)}" if mutable_positions_in_frag else "none"
        print(f"Fragment {frag['num']}: slice [{frag['slice_start']}:{frag['slice_end']}], length {frag_len}, mutated aa {mutable_range}")

    # --- 检查覆盖完整性并输出结果 ---
    missing_positions = all_positions - mutated_positions
    if missing_positions:
        print(f"\nWarning: Missing mutations for amino acid positions: {sorted(list(missing_positions))}")
    else:
        print("\nSuccess: All amino acid positions (from 2 to the end) are covered by the design.")

    try:
        with open(output_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Name', 'Sequence', 'Length'])
            writer.writerows(variants)
        print(f"Generated {len(variants)} oligos. Saved to {output_filename}")
    except IOError as e:
        raise IOError(f"Failed to write to {output_filename}: {e}")

if __name__ == "__main__":
    main()
