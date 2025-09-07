# Saturation Mutagenesis Library Generator

## Overview
This Python script (`Saturation_mutation_NNK.py`) generates a saturation mutagenesis library with overlapping oligonucleotide (oligo) sequences for a given wild-type (WT) DNA sequence. It is designed to create oligos for each amino acid position (except the start codon) mutated to 19 alternative amino acids (excluding the WT amino acid) using *E. coli* optimized codons. The script ensures overlaps between fragments for downstream assembly (e.g., Gibson assembly) and includes safety zones to prevent mutations in overlap regions.

## Features
- Generates oligos with a target length (default: 300 bp).
- Maintains 60 bp overlaps (2 × 30 bp safe zones) between fragments to ensure assembly compatibility.
- Skips the WT amino acid at each position, generating 19 mutants per position.
- Includes upstream and downstream primer sequences for PCR amplification.
- Validates input DNA and primer sequences (only A, T, C, G allowed).
- Outputs a CSV file with oligo names, sequences, and lengths.
- Warns if any amino acid positions are not covered or if oligo lengths exceed the target.

## Prerequisites
- Python 3.6+
- Required libraries: `csv`, `math`, `argparse` (standard libraries)

## Installation
No additional installation is required beyond Python. Clone or download the script to your local machine.

```bash
git clone <repository_url>
cd <repository_directory>
```

## Usage
Run the script from the command line with the following arguments:

```bash
python Saturation_mutation_NNK_v1.py --wt_DNA <wild_type_dna_sequence> [options]
```

### Command-Line Arguments
| Argument         | Description                                                                 | Default Value                       |
|------------------|-----------------------------------------------------------------------------|-------------------------------------|
| `--wt_DNA`      | Wild-type DNA sequence (required, must contain only A, T, C, G, length divisible by 3). | None                                |
| `--up_primer`   | Upstream primer sequence.                                                   | `ACAATTCTGCCTAGGAGATCT`         |
| `--down_primer` | Downstream primer sequence.                                                 | `TGACATCTGtagtgcaACAAG`         |
| `--oligo_len`   | Target synthesis length for each oligo (bp).                                | 300                                 |
| `--safe_len`    | Safe zone length at each fragment end (bp, overlap = 2 × safe_len).        | 30                                  |
| `--output`      | Output CSV filename.                                                       | `saturation_library_with_overlaps.csv` |

### Example
```bash
python Saturation_mutation_NNK_v1.py --wt_DNA <your_wt_dna_sequence> --output s9_NNK_oligo.csv
```

### Output
- **CSV File**: Contains columns `Name`, `Sequence`, `Length`. Each row represents an oligo with a unique mutation (e.g., `Frag1_Pos2_A` for fragment 1, position 2, mutated to alanine).
- **Console Output**: Reports fragment details (slice, length, mutated positions), warnings for missing positions or length violations, and total oligos generated.

Example console output:
```
Fragment 1: slice [0:279], length 279, mutated aa 2 to 83
Fragment 2: slice [219:519], length 300, mutated aa 84 to 163
Fragment 3: slice [459:759], length 300, mutated aa 164 to 243
...
All amino acid positions (except start codon) covered.
Generated 8588 oligos. Saved to s9_NNK_oligo.csv
```

## Implementation Details
- **Fragmentation**: The WT DNA is divided into fragments with a step size adjusted for the first fragment (219 bp gene part + 21 bp primer = 240 bp) to ensure 60 bp overlaps (2 × 30 bp safe zones). Subsequent fragments use 300 bp with 60 bp overlaps.
- **Mutation Strategy**: For each mutable position, the script identifies the WT amino acid and generates mutants for the 19 other amino acids using *E. coli* optimized codons.
- **Safety Zones**: Mutations are skipped in the first 30 bp and last 30 bp of each fragment (except specific cases) to maintain WT sequence in overlaps for assembly.
- **Optimization for NNK**: Only 19 mutants per position are generated, excluding the WT amino acid, reducing oligo count compared to full 20-amino-acid saturation.

## Notes
- The first fragment’s total length may exceed 300 bp due to the upstream primer (e.g., 321 bp with default primer). Check synthesis platform limits.
- The script assumes the WT DNA does not include stop codons in the input sequence.
- Overlaps are designed for assembly methods like Gibson assembly, ensuring WT sequence in overlap regions.
- If the WT sequence is short, some positions may fall in safe zones and be reported as missing.

## Troubleshooting
- **Script hangs**: Likely due to large WT DNA or excessive oligo generation. Try reducing `oligo_len` or increasing `safe_len` to reduce fragment count.
- **Missing positions**: Increase `oligo_len` or reduce `safe_len` to cover more positions, but ensure sufficient overlap for assembly.
- **Length warnings**: Adjust `oligo_len` to accommodate primer lengths or check synthesis constraints.

## License
This project is licensed under the MIT License.

## Contact
For issues or suggestions, please open an issue on the repository or contact the maintainer.
