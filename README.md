# Saturation Mutagenesis Library Generator

This Python script generates a saturation mutagenesis library for a given wild-type DNA sequence. It creates oligo sequences for each fragment, ensuring that each amino acid position (except the start codon and stop codon) is mutated to all 20 standard amino acids. The script is designed to support fragment assembly with overlapping regions, suitable for PCR amplification and subsequent assembly (e.g., Gibson Assembly or Golden Gate with BsaI).

## Features
- **Input Flexibility**: Accepts a wild-type DNA sequence and customizable parameters via command-line arguments.
- **Fragmentation**: Automatically splits the input sequence into fragments of specified length (default 300 bp), with 30 bp safe zones at each end to serve as PCR primer binding sites or overlap regions for assembly.
- **Saturation Mutagenesis**: Generates oligos for each position (except start codon) to cover all 20 amino acids using optimized codons for *E. coli*.
- **Overlap Design**: Ensures 60 bp overlaps between fragments (30 bp from the end of one fragment + 30 bp from the start of the next) for reliable assembly.
- **Output**: Saves results to a CSV file with oligo names, sequences, and lengths.

## Prerequisites
- **Python**: Version 3.6 or higher.
- **Required Libraries**: 
  - `csv` (built-in)
  - `math` (built-in)
  - `argparse` (built-in)

No external dependencies like Biopython are required.

## Installation
1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-name>
   ```
2. Ensure Python 3 is installed:
   ```bash
   python3 --version
   ```

## Usage
Run the script from the command line with the required and optional parameters.

### Command
```bash
python3 script.py --wt_DNA <wild-type-DNA-sequence> [--up_primer <upstream-primer>] [--down_primer <downstream-primer>] [--oligo_len <length>] [--safe_len <length>] [--output <filename>]
```

### Parameters
- `--wt_DNA` (required): The wild-type DNA sequence (string, must be divisible by 3 to maintain codon frame).
- `--up_primer` (optional): Upstream primer sequence for the first fragment (default: `ACAATTCTGCCTAGGAGATCT`, 21 bp).
- `--down_primer` (optional): Downstream primer sequence for the last fragment (default: `TGACATCTGtagtgcaACAAG`, 21 bp).
- `--oligo_len` (optional): Target oligo synthesis length (default: 300 bp).
- `--safe_len` (optional): Length of the safe zone at each fragment's start/end, not mutated, used for PCR primers/overlaps (default: 30 bp).
- `--output` (optional): Output CSV filename (default: `saturation_library_with_overlaps.csv`).

### Example
```bash
python3 script.py --wt_DNA "ATGACTATCGCT" --up_primer "GGTCTCAACAATTCT" --down_primer "ACAAGTGAGACC" --oligo_len 300 --safe_len 30 --output my_library.csv
```
This generates a CSV file (`my_library.csv`) with oligos for the input sequence, adding the specified primers to the first and last fragments, with 30 bp safe zones and 300 bp target length.

## Output
The script generates a CSV file with the following columns:
- **Name**: Oligo description (e.g., `Frag1_Pos2_GCG` for fragment 1, position 2, mutated to alanine).
- **Sequence**: The oligo sequence, including primers for first/last fragments.
- **Length**: Length of the oligo in base pairs.

Example output in `my_library.csv`:
```csv
Name,Sequence,Length
Frag1_Pos2_GCG,ACAATTCTGCCTAGGAGATCTATGGCGATCGCT,33
Frag1_Pos2_CGC,ACAATTCTGCCTAGGAGATCTATGCGCATCGCT,33
...
```

## Design Details
- **Fragmentation**: The input sequence is split into fragments based on `step_size = oligo_len - 2*safe_len` (default: 300 - 2*30 = 240 bp). The number of fragments is calculated as `ceil(len(wt_dna) / step_size)`.
- **Overlaps**: Adjacent fragments overlap by 60 bp (30 bp from the end of one fragment + 30 bp from the start of the next), ensuring robust assembly. The first fragment has a shorter overlap (e.g., 39 bp) due to the upstream primer.
- **Safe Zones**: Each fragment has 30 bp at both ends (except first/last adjusted for primers) that are not mutated, serving as PCR primer binding sites or assembly overlaps.
- **Mutagenesis**: Excludes the start codon (position 1) and uses *E. coli*-optimized codons for 20 amino acids (no stop codon). Each mutable position generates 20 oligos.
- **Assembly**: Designed for PCR amplification followed by Gibson Assembly or Golden Gate (if BsaI sites are added to primers, e.g., `GGTCTCA` for upstream).

## Notes
- **Sequence Validation**: Ensure the input `wt_DNA` contains only valid DNA bases (A, T, C, G). Invalid characters (e.g., 'D', 'V') will cause synthesis errors.
- **Length Constraints**: If a fragment exceeds `oligo_len` after adding primers, the script adjusts the slice to fit (e.g., first fragment trims to 279 bp + 21 bp primer = 300 bp).
- **BsaI Support**: To use BsaI for Golden Gate assembly, prepend `GGTCTCA` to `--up_primer` and append `GAGACC` to `--down_primer` with appropriate overhangs.
- **NNK Alternative**: The current script generates one oligo per amino acid (20 per position). For a smaller library, modify to use NNK degenerate codons (1 oligo per position covering 20 AA + stop).
- **Testing**: Use a short wt_DNA sequence (e.g., 600 bp) to verify output before running on a long sequence (e.g., 1368 bp).

## Example Output for 1368 bp Sequence
For a 1368 bp sequence (456 amino acids), the script typically generates:
- **Fragments**: ~6 (depending on step_size).
- **Oligos**: ~8960 (448 mutable positions × 20 amino acids).
- **Overlap**: 60 bp between fragments, except ~39 bp for first fragment due to primer.
- **File**: CSV with ~8960 rows, each with a unique oligo sequence.

## Troubleshooting
- **Invalid Sequence**: If the wt_DNA contains non-DNA characters, replace them (e.g., 'D' → 'C', 'V' → 'G') before running.
- **Too Many Oligos**: Switch to NNK degeneracy to reduce the library size (~450 oligos).
- **Short Last Fragment**: If the last fragment is too short (<150 bp), check synthesis provider requirements.
- **Assembly Issues**: Ensure primers include enzyme sites (e.g., BsaI) if needed, and verify overlap lengths (60 bp is robust).

## License
MIT License. See `LICENSE` file for details.

## Contact
For issues or feature requests, open a GitHub issue or contact the repository maintainer.
