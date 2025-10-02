from collections import Counter
from pathlib import Path

FILE = "fasta.txt"  # FASTA file

def read_fasta_sequence(path: str) -> str:
    """Read a FASTA file and return the concatenated sequence (uppercase, no spaces)."""
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):  # skip headers & empty lines
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()

def main():
    if not Path(FILE).exists():
        print(f"File not found: {FILE}")
        return

    seq = read_fasta_sequence(FILE)
    if not seq:
        print("No sequence found (empty file or only headers).")
        return

    counts = Counter(seq)
    total = sum(counts.values())

    # Print alphabet (unique symbols)
    canonical = ["A", "C", "G", "T"]
    others = sorted([b for b in counts.keys() if b not in canonical])
    ordered = [b for b in canonical if b in counts] + others

    print("Alphabet (unique symbols):")
    print(", ".join(ordered))
    print()

    print("Relative frequencies (symbol : count / total):")
    for b in ordered:
        print(f"{b}: {counts[b]} / {total} = {counts[b]/total:.6f}")
    print()

    # GC content (only over A/C/G/T bases)
    acgt_total = sum(counts.get(b, 0) for b in canonical)
    gc = counts.get("G", 0) + counts.get("C", 0)
    if acgt_total > 0:
        gc_content = gc / acgt_total
        print(f"GC content (over A/C/G/T only): {gc} / {acgt_total} = {gc_content:.6%}")
    else:
        print("GC content: N/A (no A/C/G/T found)")

if __name__ == "__main__":
    main()
