from itertools import product

ALPHABET = ["A", "C", "G", "T"]

def all_combinations(k: int):
    return ["".join(p) for p in product(ALPHABET, repeat=k)]

def counts(seq: str, k: int) -> dict:
    seq = seq.upper()
    combos = all_combinations(k)
    count = {c: 0 for c in combos}
    for i in range(len(seq) - k + 1):
        aux = seq[i:i+k]
        if aux in count:
            count[aux] += 1
    return count

def percentages(seq: str, k: int) -> dict:
    cnt = counts(seq, k)
    total = max(len(seq) - k + 1, 1)
    return {c: (v / total) * 100.0 for c, v in cnt.items()}

def main():
    S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
    print(f"Sequence: {S}")
    print(f"Length: {len(S)}\n")

    for k in [2, 3]:
        perc = percentages(S, k)
        print(f"--- k = {k} ---")
        print("-" * 30)
        for p in sorted(perc.keys()):
            print(f"{p}\t{perc[p]:6.2f}%")
        print()

if __name__ == "__main__":
    main()
