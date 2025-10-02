seq = "ATTGCCCCGAAT"
alph = set(seq)

for base in alph:
    freq = seq.count(base) / len(seq)
    print(base, freq)
