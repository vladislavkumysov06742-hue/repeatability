def read_reference(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            seq.append(line.strip().upper())
    return ''.join(seq)
