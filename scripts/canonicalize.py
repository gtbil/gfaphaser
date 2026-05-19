import sys
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
for line in sys.stdin:
    parts = line.rstrip("\n").rstrip("\t").split("\t")
    if len(parts) < 2:
        continue
    k = parts[0].upper()
    n = parts[1].split(":", 1)[0]
    rc = k.translate(COMP)[::-1]
    ck = k if k <= rc else rc
    sys.stdout.write(f"{ck}\t{n}\n")
