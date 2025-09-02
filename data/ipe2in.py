import os
import xml.etree.ElementTree as ET

def extract_chains(filename):
    tree = ET.parse(filename)
    root = tree.getroot()

    # only consider paths inside <page>, skip <symbol>
    for page in root.iter('page'):
        for path in page.iter('path'):
            coords = []
            tokens = path.text.strip().split()
            i = 0
            while i < len(tokens):
                tok = tokens[i]

                # commands we know
                if tok in ('m', 'l', 'h'):
                    i += 1
                    continue

                # try parsing coordinates
                try:
                    x = float(tok)
                    y = float(tokens[i+1])
                    coords.append((x, y))
                    i += 2
                except (ValueError, IndexError):
                    i += 1

            if coords:
                yield coords

def write_chains(chains, outname):
    with open(outname, "w") as f:
        for chain in chains:
            for (x, y) in chain:
                f.write(f"{x:.6f} {y:.6f}\n")
            f.write("\n")

def main():
    folder = os.path.join(os.getcwd(), "data/IPE")
    for fname in os.listdir(folder):
        if fname.lower().endswith(".ipe"):
            inpath = os.path.join(folder, fname)

            subdir = os.path.join(folder, os.path.splitext(fname)[0])
            os.makedirs(subdir, exist_ok=True)

            outpath = os.path.join(subdir, "data.in")

            try:
                chains = list(extract_chains(inpath))
            except ET.ParseError as e:
                print(f"Skipping {fname}: parse error ({e})")
                continue

            if chains:
                write_chains(chains, outpath)
                print(f"Converted {fname} -> {os.path.relpath(outpath, folder)}")
            else:
                print(f"No paths found in {fname}")

if __name__ == "__main__":
    main()
