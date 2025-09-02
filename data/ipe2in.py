import os
import xml.etree.ElementTree as ET

def extract_chains(filename):
    tree = ET.parse(filename)
    root = tree.getroot()

    # only consider paths inside <page>, skip <symbol>
    for page in root.iter('page'):
        for path in page.iter('path'):
            if path.text is None:
                continue

            tokens = path.text.strip().split()
            chains = []
            current = []
            start_pt = None
            numbuf = []  # holds numbers seen since last command
            i = 0
            while i < len(tokens):
                tok = tokens[i]

                if tok in ('m', 'l'):
                    # expect the last two numbers in number buffer to be the coordinate
                    if len(numbuf) >= 2:
                        try:
                            x = float(numbuf[-2])
                            y = float(numbuf[-1])
                        except ValueError:
                            x = y = None
                        # clear number buffer for next coord(s)
                        numbuf.clear()

                        if x is not None:
                            if tok == 'm':
                                # starting a new subpath; finalize any open one
                                if current:
                                    chains.append(current)
                                current = [(x, y)]
                                start_pt = (x, y)
                            else:  # 'l'
                                if not current:
                                    # robustness: lineto without moveto -> start here
                                    current = [(x, y)]
                                    start_pt = (x, y)
                                else:
                                    current.append((x, y))
                    else:
                        # command without enough numbers; ignore and reset buffer
                        numbuf.clear()
                        print (f"Warning: '{tok}' command without enough coordinates in {filename}")
                    i += 1
                    continue

                if tok == 'h':
                    # close current subpath: append start point if not already equal
                    if current:
                        if start_pt is not None and current[-1] != start_pt:
                            current.append(start_pt)
                        chains.append(current)
                        current = []
                        start_pt = None
                    numbuf.clear()
                    i += 1
                    continue

                # otherwise, try to parse as number token
                try:
                    float(tok)
                    numbuf.append(tok)
                except ValueError:
                    # unknown token/command – discard any pending numbers to avoid mispairing
                    numbuf.clear()
                i += 1

            if current:
                chains.append(current)

            # yield all collected chains for this <path>
            for ch in chains:
                if ch:
                    yield ch

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
