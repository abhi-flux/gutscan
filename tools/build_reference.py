"""Build gutscan/data/reference.fa.gz from the Greengenes 13_8 97% OTU set.

Source: PyPI package `qiime-default-reference` (Greengenes 13_8, CC BY-SA 3.0).
Usage:  python tools/build_reference.py
Selection: every sequence of the target genera/families (capped per genus so the
database stays small), plus up to N decoy sequences per non-target genus so that
unrelated bacteria are recognised as "other" instead of being forced onto a target.
"""
import gzip, os, random, re, subprocess, sys, tempfile, tarfile, glob

OUT = os.path.join(os.path.dirname(__file__), "..", "gutscan", "data", "reference.fa.gz")
TARGET_GENERA = {"Faecalibacterium", "Roseburia", "Eubacterium", "Akkermansia", "Bacteroides",
                 "Lactobacillus", "Bifidobacterium", "Enterococcus", "Streptococcus", "Escherichia",
                 "Bacillus", "Clostridium", "Coprococcus", "Dialister", "Prevotella", "Desulfovibrio",
                 "Eggerthella", "Ruminococcus", "Klebsiella", "Shigella", "Salmonella"}
TARGET_FAMILIES = {"Lachnospiraceae", "Enterobacteriaceae"}
CAP_TARGET_GENUS = 150     # per genus inside target set
CAP_TARGET_SPECIES = 60    # per labelled species (keeps rare species, e.g. C. difficile, fully)
CAP_DECOY_GENUS = 8
CAP_DECOY_FAMILY_UNLABELLED = 25
MIN_LEN = 900
random.seed(7)


def clean(x):
    return re.sub(r"^[kpcofgs]__", "", x.strip()).strip("[]")


def find_pkg():
    d = tempfile.mkdtemp()
    subprocess.check_call([sys.executable, "-m", "pip", "download", "qiime-default-reference",
                           "--no-deps", "--no-binary", ":all:", "-d", d, "-q"])
    tgz = glob.glob(os.path.join(d, "*.tar.gz"))[0]
    with tarfile.open(tgz) as t:
        t.extractall(d)
    return glob.glob(os.path.join(d, "*", "qiime_default_reference", "gg_13_8_otus"))[0]


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else find_pkg()
    tax = {}
    with gzip.open(os.path.join(base, "taxonomy", "97_otu_taxonomy.txt.gz"), "rt") as f:
        for line in f:
            i, t = line.rstrip("\n").split("\t")
            parts = [clean(p) for p in t.split(";")]
            parts += [""] * (7 - len(parts))
            tax[i] = parts[1:7]          # phylum..species
    seqs = {}
    with gzip.open(os.path.join(base, "rep_set", "97_otus.fasta.gz"), "rt") as f:
        cur = None; buf = []
        for line in f:
            if line[0] == ">":
                if cur: seqs[cur] = "".join(buf)
                cur = line[1:].strip(); buf = []
            else:
                buf.append(line.strip())
        if cur: seqs[cur] = "".join(buf)

    groups = {}
    for i, lin in tax.items():
        s = seqs.get(i)
        if not s or len(s) < MIN_LEN or "N" in s[:50]:
            continue
        fam, gen, sp = lin[3], lin[4], lin[5]
        is_t = gen in TARGET_GENERA or fam in TARGET_FAMILIES
        if is_t:
            key = ("T", gen or fam, sp if sp else "")
        elif gen:
            key = ("D", gen, "")
        else:
            key = ("DF", fam or lin[2] or lin[0], "")
        groups.setdefault(key, []).append(i)

    keep = []
    # target: cap per genus, but keep every species-labelled group separately
    tgen = {}
    for (kind, name, sp), ids in groups.items():
        if kind == "T":
            tgen.setdefault(name, {}).setdefault(sp, []).extend(ids)
    for name, spd in tgen.items():
        chosen = []
        for sp, ids in spd.items():
            random.shuffle(ids)
            chosen += ids[:CAP_TARGET_SPECIES] if sp else []
        rest = [i for sp, ids in spd.items() if not sp for i in ids]
        random.shuffle(rest)
        chosen += rest[:max(0, CAP_TARGET_GENUS - len(chosen))]
        keep += chosen
    for (kind, name, sp), ids in groups.items():
        if kind == "D":
            random.shuffle(ids); keep += ids[:CAP_DECOY_GENUS]
        elif kind == "DF":
            random.shuffle(ids); keep += ids[:CAP_DECOY_FAMILY_UNLABELLED]
    keep = sorted(set(keep), key=lambda x: int(x))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with gzip.open(OUT, "wt", compresslevel=9) as o:
        for i in keep:
            o.write(">%s|%s\n%s\n" % (i, "|".join(tax[i]), seqs[i]))
    print("wrote", len(keep), "sequences ->", os.path.abspath(OUT), os.path.getsize(OUT) // 1024, "KB")


if __name__ == "__main__":
    main()
