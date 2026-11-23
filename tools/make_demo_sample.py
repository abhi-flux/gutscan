"""Create a SYNTHETIC demo stool FASTA with known composition (for testing the app).
Reads are 250-bp windows from reference sequences with sequencing errors, mixed strands,
plus some deliberately messy records (lower-case, N-runs, too short, protein, random DNA)."""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from validate import window, mutate
from gutscan.refdb import RefDB, revcomp, _read_fasta_ref, data_dir

random.seed(2024)
ids, lins, seqs = [], [], []
_read_fasta_ref(os.path.join(data_dir(), "reference.fa.gz"), ids, lins, seqs)
COMP = {  # genus (or species) -> number of reads
    ("g", "Bacteroides"): 2500, ("g", "Faecalibacterium"): 1200, ("g", "Prevotella"): 500, ("g", "Roseburia"): 400,
    ("g", "Akkermansia"): 300, ("g", "Bifidobacterium"): 300, ("g", "Lactobacillus"): 100,
    ("g", "Streptococcus"): 150, ("g", "Desulfovibrio"): 50, ("s", "gnavus"): 100, ("g", "Escherichia"): 20,
    ("g", "Klebsiella"): 120, ("f", "Lachnospiraceae"): 1500,
}
pool = {}
for i, l in enumerate(lins):
    pool.setdefault(("g", l[4]), []).append(i); pool.setdefault(("f", l[3]), []).append(i); pool.setdefault(("s", l[5]), []).append(i)
recs = []
for key, n in COMP.items():
    for _ in range(n):
        i = random.choice(pool[key]); recs.append(window(seqs[i], 250))
others = [i for i, l in enumerate(lins) if l[4] in ("Staphylococcus", "Pseudomonas", "Veillonella", "Blautia", "Alistipes", "Parabacteroides", "Sutterella", "Haemophilus", "Collinsella")]
for _ in range(2000):
    recs.append(window(seqs[random.choice(others)], 250))
random.shuffle(recs)
out = []
for k, r in enumerate(recs):
    out.append(f">read_{k}\n{r}")
# messy records
out.append(">lowercase_read\n" + window(seqs[pool[('g','Bacteroides')][0]], 250).lower())
out.append(">with_gaps_and_numbers\n" + " ".join(window(seqs[pool[('g','Faecalibacterium')][0]], 250)[i:i+10] for i in range(0, 250, 10)))
out.append(">too_short\nACGTACGTACGTACGTACGTAGCTAGCTAGCATCG")
out.append(">many_N\n" + "N" * 120 + window(seqs[pool[('g','Roseburia')][0]], 250)[:100])
out.append(">protein\nMKVLAAGIVGLLLAQPSTAWAEEKKRHHHSSSSDDDFFFGGGWWW" * 5)
out.append(">polyA\n" + "A" * 300)
for j in range(100):
    out.append(f">random_dna_{j}\n" + "".join(random.choice("ACGT") for _ in range(260)))
out.append(">abundant_asv;size=500;\n" + window(seqs[pool[('g','Dialister')][0]], 250))
open(os.path.join(os.path.dirname(__file__), "..", "sample_stool_DEMO_synthetic.fasta"), "w").write("\n".join(out) + "\n")
print("records:", len(out))
