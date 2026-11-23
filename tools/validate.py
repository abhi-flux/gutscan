"""Benchmark GutScan on simulated reads.

Mode A  (held-out, strict): 15% of reference sequences are removed from the database and
        reads simulated from them. Each held-out sequence is a *different 97% OTU* from
        everything left in the DB, so this is a harsh "novel strain" test.
Mode B  (seen strains): reads simulated from sequences that ARE in the database.
Both add sequencing errors (substitutions + indels) and random strand.
Truth = Greengenes labels of the source sequence, matched with the same target rules.
"""
import os, random, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from collections import defaultdict
from gutscan.refdb import RefDB, revcomp, _read_fasta_ref, data_dir
from gutscan import classify as C, targets as T

random.seed(11)


def mutate(s, sub=0.006, indel=0.001):
    out = []
    for ch in s:
        r = random.random()
        if r < sub:
            out.append(random.choice([b for b in "ACGT" if b != ch]))
        elif r < sub + indel / 2:
            continue
        elif r < sub + indel:
            out.append(ch); out.append(random.choice("ACGT"))
        else:
            out.append(ch)
    return "".join(out)


def window(s, L):
    a = random.randint(min(200, len(s) - L - 1), max(min(200, len(s) - L - 1), len(s) - L - 100))
    r = mutate(s[a:a + L])
    return r if random.random() < 0.5 else revcomp(r)


def truth(lin):
    fam, gen, sp = lin[3], lin[4], lin[5]
    return {t.name for t in T.TARGETS if T.matches(t, fam, gen, gen, sp)}


def pred(call):
    return {t.name for t in T.TARGETS if T.matches(t, call.family, call.genus, call.species_genus, call.species)} \
        if call.level != "none" else set()


def evaluate(db, items, label):
    tp, fp, fn = defaultdict(int), defaultdict(int), defaultdict(int)
    unclassified_target = 0; n_target_reads = 0; n = 0; t0 = time.time()
    for seq, lin in items:
        c = C.classify_read(db, seq)
        tr, pr = truth(lin), pred(c)
        n += 1
        if tr:
            n_target_reads += 1
            if c.level == "none": unclassified_target += 1
        for x in tr & pr: tp[x] += 1
        for x in pr - tr: fp[x] += 1
        for x in tr - pr: fn[x] += 1
    print(f"\n=== {label}: {n} reads ({n_target_reads} from target taxa), {n/(time.time()-t0):.0f} reads/s")
    print(f"{'target':34s}{'TP':>6s}{'FP':>6s}{'FN':>6s}{'precision':>10s}{'recall':>8s}")
    rows = []
    for t in T.TARGETS:
        a, b, c_ = tp[t.name], fp[t.name], fn[t.name]
        if a + b + c_ == 0: continue
        pr_ = a / (a + b) if a + b else float("nan"); rc = a / (a + c_) if a + c_ else float("nan")
        nm = ("  " + t.name) if t.rank == "species" else t.name
        print(f"{nm:34s}{a:6d}{b:6d}{c_:6d}{pr_:10.3f}{rc:8.3f}")
    TP, FP, FN = sum(tp[t.name] for t in T.TARGETS if t.rank != "species"), sum(fp[t.name] for t in T.TARGETS if t.rank != "species"), sum(fn[t.name] for t in T.TARGETS if t.rank != "species")
    print(f"{'ALL genus/family rows':34s}{TP:6d}{FP:6d}{FN:6d}{TP/(TP+FP):10.3f}{TP/(TP+FN):8.3f}")
    return TP, FP, FN


def main():
    ids, lins, seqs = [], [], []
    _read_fasta_ref(os.path.join(data_dir(), "reference.fa.gz"), ids, lins, seqs)
    idx = list(range(len(ids))); random.shuffle(idx)
    hold = set(idx[: int(0.15 * len(idx))])
    keep = [i for i in range(len(ids)) if i not in hold]
    db = RefDB([ids[i] for i in keep], [lins[i] for i in keep], [seqs[i] for i in keep])
    # held-out: all target-taxon refs (up to 700) + 700 random non-target refs
    ht = [i for i in hold if truth(lins[i])]; hn = [i for i in hold if not truth(lins[i])]
    random.shuffle(ht); random.shuffle(hn)
    items = []
    for i in ht[:700] + hn[:700]:
        for L in (250, 450):
            items.append((window(seqs[i], L), lins[i]))
    A = evaluate(db, items, "A. HELD-OUT strains (harsh)")
    # seen strains
    st = [i for i in keep if truth(lins[i])]; sn = [i for i in keep if not truth(lins[i])]
    random.shuffle(st); random.shuffle(sn)
    items = [(window(seqs[i], L), lins[i]) for i in st[:700] + sn[:700] for L in (250,)]
    B = evaluate(db, items, "B. SEEN strains (reads with sequencing error)")
    # junk
    junk = ["".join(random.choice("ACGT") for _ in range(300)) for _ in range(200)]
    wrong = sum(1 for s in junk if C.classify_read(db, s).level != "none")
    print(f"\nC. random DNA: {wrong}/200 wrongly classified")
    pa, pb, rb = A[0] / (A[0] + A[1]), B[0] / (B[0] + B[1]), B[0] / (B[0] + B[2])
    ok = pa >= 0.95 and pb >= 0.98 and rb >= 0.90 and wrong == 0
    print(f"\nACCURACY GATE (held-out precision>=0.95, seen precision>=0.98, seen recall>=0.90, junk=0): {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
