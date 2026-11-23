"""Reference database + k-mer prefilter index."""
import gzip, os, sys
import numpy as np

K = 11
SAMPLE = 3                 # keep ~1/3 of canonical k-mers (consistent for refs and reads)
MAX_POSTING = 2500         # k-mers present in more refs than this are uninformative
RANKS = ("phylum", "class", "order", "family", "genus", "species")

_LUT = np.full(256, 4, np.uint8)
for i, c in enumerate("ACGT"):
    _LUT[ord(c)] = i
    _LUT[ord(c.lower())] = i


def data_dir():
    if hasattr(sys, "_MEIPASS"):                      # PyInstaller bundle
        return os.path.join(sys._MEIPASS, "gutscan", "data")
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def encode(seq: str) -> np.ndarray:
    return _LUT[np.frombuffer(seq.encode("ascii", "replace"), np.uint8)]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def sampled_kmers(codes: np.ndarray) -> np.ndarray:
    n = len(codes) - K + 1
    if n <= 0:
        return np.empty(0, np.uint32)
    valid = np.ones(n, bool)
    fwd = np.zeros(n, np.uint32)
    rev = np.zeros(n, np.uint32)
    for j in range(K):
        c = codes[j:j + n]
        valid &= c < 4
        c32 = (c & 3).astype(np.uint32)
        fwd = (fwd << np.uint32(2)) | c32
        rev |= (np.uint32(3) - c32) << np.uint32(2 * j)
    canon = np.minimum(fwd, rev)[valid]
    h = (canon * np.uint32(2654435761)) >> np.uint32(16)
    return canon[(h % np.uint32(SAMPLE)) == 0]


class RefDB:
    def __init__(self, ids, lineages, seqs):
        self.ids = ids
        self.lineages = lineages      # list of 6-tuples phylum..species
        self.seqs = seqs
        self._build_index()

    # ---------------------------------------------------------------- loading
    @classmethod
    def load(cls, path=None, extra=None, progress=None):
        path = path or os.path.join(data_dir(), "reference.fa.gz")
        ids, lins, seqs = [], [], []
        _read_fasta_ref(path, ids, lins, seqs)
        db_genus_lin = {}
        for l in lins:
            if l[4] and l[4] not in db_genus_lin:
                db_genus_lin[l[4]] = l[:4]
        bundled = os.path.join(data_dir(), "ncbi_refseq_16S.fasta.gz")
        for f in (bundled, extra):
            if f and os.path.exists(f):
                _read_ncbi_style(f, db_genus_lin, ids, lins, seqs)
        return cls(ids, lins, seqs)

    # ------------------------------------------------------------------ index
    def _build_index(self):
        chunks_keys, chunks_ids = [], []
        n = len(self.seqs)
        B = 500
        for start in range(0, n, B):
            sl = self.seqs[start:start + B]
            lens = np.fromiter((len(s) for s in sl), np.int64, len(sl))
            joined = "N".join(sl)
            codes = encode(joined)
            km_pos_ref = np.repeat(np.arange(start, start + len(sl), dtype=np.int64), lens + 1)[:len(codes)]
            # recompute k-mers on whole chunk but need per-position ref id
            nk = len(codes) - K + 1
            valid = np.ones(nk, bool)
            fwd = np.zeros(nk, np.uint32); rev = np.zeros(nk, np.uint32)
            for j in range(K):
                c = codes[j:j + nk]
                valid &= c < 4
                c32 = (c & 3).astype(np.uint32)
                fwd = (fwd << np.uint32(2)) | c32
                rev |= (np.uint32(3) - c32) << np.uint32(2 * j)
            canon = np.minimum(fwd, rev)
            h = (canon * np.uint32(2654435761)) >> np.uint32(16)
            keep = valid & ((h % np.uint32(SAMPLE)) == 0)
            chunks_keys.append(canon[keep])
            chunks_ids.append(km_pos_ref[:nk][keep].astype(np.uint64))
        keys = np.concatenate(chunks_keys).astype(np.uint64)
        rid = np.concatenate(chunks_ids)
        comb = np.unique((keys << np.uint64(20)) | rid)
        keys = (comb >> np.uint64(20)).astype(np.uint32)
        rid = (comb & np.uint64((1 << 20) - 1)).astype(np.int32)
        ukeys, first, counts = np.unique(keys, return_index=True, return_counts=True)
        ok = counts <= MAX_POSTING
        self.ukeys = ukeys[ok]
        self.starts = first[ok]
        self.counts = counts[ok]
        self.postings = rid

    # ------------------------------------------------------------------ query
    def candidates(self, seq: str, top=12):
        """Return list of (ref_index, votes) best-first."""
        q = sampled_kmers(encode(seq))
        if len(q) == 0:
            return []
        q = np.unique(q)
        pos = np.searchsorted(self.ukeys, q)
        pos[pos >= len(self.ukeys)] = 0
        hit = self.ukeys[pos] == q
        pos = pos[hit]
        if len(pos) == 0:
            return []
        parts = [self.postings[self.starts[p]:self.starts[p] + self.counts[p]] for p in pos]
        allr = np.concatenate(parts)
        refs, votes = np.unique(allr, return_counts=True)
        order = np.argsort(-votes, kind="stable")[:top]
        return [(int(refs[i]), int(votes[i])) for i in order]

    def count_for(self, genus="", family="", species=""):
        n = 0
        for l in self.lineages:
            if genus and l[4] != genus: continue
            if family and l[3] != family: continue
            if species and l[5] != species: continue
            n += 1
        return n


def _read_fasta_ref(path, ids, lins, seqs):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        cur, buf = None, []
        def flush():
            if cur is None: return
            parts = cur.split("|")
            lin = tuple((parts[1:7] + [""] * 6)[:6])
            ids.append(parts[0]); lins.append(lin); seqs.append("".join(buf).upper())
        for line in f:
            if line[0] == ">":
                flush(); cur, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        flush()


# NCBI renamed many genera after the Greengenes 2013 taxonomy. Map new -> (old genus used by GutScan, family).
SYNONYMS = {
    "Lacticaseibacillus": ("Lactobacillus", "Lactobacillaceae"), "Lactiplantibacillus": ("Lactobacillus", "Lactobacillaceae"),
    "Limosilactobacillus": ("Lactobacillus", "Lactobacillaceae"), "Ligilactobacillus": ("Lactobacillus", "Lactobacillaceae"),
    "Lentilactobacillus": ("Lactobacillus", "Lactobacillaceae"), "Levilactobacillus": ("Lactobacillus", "Lactobacillaceae"),
    "Latilactobacillus": ("Lactobacillus", "Lactobacillaceae"), "Companilactobacillus": ("Lactobacillus", "Lactobacillaceae"),
    "Clostridioides": ("Clostridium", "Peptostreptococcaceae"),
    "Agathobacter": ("Eubacterium", "Lachnospiraceae"),
    "Mediterraneibacter": ("Ruminococcus", "Lachnospiraceae"),
    "Phocaeicola": ("Bacteroides", "Bacteroidaceae"),
    "Segatella": ("Prevotella", "Prevotellaceae"),
}
SPECIES_FIX = {("Eubacterium", "rectalis"): "rectale"}


def parse_ncbi_header(hdr):
    """'NR_1.1 Lacticaseibacillus rhamnosus strain X 16S ...' -> (acc, genus, species, family_override)."""
    toks = hdr.replace("[", "").replace("]", "").split()
    if len(toks) < 3 or toks[1] in ("Candidatus", "uncultured", "Uncultured"):
        return None
    gen, sp = toks[1], toks[2]
    if sp in ("sp.", "sp", "bacterium", "cf.", "aff."):
        sp = ""
    fam = ""
    if "subsp." in toks:                       # e.g. Bifidobacterium longum subsp. infantis
        k = toks.index("subsp.")
        if k + 1 < len(toks) and toks[k + 1] == "infantis":
            sp = "infantis"
    if gen in SYNONYMS:
        gen, fam = SYNONYMS[gen]
    sp = SPECIES_FIX.get((gen, sp), sp)
    return toks[0], gen, sp, fam


def _read_ncbi_style(path, genus_lin, ids, lins, seqs):
    """Extra references with NCBI headers:  >NR_037122.1 Lactobacillus rhamnosus strain ..."""
    from .cleaning import read_records
    seen = set(ids)
    for hdr, seq in read_records(path):
        r = parse_ncbi_header(hdr)
        if not r:
            continue
        acc, gen, sp, fam = r
        if acc in seen:
            continue
        s = "".join(c for c in seq.upper().replace("U", "T") if c in "ACGTN")
        if len(s) < 900 or s.count("N") / len(s) > 0.01:
            continue
        up = genus_lin.get(gen, ("", "", "", ""))
        seen.add(acc)
        ids.append(acc); lins.append((up[0], up[1], up[2], fam or up[3], gen, sp)); seqs.append(s)
