"""Read classification: k-mer prefilter -> exact alignment -> identity-based taxonomy.

Thresholds follow the widely used 16S rRNA identity cut-offs (Yarza et al. 2014):
  species >= 98.7 %, genus >= 94.5 %, family >= 86.5 %.
Everything below LOW_FLOOR is reported as "unclassified" rather than forced
onto the closest target.
"""
from collections import Counter
from dataclasses import dataclass, field
import edlib

from .refdb import RefDB, revcomp
from . import targets as T

SPECIES_ID, GENUS_ID, FAMILY_ID = 0.987, 0.945, 0.865
NEAR_BEST = 0.005          # hits within 0.5 % of the best form the consensus
CONSENSUS = 0.80           # share of near-best hits that must agree on a label
SPECIES_CLASH = 0.003      # a different species label this close to best blocks a species call
N_CANDIDATES = 12


@dataclass
class ReadCall:
    identity: float = 0.0
    level: str = "none"              # species | genus | family | none
    family: str = ""
    genus: str = ""
    phylum: str = ""
    species_genus: str = ""
    species: str = ""
    best_ref: str = ""


@dataclass
class Params:
    min_reads: int = 3               # reads needed to call a microbe "detected"
    min_abundance: float = 0.0001    # ... and at least this fraction of the sample (0.01 %)


def _align(q, t):
    r = edlib.align(q, t, mode="HW", task="distance")
    return r["editDistance"]


def classify_read(db: RefDB, seq: str) -> ReadCall:
    cands = db.candidates(seq, N_CANDIDATES)
    if not cands:
        return ReadCall()
    # decide strand using the top candidate
    top = db.seqs[cands[0][0]]
    rc = revcomp(seq)
    q = seq if _align(seq, top) <= _align(rc, top) else rc
    L = len(q)
    hits = []
    for ri, _ in cands:
        ed = _align(q, db.seqs[ri])
        hits.append((1.0 - ed / L, ri))
    hits.sort(key=lambda x: -x[0])
    best = hits[0][0]
    call = ReadCall(identity=best, best_ref=db.ids[hits[0][1]])
    if best < FAMILY_ID:
        return call
    near = [ri for idt, ri in hits if idt >= best - NEAR_BEST]

    def consensus(rank_idx):
        c = Counter(db.lineages[ri][rank_idx] for ri in near)
        lab, n = c.most_common(1)[0]
        return lab if lab and n / len(near) >= CONSENSUS else ""

    call.phylum = consensus(0)
    call.family = consensus(3)
    call.level = "family"
    if best >= GENUS_ID:
        call.genus = consensus(4)
        if call.genus:
            call.level = "genus"
    if best >= SPECIES_ID:
        blri = hits[0][1]
        sp, gen = db.lineages[blri][5], db.lineages[blri][4]
        if sp and gen:
            clash = any(db.lineages[ri][5] and db.lineages[ri][5] != sp and
                        (db.lineages[ri][4] or gen) == gen
                        for idt, ri in hits[1:] if idt >= best - SPECIES_CLASH)
            if not clash:
                call.species, call.species_genus = sp, gen
                if call.level == "family" and not call.genus:
                    call.genus = gen if db.lineages[blri][3] == call.family else ""
                    if call.genus:
                        call.level = "genus"
    return call


# --------------------------------------------------------------------------
@dataclass
class TargetResult:
    target: T.Target
    reads: int = 0               # reads matching this taxon (inclusive of sub-groups)
    exclusive: int = 0           # reads whose most specific target is this one
    identity_sum: float = 0.0
    ref_count: int = 0
    pct_sample: float = 0.0      # % of all cleaned reads
    pct_classified: float = 0.0  # % of reads classified at family level or better
    mean_identity: float = 0.0
    status: str = "Not detected"


@dataclass
class Result:
    targets: list
    n_total: int = 0
    n_classified: int = 0
    n_unclassified: int = 0
    n_other: int = 0
    other_genera: Counter = field(default_factory=Counter)
    other_families: Counter = field(default_factory=Counter)
    unique_seqs: int = 0
    cancelled: bool = False


def run_classification(db: RefDB, uniq: dict, params: Params = None, progress=None, cancel=None) -> Result:
    p = params or Params()
    res = Result(targets=[TargetResult(t, ref_count=_ref_count(db, t)) for t in T.TARGETS])
    by_name = {r.target.name: r for r in res.targets}
    res.unique_seqs = len(uniq)
    items = sorted(uniq.items(), key=lambda kv: -kv[1])
    for n, (seq, w) in enumerate(items, 1):
        if cancel and cancel():
            res.cancelled = True
            break
        c = classify_read(db, seq)
        res.n_total += w
        if c.level == "none":
            res.n_unclassified += w
        else:
            res.n_classified += w
            matched = [by_name[t.name] for t in T.TARGETS
                       if T.matches(t, c.family, c.genus, c.species_genus, c.species)]
            if not matched:
                res.n_other += w
                res.other_genera[c.genus or f"({c.family or c.phylum or 'unresolved'})"] += w
                if c.family:
                    res.other_families[c.family] += w
            else:
                for r in matched:
                    r.reads += w
                    r.identity_sum += c.identity * w
                spec = [r for r in matched if r.target.rank == "species"]
                gen = [r for r in matched if r.target.rank == "genus"]
                (spec or gen or matched)[0].exclusive += w
        if progress and (n % 200 == 0 or n == len(items)):
            progress(n, len(items))
    for r in res.targets:
        if res.n_total:
            r.pct_sample = 100.0 * r.reads / res.n_total
        if res.n_classified:
            r.pct_classified = 100.0 * r.reads / res.n_classified
        if r.reads:
            r.mean_identity = 100.0 * r.identity_sum / r.reads
        frac = r.reads / res.n_total if res.n_total else 0
        if r.reads >= p.min_reads and frac >= p.min_abundance:
            r.status = "Detected"
        elif r.reads > 0:
            r.status = "Trace (below threshold)"
    return res


def _ref_count(db, t):
    if t.rank == "genus": return db.count_for(genus=t.genus)
    if t.rank == "family": return db.count_for(family=t.family)
    return db.count_for(genus=t.genus, species=t.species)
