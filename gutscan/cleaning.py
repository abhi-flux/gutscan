"""FASTA/FASTQ reading and sequence cleaning."""
import gzip, re
from collections import Counter
from dataclasses import dataclass, field

IUPAC_AMBIG = set("RYSWKMBDHVN")
_DNA_OK = set("ACGTUN") | IUPAC_AMBIG
_SIZE = re.compile(r"(?:^|[;\s_|])size=(\d+)")
_STRIP = re.compile(r"[\s\-\.\*\d]")


@dataclass
class CleanParams:
    min_len: int = 150
    max_len: int = 2500
    max_n_frac: float = 0.03
    use_size_annotations: bool = True


@dataclass
class CleanStats:
    reads_in: int = 0
    reads_kept: int = 0
    unique_kept: int = 0
    dropped: Counter = field(default_factory=Counter)
    weighted_in: int = 0
    weighted_kept: int = 0


def _open(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def read_records(path):
    """Yield (header, sequence) from FASTA or FASTQ (plain or .gz)."""
    with _open(path) as f:
        first = ""
        for line in f:
            if line.strip():
                first = line
                break
        if not first:
            return
        if first.startswith("@"):                      # FASTQ
            hdr = first[1:].strip()
            while True:
                seq = f.readline()
                if not seq:
                    break
                f.readline(); f.readline()
                yield hdr, seq.strip()
                nxt = f.readline()
                while nxt and not nxt.strip():
                    nxt = f.readline()
                if not nxt:
                    break
                hdr = nxt[1:].strip()
            return
        if not first.startswith(">"):
            raise ValueError("This does not look like a FASTA file (first line should start with '>').")
        hdr, buf = first[1:].strip(), []
        for line in f:
            if line.startswith(">"):
                yield hdr, "".join(buf)
                hdr, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        yield hdr, "".join(buf)


def clean_sequence(seq, p: CleanParams):
    """Return (clean_seq or None, reason_if_dropped)."""
    s = _STRIP.sub("", seq.upper()).replace("U", "T")
    if not s:
        return None, "empty sequence"
    bad = sum(1 for c in set(s) if c not in _DNA_OK)
    if bad and sum(s.count(c) for c in set(s) if c not in _DNA_OK) / len(s) > 0.10:
        return None, "not DNA (protein or text)"
    s = "".join(c if c in "ACGT" else "N" for c in s)   # IUPAC / junk -> N
    s = s.strip("N")
    if len(s) < p.min_len:
        return None, f"too short (<{p.min_len} bp)"
    if len(s) > p.max_len:
        return None, f"too long (>{p.max_len} bp - genome/contig, not a 16S read)"
    if s.count("N") / len(s) > p.max_n_frac:
        return None, "too many ambiguous bases (N)"
    top = max(s.count(b) for b in "ACGT") / len(s)
    if top > 0.60:
        return None, "low complexity"
    if re.search(r"([ACGT])\1{29,}", s):
        return None, "long homopolymer (artifact)"
    return s, ""


def clean_file(path, params: CleanParams = None, progress=None):
    """Returns (unique dict {seq: count}, stats)."""
    p = params or CleanParams()
    st = CleanStats()
    uniq = {}
    for n, (hdr, seq) in enumerate(read_records(path), 1):
        w = 1
        if p.use_size_annotations:
            m = _SIZE.search(hdr)
            if m:
                w = max(1, int(m.group(1)))
        st.reads_in += 1
        st.weighted_in += w
        s, why = clean_sequence(seq, p)
        if s is None:
            st.dropped[why] += w
            continue
        uniq[s] = uniq.get(s, 0) + w
        st.reads_kept += 1
        st.weighted_kept += w
        if progress and n % 5000 == 0:
            progress(n)
    st.unique_kept = len(uniq)
    if st.reads_in == 0:
        raise ValueError("No sequences found in this file.")
    return uniq, st
