"""Print how many reference sequences each target has; warn about empty species rows."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from gutscan.refdb import RefDB, data_dir
from gutscan import targets as T
db = RefDB.load()
print("total reference sequences:", len(db.seqs))
has_ncbi = os.path.exists(os.path.join(data_dir(), "ncbi_refseq_16S.fasta.gz"))
print("NCBI RefSeq file bundled:", has_ncbi)
empty = []
for t in T.TARGETS:
    n = db.count_for(genus=t.genus, species=t.species) if t.rank == "species" else (
        db.count_for(genus=t.genus) if t.rank == "genus" else db.count_for(family=t.family))
    print(f"  {t.name:36s}{n:6d}")
    if n == 0: empty.append(t.name)
if empty:
    print("WARNING: no reference sequences for:", ", ".join(empty))
if os.environ.get("REQUIRE_NCBI") and not has_ncbi:
    raise SystemExit("NCBI RefSeq file missing")
