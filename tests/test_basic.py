import os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from gutscan.cleaning import clean_sequence, CleanParams
from gutscan.refdb import RefDB, revcomp
from gutscan import classify as C

_db = None
def db():
    global _db
    if _db is None: _db = RefDB.load()
    return _db


def test_cleaning():
    p = CleanParams()
    s = "acgu" * 60
    assert clean_sequence(s, p)[0] == "ACGT" * 60
    assert clean_sequence("ACGT" * 5, p)[0] is None
    assert clean_sequence("MKVLAAGIVGLLLAQPSTAWAEEKKRHH" * 8, p)[1].startswith("not DNA")
    assert clean_sequence("A" * 300, p)[0] is None


def test_known_genus_both_strands():
    d = db()
    i = next(k for k, l in enumerate(d.lineages) if l[4] == "Faecalibacterium")
    seq = d.seqs[i][300:600]
    for q in (seq, revcomp(seq)):
        c = C.classify_read(d, q)
        assert c.genus == "Faecalibacterium" and c.identity > 0.99


def test_random_dna_unclassified():
    random.seed(1)
    d = db()
    for _ in range(20):
        s = "".join(random.choice("ACGT") for _ in range(300))
        assert C.classify_read(d, s).level == "none"


def test_ncbi_header_renames():
    from gutscan.refdb import parse_ncbi_header as P
    assert P("NR_1.1 Lacticaseibacillus rhamnosus strain ATCC 7469 16S ribosomal RNA") == ("NR_1.1", "Lactobacillus", "rhamnosus", "Lactobacillaceae")
    assert P("NR_2.1 Clostridioides difficile strain X 16S")[1:] == ("Clostridium", "difficile", "Peptostreptococcaceae")
    assert P("NR_3.1 Agathobacter rectalis strain X 16S")[1:3] == ("Eubacterium", "rectale")
    assert P("NR_4.1 Bifidobacterium longum subsp. infantis ATCC 15697 16S")[2] == "infantis"
    assert P("NR_5.1 Candidatus Foo bar 16S") is None
