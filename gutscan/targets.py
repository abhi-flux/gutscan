"""The microbes GutScan is trained to look for (taken from the project brief).

A Target is matched against the taxonomy GutScan assigns to each read:
  * rank "genus"   -> read's genus equals `genus`
  * rank "family"  -> read's family equals `family`
  * rank "species" -> read's genus AND species equal `genus` / `species`
Species rows are *sub-rows* of a genus row (parent=...). 16S cannot always
separate species, so species calls are only made at >=98.7 % identity and are
labelled "putative" in the app.
"""
from dataclasses import dataclass

SCFA = "SCFA producers"
NEURO = "Neurotransmitter / neuroactive producers"
MENTAL = "Linked to mental-health & neurological conditions"
HARM = "Potentially harmful / pro-inflammatory"

GROUPS = [SCFA, NEURO, MENTAL, HARM]


@dataclass(frozen=True)
class Target:
    name: str
    group: str
    rank: str                 # genus | family | species
    genus: str = ""
    species: str = ""
    family: str = ""
    role: str = ""
    note: str = ""
    parent: str = ""          # name of parent row (species rows only)


def _g(name, group, role, note=""):
    return Target(name, group, "genus", genus=name, role=role, note=note)


def _s(genus, species, group, parent, role="", note=""):
    return Target(f"{genus} {species}", group, "species", genus=genus,
                  species=species, parent=parent, role=role, note=note)


TARGETS = [
    # ---- SCFA producers -------------------------------------------------
    _g("Faecalibacterium", SCFA, "Butyrate producer; anti-inflammatory; supports gut barrier; depleted in depression"),
    _s("Faecalibacterium", "prausnitzii", SCFA, "Faecalibacterium"),
    _g("Roseburia", SCFA, "Butyrate producers; linked to reduced neuroinflammation"),
    _g("Eubacterium", SCFA, "Butyrate production (E. rectale), mucosal health",
       "E. rectale is also reported as a putative species row (needs the NCBI RefSeq references)."),
    _s("Eubacterium", "rectale", SCFA, "Eubacterium",
       note="Now classified as Agathobacter rectalis by NCBI; GutScan keeps the older name."),
    _g("Akkermansia", SCFA, "Mucin degrader; strengthens gut barrier, lowers systemic inflammation; altered in Parkinson's/Alzheimer's"),
    _s("Akkermansia", "muciniphila", SCFA, "Akkermansia"),
    _g("Bacteroides", SCFA, "Propionate and acetate producers; some also make GABA"),
    # ---- Neurotransmitter producers ------------------------------------
    _g("Lactobacillus", NEURO, "GABA, acetylcholine, vagus-mediated effects on anxiety/stress (L. rhamnosus, plantarum, reuteri)"),
    _s("Lactobacillus", "reuteri", NEURO, "Lactobacillus"),
    _s("Lactobacillus", "plantarum", NEURO, "Lactobacillus"),
    _s("Lactobacillus", "rhamnosus", NEURO, "Lactobacillus",
       note="16S cannot reliably separate L. rhamnosus from L. casei/paracasei."),
    _g("Bifidobacterium", NEURO, "GABA, tryptophan availability, reduced HPA-axis stress response (B. longum, breve, infantis)"),
    _s("Bifidobacterium", "longum", NEURO, "Bifidobacterium"),
    _s("Bifidobacterium", "breve", NEURO, "Bifidobacterium"),
    _s("Bifidobacterium", "infantis", NEURO, "Bifidobacterium",
       note="B. infantis is a subspecies of B. longum; identical by 16S."),
    _g("Enterococcus", NEURO, "Serotonin and dopamine precursors"),
    _g("Streptococcus", NEURO, "Serotonin and dopamine precursors"),
    _g("Escherichia", NEURO, "Norepinephrine and dopamine",
       "Escherichia/Shigella/Citrobacter are near-identical by 16S; often only resolvable to Enterobacteriaceae."),
    _s("Escherichia", "coli", NEURO, "Escherichia"),
    _g("Bacillus", NEURO, "Norepinephrine and dopamine"),
    _g("Clostridium", NEURO, "Spore-formers; stimulate enterochromaffin cells (~90% of body serotonin)",
       "Clostridium is polyphyletic - this row includes every read labelled Clostridium."),
    # ---- Mental-health / neurological genera ----------------------------
    _g("Coprococcus", MENTAL, "Consistently depleted in depression; tied to quality of life; dopamine metabolism"),
    _g("Dialister", MENTAL, "Consistently depleted in depression; tied to quality of life"),
    _g("Prevotella", MENTAL, "Elevated in Parkinson's in some studies; context-dependent"),
    Target("Lachnospiraceae", MENTAL, "family", family="Lachnospiraceae",
           role="SCFA producers; reduced in several neuropsychiatric conditions",
           note="Family row: includes Roseburia, Coprococcus, Eubacterium, Ruminococcus gnavus etc."),
    # ---- Potentially harmful --------------------------------------------
    _g("Desulfovibrio", HARM, "Hydrogen sulfide and LPS; linked to Parkinson's and autism"),
    Target("Enterobacteriaceae", HARM, "family", family="Enterobacteriaceae",
           role="LPS-driven inflammation; increased gut and blood-brain barrier permeability (e.g. Klebsiella, Escherichia)",
           note="Family row: includes Escherichia, Klebsiella, etc."),
    _s("Clostridium", "difficile", HARM, "Clostridium", "Neurotoxic metabolites such as p-cresol",
       "Only 1 reference sequence in the bundled database - add NCBI RefSeq sequences for better coverage."),
    _g("Eggerthella", HARM, "Elevated in depression; linked to inflammation"),
    _s("Eggerthella", "lenta", HARM, "Eggerthella"),
    _g("Ruminococcus", HARM, "Inflammatory polysaccharides (R. gnavus); linked to mood disorders and autism",
       "Genus row includes non-gnavus Ruminococcus."),
    _s("Ruminococcus", "gnavus", HARM, "Ruminococcus"),
]

TARGET_BY_NAME = {t.name: t for t in TARGETS}


def matches(t: Target, family: str, genus: str, species_genus: str, species: str) -> bool:
    if t.rank == "genus":
        return genus == t.genus or (t.genus == "Escherichia" and species_genus == "Escherichia")
    if t.rank == "family":
        return family == t.family
    return species_genus == t.genus and species == t.species
