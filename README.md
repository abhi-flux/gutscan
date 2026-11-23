# GutScan - stool-sample 16S screener (desktop app)

Drop in a stool-sample FASTA -> GutScan cleans it -> aligns every sequence against a reference database ->
tells you which of the gut-brain-axis microbes are present and how much (% of the sample).

Target microbes: Faecalibacterium (F. prausnitzii), Roseburia, Eubacterium, Akkermansia (A. muciniphila), Bacteroides,
Lactobacillus (L. reuteri/plantarum/rhamnosus), Bifidobacterium (B. longum/breve/infantis), Enterococcus, Streptococcus,
Escherichia (E. coli), Bacillus, Clostridium, Coprococcus, Dialister, Prevotella, Lachnospiraceae, Desulfovibrio,
Enterobacteriaceae, C. difficile, Eggerthella (E. lenta), Ruminococcus (R. gnavus).

## Input it expects
16S rRNA amplicon sequences (reads, ASVs or OTUs) in FASTA (FASTQ and .gz also work).
NOT suitable: whole-genome / shotgun metagenome files.

## Run it (Windows, cmd.exe) - one command at a time
    pip install -r requirements.txt
    python main.py

Command-line (no window):  `python -m gutscan.cli sample.fasta --csv result.csv`

## Get a double-click .exe (no Python needed) - built on GitHub, no local tools
1. Put this folder in a GitHub repository (upload all files, including `.github`).
2. GitHub -> Actions -> "Build GutScan (Windows)" -> Run workflow.
3. When it finishes, download the `GutScan-windows` artifact, unzip, run `GutScan\GutScan.exe`.
The build first downloads real NCBI RefSeq 16S sequences for every target taxon, prints how many references each microbe has,
runs the accuracy benchmark (the build FAILS if accuracy drops below the gate), runs the tests, then builds the .exe.
The benchmark report is included in the download. See "Not verified" below.

## How it works
1. **Cleaning** - upper-case, U->T, strip gaps/digits, IUPAC -> N, trim N ends; drop too short (<150 bp), too long (>2500 bp),
   >3 % N, low-complexity, homopolymer artifacts, non-DNA; merge identical sequences (`;size=N` honoured).
2. **Search** - sampled canonical 11-mer index finds the 12 most similar reference sequences (either strand).
3. **Exact alignment** - each candidate is aligned end-to-end (edlib); % identity = 1 - edits / read length.
4. **Taxonomy rules** - species hint >= 98.7 %, genus >= 94.5 %, family >= 86.5 % identity (standard 16S cut-offs);
   the label must be shared by >= 80 % of hits within 0.5 % of the best, else it is lowered a rank. < 86.5 % = Unclassified.
5. **Presence call** - "Detected" = at least N reads (default 3) AND >= 0.01 % of the sample; fewer reads = "Trace".

## Measured accuracy (simulated reads, see docs/validation_results.txt; reproduce with `python tools/validate.py`)
| Test | Precision | Recall |
|---|---|---|
| Seen strains + sequencing errors (250 bp) | 99.9 % | 97.1 % |
| Held-out strains (never in the database; harsh) | 98.2 % | 82.6 % |
| Random DNA | 0 / 200 wrongly classified | |
Misses are mostly short reads that cannot separate close genera (e.g. Roseburia vs Blautia); those are reported at
family level instead of guessing. These numbers are from *simulated* data labelled with the database's own taxonomy -
they are not a substitute for testing on real, independently-identified samples.

## Limits (please read)
- % = relative abundance of 16S reads, not absolute counts; 16S copy number differences are not corrected.
- Species rows are "putative": 16S cannot reliably separate E. coli / Shigella, L. rhamnosus / L. casei, B. longum / B. infantis.
  Reference has no L. rhamnosus, B. infantis or E. rectale labels, and only 1 C. difficile sequence - so these rows cannot be called yet.
  Fix: put NCBI RefSeq 16S FASTA (headers like `>NR_037122.1 Lactobacillus rhamnosus ...`) at
  `C:\Users\<you>\GutScan\extra_references.fasta` (or next to the exe). It is loaded automatically.
- Reference: Greengenes 13_8, 97 % OTUs (2013 taxonomy; some labels are outdated, e.g. Clostridium is polyphyletic).
  Licence CC BY-SA 3.0 (`gutscan/data/LICENSE_greengenes.txt`). I could not download NCBI/SILVA from my build environment.
- Research/educational use only - not a medical diagnostic.

## Not verified
- The NCBI download step (NCBI is blocked from my environment, so it has never run live; loader/renaming logic is unit-tested with fake data).
- The Windows build (PyInstaller) - I can only run Linux here. The app itself was run and screenshotted headless on Linux.
- Accuracy on real stool data.

## Files
`gutscan/` engine + GUI, `gutscan/data/reference.fa.gz` database (24,926 sequences), `tools/build_reference.py` rebuilds it,
`tools/validate.py` benchmark, `tools/make_demo_sample.py` + `sample_stool_DEMO_synthetic.fasta` (synthetic demo file).
