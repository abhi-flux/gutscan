"""Download real 16S rRNA RefSeq sequences (NCBI) for every target taxon.

Writes gutscan/data/ncbi_refseq_16S.fasta.gz, which GutScan loads automatically on top of
the bundled Greengenes set. Run by the GitHub Actions build; you can also run it on your own PC:
    python tools/fetch_refseq.py
Optional: set NCBI_API_KEY (free from your NCBI account) for faster downloads.
"""
import gzip, json, os, sys, time, urllib.parse, urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = os.path.join(os.path.dirname(__file__), "..", "gutscan", "data", "ncbi_refseq_16S.fasta.gz")
ORGANISMS = ["Faecalibacterium", "Roseburia", "Eubacterium", "Agathobacter", "Akkermansia", "Bacteroides", "Phocaeicola",
             "Lactobacillus", "Lacticaseibacillus", "Lactiplantibacillus", "Limosilactobacillus", "Ligilactobacillus",
             "Bifidobacterium", "Enterococcus", "Streptococcus", "Escherichia", "Shigella", "Klebsiella", "Salmonella",
             "Bacillus", "Clostridium", "Clostridioides", "Coprococcus", "Dialister", "Prevotella", "Segatella",
             "Desulfovibrio", "Eggerthella", "Ruminococcus", "Mediterraneibacter", "Lachnospiraceae", "Enterobacteriaceae"]
KEY = os.environ.get("NCBI_API_KEY", "").strip()
DELAY = 0.12 if KEY else 0.4


def call(endpoint, **params):
    params.update(tool="gutscan", email="gutscan-build@users.noreply.github.com")
    if KEY:
        params["api_key"] = KEY
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(6):
        try:
            time.sleep(DELAY)
            with urllib.request.urlopen(urllib.request.Request(BASE + endpoint, data=data), timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            print(f"   retry {attempt + 1} ({e})", flush=True)
            time.sleep(2 ** attempt)
    raise SystemExit(f"NCBI request failed repeatedly: {endpoint}")


def main():
    seqs = {}
    for org in ORGANISMS:
        term = f'"{org}"[Organism] AND "16S ribosomal RNA"[Title] AND srcdb_refseq[PROP] AND biomol_rrna[PROP]'
        js = json.loads(call("esearch.fcgi", db="nuccore", term=term, retmax=0, usehistory="y", retmode="json"))["esearchresult"]
        n = int(js.get("count", 0))
        print(f"{org}: {n} RefSeq 16S records", flush=True)
        for start in range(0, n, 200):
            txt = call("efetch.fcgi", db="nuccore", query_key=js["querykey"], WebEnv=js["webenv"],
                       retstart=start, retmax=200, rettype="fasta", retmode="text")
            for rec in txt.split(">")[1:]:
                head, _, body = rec.partition("\n")
                acc = head.split()[0]
                if "16S" in head and acc not in seqs:
                    seqs[acc] = (head, "".join(body.split()))
    print("total unique sequences:", len(seqs))
    if len(seqs) < 1500:
        raise SystemExit("Too few sequences downloaded - NCBI query/format problem. Failing the build on purpose.")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with gzip.open(OUT, "wt", compresslevel=9) as o:
        for head, s in seqs.values():
            o.write(f">{head}\n{s}\n")
    print("wrote", os.path.abspath(OUT), os.path.getsize(OUT) // 1024, "KB")


if __name__ == "__main__":
    main()
