import csv
from . import targets as T


def to_text(res, st):
    L = []
    L.append(f"Sequences in file: {st.reads_in}  | kept after cleaning: {st.reads_kept}  | "
             f"unique: {st.unique_kept}")
    for why, n in st.dropped.most_common():
        L.append(f"  dropped {n}: {why}")
    L.append(f"Classified (>=86.5% identity): {res.n_classified} / {res.n_total}   "
             f"unclassified: {res.n_unclassified}   other (non-target): {res.n_other}")
    L.append("")
    L.append(f"{'Microbe':34s}{'Status':26s}{'Reads':>8s}{'% sample':>10s}{'Ident%':>8s}{'Refs':>6s}")
    for g in T.GROUPS:
        L.append(f"-- {g}")
        for r in res.targets:
            if r.target.group != g: continue
            nm = ("   " + r.target.name) if r.target.rank == "species" else r.target.name
            L.append(f"{nm:34s}{r.status:26s}{r.reads:8d}{r.pct_sample:10.3f}"
                     f"{r.mean_identity:8.1f}{r.ref_count:6d}")
    L.append("")
    L.append("Top other genera: " + ", ".join(f"{k} {v}" for k, v in res.other_genera.most_common(10)))
    return "\n".join(L)


def to_csv(res, st, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["microbe", "rank", "group", "status", "reads", "percent_of_sample",
                    "percent_of_classified", "mean_identity_percent", "reference_sequences", "role", "note"])
        for r in res.targets:
            t = r.target
            w.writerow([t.name, t.rank, t.group, r.status, r.reads, f"{r.pct_sample:.4f}",
                        f"{r.pct_classified:.4f}", f"{r.mean_identity:.2f}", r.ref_count, t.role, t.note])
        w.writerow([])
        w.writerow(["summary", "total_reads", res.n_total, "classified", res.n_classified,
                    "unclassified", res.n_unclassified, "other_non_target", res.n_other])
        for why, n in st.dropped.items():
            w.writerow(["dropped_in_cleaning", why, n])
        for g, n in res.other_genera.most_common(50):
            w.writerow(["other_genus", g, n])
