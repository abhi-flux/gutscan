"""GutScan desktop window (PySide6)."""
import os, sys, traceback
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QBrush, QPen
from PySide6.QtWidgets import (QApplication, QFileDialog, QHBoxLayout, QLabel, QMainWindow, QMessageBox,
                               QProgressBar, QPushButton, QSpinBox, QTabWidget, QTextEdit, QTreeWidget,
                               QTreeWidgetItem, QVBoxLayout, QWidget, QScrollArea, QHeaderView, QFrame)

from . import __version__, targets as T
from .cleaning import clean_file, CleanParams
from .refdb import RefDB
from .classify import run_classification, Params
from .report import to_csv

ACCENT = "#f2c200"
STYLE = f"""
QWidget {{ background:#0e0f12; color:#d7d9de; font-size:13px; }}
QLabel#title {{ color:{ACCENT}; font-size:22px; font-weight:700; letter-spacing:3px; }}
QLabel#sub {{ color:#7c8089; }}
QFrame#card {{ background:#16181d; border:1px solid #262a31; border-radius:6px; }}
QLabel#cardnum {{ color:{ACCENT}; font-size:20px; font-weight:700; background:transparent; }}
QLabel#cardlbl {{ color:#7c8089; font-size:11px; background:transparent; }}
QPushButton {{ background:#1d2026; border:1px solid #333842; border-radius:5px; padding:7px 16px; }}
QPushButton:hover {{ border-color:{ACCENT}; }}
QPushButton#go {{ background:{ACCENT}; color:#0e0f12; font-weight:700; border:none; }}
QPushButton#go:disabled {{ background:#4a4326; color:#8a8467; }}
QPushButton:disabled {{ color:#5b5f68; }}
QTreeWidget, QTextEdit {{ background:#121418; border:1px solid #262a31; alternate-background-color:#15171c; }}
QTreeWidget::item {{ padding:3px; }}
QTreeWidget::item:selected {{ background:#2a2614; color:{ACCENT}; }}
QHeaderView::section {{ background:#1a1d23; color:{ACCENT}; border:none; padding:5px; font-weight:600; }}
QTabWidget::pane {{ border:1px solid #262a31; }}
QTabBar::tab {{ background:#16181d; padding:8px 18px; border:1px solid #262a31; border-bottom:none; }}
QTabBar::tab:selected {{ color:{ACCENT}; background:#1d2026; }}
QProgressBar {{ border:1px solid #333842; border-radius:4px; text-align:center; height:16px; background:#16181d; }}
QProgressBar::chunk {{ background:{ACCENT}; }}
QSpinBox {{ background:#16181d; border:1px solid #333842; padding:3px; }}
"""
EXTRA_REFS = [os.path.join(os.path.expanduser("~"), "GutScan", "extra_references.fasta"),
              os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "extra_references.fasta")]

ABOUT = f"""GutScan {__version__} - 16S rRNA screener

WHAT IT DOES
Reads a stool-sample FASTA (16S rRNA amplicon reads / ASVs / OTUs), cleans it, aligns every distinct sequence against a
reference database and reports which of the listed gut-brain-axis microbes are present and what share of the sample they make up.

HOW IT DECIDES (accuracy rules)
 - Each sequence is aligned (both strands) to its closest reference sequences; % identity is computed over the whole read.
 - Genus call needs >= 94.5 % identity, family call >= 86.5 %, species hint >= 98.7 %  (standard 16S cut-offs).
 - The label must be agreed by >= 80 % of the near-best hits; otherwise the call is lowered to family. It never guesses.
 - Below 86.5 % identity a sequence is "Unclassified" - it is NOT forced onto the closest target.
 - A microbe is "Detected" only with >= N reads (default 3) AND >= 0.01 % of the sample. Fewer = "Trace".

LIMITS YOU SHOULD KNOW
 - % values are RELATIVE abundance of 16S reads (not absolute cell counts); 16S copy-number differences are not corrected.
 - 16S cannot reliably separate some species: Escherichia vs Shigella, L. rhamnosus vs L. casei, B. longum vs B. infantis.
   Species rows are therefore "putative"; genus/family rows are the dependable result.
 - Short reads (V4 etc.) may not separate close genera inside Lachnospiraceae; those reads are reported at family level.
 - Whole-genome / shotgun metagenome FASTA files are not suitable (too long/too much data) - this is a 16S tool.
 - Reference = Greengenes 13_8 (97 % OTUs, CC BY-SA 3.0), 2013 taxonomy. 'Eubacterium rectale' has no species label in it.
 - To add more references (e.g. NCBI RefSeq 16S), put a FASTA with NCBI-style headers (>NR_xxxx Genus species ...) in:
   {EXTRA_REFS[0]}
 - Research / educational use only. This is NOT a medical diagnostic.
"""


class Card(QFrame):
    def __init__(self, label):
        super().__init__(); self.setObjectName("card")
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 8, 12, 8)
        self.num = QLabel("-"); self.num.setObjectName("cardnum")
        lbl = QLabel(label); lbl.setObjectName("cardlbl")
        lay.addWidget(self.num); lay.addWidget(lbl)

    def set(self, v): self.num.setText(str(v))


class Bars(QWidget):
    def __init__(self):
        super().__init__(); self.rows = []; self.setMinimumHeight(60)

    def set_rows(self, rows):
        self.rows = rows; self.setMinimumHeight(30 + 26 * len(rows)); self.update()

    def paintEvent(self, _):
        p = QPainter(self); p.setRenderHint(QPainter.Antialiasing)
        if not self.rows:
            p.setPen(QColor("#7c8089")); p.drawText(20, 30, "Run an analysis to see the composition."); return
        w = self.width(); lx = 230; mx = max(r[1] for r in self.rows) or 1
        for i, (name, pct, col) in enumerate(self.rows):
            y = 12 + i * 26
            p.setPen(QColor("#d7d9de")); p.drawText(10, y, lx - 20, 20, Qt.AlignVCenter | Qt.AlignRight, name)
            bw = max(2, (w - lx - 90) * pct / mx)
            p.setPen(Qt.NoPen); p.setBrush(QBrush(QColor(col))); p.drawRoundedRect(lx, y + 3, bw, 14, 3, 3)
            p.setPen(QColor("#d7d9de")); p.drawText(int(lx + bw + 8), y, 80, 20, Qt.AlignVCenter, f"{pct:.2f} %")


class Worker(QThread):
    progress = Signal(str, int, int)
    done = Signal(object, object)
    failed = Signal(str)

    def __init__(self, path, clean_p, params, db):
        super().__init__(); self.path, self.clean_p, self.params, self.db = path, clean_p, params, db
        self._cancel = False

    def cancel(self): self._cancel = True

    def run(self):
        try:
            self.progress.emit("Cleaning sequences...", 0, 0)
            uniq, st = clean_file(self.path, self.clean_p)
            if not uniq:
                raise ValueError("No usable sequences left after cleaning. See the 'Cleaning log' for why reads were dropped.")
            if self.db[0] is None:
                self.progress.emit("Loading reference database (first run only)...", 0, 0)
                extra = next((p for p in EXTRA_REFS if os.path.exists(p)), None)
                self.db[0] = RefDB.load(extra=extra)
            res = run_classification(self.db[0], uniq, self.params,
                                     progress=lambda a, b: self.progress.emit("Classifying sequences...", a, b),
                                     cancel=lambda: self._cancel)
            self.done.emit(res, st)
        except Exception as e:
            self.failed.emit(f"{e}")


class Main(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"GutScan {__version__} - stool 16S microbe screener"); self.resize(1120, 760)
        self.setAcceptDrops(True)
        self.path = None; self.res = None; self.st = None; self.db = [None]; self.worker = None
        root = QWidget(); self.setCentralWidget(root); lay = QVBoxLayout(root)
        t = QLabel("GUTSCAN"); t.setObjectName("title")
        s = QLabel("Stool-sample 16S rRNA screener  |  gut-brain-axis microbes"); s.setObjectName("sub")
        lay.addWidget(t); lay.addWidget(s)

        row = QHBoxLayout()
        self.fl = QLabel("No file selected - drop a FASTA file here or click Browse"); self.fl.setWordWrap(True)
        b = QPushButton("Browse FASTA..."); b.clicked.connect(self.browse)
        row.addWidget(self.fl, 1); row.addWidget(b); lay.addLayout(row)

        row = QHBoxLayout()
        row.addWidget(QLabel("Min read length (bp):")); self.minlen = QSpinBox(); self.minlen.setRange(50, 1500); self.minlen.setValue(150); row.addWidget(self.minlen)
        row.addSpacing(16); row.addWidget(QLabel("Min reads to call 'Detected':")); self.minreads = QSpinBox(); self.minreads.setRange(1, 1000); self.minreads.setValue(3); row.addWidget(self.minreads)
        row.addStretch(1)
        self.go = QPushButton("ANALYSE"); self.go.setObjectName("go"); self.go.setEnabled(False); self.go.clicked.connect(self.start); row.addWidget(self.go)
        self.exp = QPushButton("Export CSV..."); self.exp.setEnabled(False); self.exp.clicked.connect(self.export); row.addWidget(self.exp)
        lay.addLayout(row)

        self.bar = QProgressBar(); self.bar.setVisible(False); lay.addWidget(self.bar)
        self.msg = QLabel(""); self.msg.setObjectName("sub"); lay.addWidget(self.msg)

        cards = QHBoxLayout(); self.c = {}
        for k, l in [("in", "SEQUENCES IN FILE"), ("kept", "KEPT AFTER CLEANING"), ("cls", "CLASSIFIED %"),
                     ("det", "TARGET MICROBES DETECTED"), ("oth", "OTHER (NON-TARGET) %"), ("unc", "UNCLASSIFIED %")]:
            self.c[k] = Card(l); cards.addWidget(self.c[k])
        lay.addLayout(cards)

        self.tabs = QTabWidget(); lay.addWidget(self.tabs, 1)
        self.tree = QTreeWidget(); self.tree.setAlternatingRowColors(True); self.tree.setRootIsDecorated(True)
        self.tree.setHeaderLabels(["Microbe", "Status", "Reads", "% of sample", "Mean identity", "Ref. seqs", "Role in gut-brain axis"])
        self.tree.header().setSectionResizeMode(QHeaderView.Interactive)
        for i, w in enumerate([250, 170, 70, 95, 105, 80, 600]): self.tree.setColumnWidth(i, w)
        self.tabs.addTab(self.tree, "Target microbes")
        self.bars = Bars(); sc = QScrollArea(); sc.setWidgetResizable(True); sc.setWidget(self.bars)
        self.tabs.addTab(sc, "Sample composition")
        self.other = QTreeWidget(); self.other.setHeaderLabels(["Other genus / group (not in target list)", "Reads", "% of sample"])
        self.other.setColumnWidth(0, 400); self.tabs.addTab(self.other, "Other bacteria")
        self.log = QTextEdit(); self.log.setReadOnly(True); self.tabs.addTab(self.log, "Cleaning log")
        ab = QTextEdit(); ab.setReadOnly(True); ab.setPlainText(ABOUT); ab.setFont(QFont("Consolas", 10)); self.tabs.addTab(ab, "How it works / limits")
        self.populate_empty()

    # ---------------------------------------------------------------- file
    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls(): e.acceptProposedAction()

    def dropEvent(self, e):
        u = e.mimeData().urls()
        if u: self.set_file(u[0].toLocalFile())

    def browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Choose stool-sample FASTA", "", "Sequence files (*.fasta *.fa *.fna *.fas *.fastq *.fq *.gz *.txt);;All files (*)")
        if p: self.set_file(p)

    def set_file(self, p):
        self.path = p; self.fl.setText(p); self.go.setEnabled(True)

    # ---------------------------------------------------------------- run
    def start(self):
        if not self.path: return
        self.go.setEnabled(False); self.exp.setEnabled(False); self.bar.setVisible(True); self.bar.setRange(0, 0)
        self.worker = Worker(self.path, CleanParams(min_len=self.minlen.value()), Params(min_reads=self.minreads.value()), self.db)
        self.worker.progress.connect(self.on_progress); self.worker.done.connect(self.on_done); self.worker.failed.connect(self.on_fail)
        self.worker.start()

    def on_progress(self, text, a, b):
        self.msg.setText(text)
        if b: self.bar.setRange(0, b); self.bar.setValue(a)
        else: self.bar.setRange(0, 0)

    def on_fail(self, e):
        self.bar.setVisible(False); self.go.setEnabled(True); self.msg.setText("")
        QMessageBox.warning(self, "GutScan", e)

    def on_done(self, res, st):
        self.bar.setVisible(False); self.go.setEnabled(True); self.exp.setEnabled(True)
        self.res, self.st = res, st; self.msg.setText("Done.")
        tot = res.n_total or 1
        self.c["in"].set(f"{st.reads_in:,}"); self.c["kept"].set(f"{st.reads_kept:,}")
        self.c["cls"].set(f"{100 * res.n_classified / tot:.1f}"); self.c["unc"].set(f"{100 * res.n_unclassified / tot:.1f}")
        self.c["oth"].set(f"{100 * res.n_other / tot:.1f}")
        self.c["det"].set(sum(1 for r in res.targets if r.status == "Detected" and r.target.rank != "species"))
        self.fill_tree(res); self.fill_bars(res); self.fill_other(res); self.fill_log(st)

    # ------------------------------------------------------------- tables
    def populate_empty(self):
        class R:  # placeholder rows before any analysis
            pass
        self.fill_tree(None)

    def fill_tree(self, res):
        self.tree.clear(); by = {r.target.name: r for r in res.targets} if res else {}
        for g in T.GROUPS:
            gi = QTreeWidgetItem([g]); f = gi.font(0); f.setBold(True); gi.setFont(0, f); gi.setForeground(0, QColor(ACCENT))
            self.tree.addTopLevelItem(gi); gi.setFirstColumnSpanned(True); gi.setExpanded(True); parents = {}
            for t in T.TARGETS:
                if t.group != g: continue
                r = by.get(t.name)
                if t.rank == "species":
                    label = f"{t.genus[0]}. {t.species}  (putative)"
                    item = QTreeWidgetItem(parents.get(t.parent, gi))
                else:
                    label = t.name + ("  (family)" if t.rank == "family" else "")
                    item = QTreeWidgetItem(gi); parents[t.name] = item
                item.setText(0, label); item.setText(6, t.role if t.rank != "species" else "")
                tip = t.note or t.role; item.setToolTip(0, tip); item.setToolTip(6, tip)
                if r:
                    status = r.status
                    if r.reads == 0 and r.ref_count == 0: status = "No reference data"
                    elif r.reads == 0 and t.rank == "species" and r.ref_count < 3: status = "Not detected (few refs)"
                    item.setText(1, status); item.setText(2, f"{r.reads:,}"); item.setText(3, f"{r.pct_sample:.3f} %")
                    item.setText(4, f"{r.mean_identity:.1f} %" if r.reads else "-"); item.setText(5, str(r.ref_count))
                    col = {"Detected": ACCENT, "Trace (below threshold)": "#e0883a"}.get(r.status, "#6b6f78")
                    item.setForeground(1, QColor(col))
                    if r.status == "Detected":
                        f = item.font(0); f.setBold(True); item.setFont(0, f)
                    else:
                        item.setForeground(0, QColor("#8a8e97"))
                for c in (2, 3, 4, 5): item.setTextAlignment(c, Qt.AlignRight | Qt.AlignVCenter)

    def fill_bars(self, res):
        tot = res.n_total or 1; rows = []
        pal = {T.SCFA: "#f2c200", T.NEURO: "#4fb0c6", T.MENTAL: "#9b8cff", T.HARM: "#e0583a"}
        for r in res.targets:
            if r.exclusive:
                nm = r.target.name + (" (other)" if r.target.rank != "species" and any(
                    x.target.parent == r.target.name and x.exclusive for x in res.targets) else "")
                rows.append((nm, 100 * r.exclusive / tot, pal[r.target.group]))
        rows.append(("Other bacteria (non-target)", 100 * res.n_other / tot, "#5b6270"))
        rows.append(("Unclassified", 100 * res.n_unclassified / tot, "#33363d"))
        rows.sort(key=lambda x: -x[1]); self.bars.set_rows([r for r in rows if r[1] > 0])

    def fill_other(self, res):
        self.other.clear(); tot = res.n_total or 1
        for g, n in res.other_genera.most_common(300):
            i = QTreeWidgetItem([g, f"{n:,}", f"{100 * n / tot:.3f} %"]); self.other.addTopLevelItem(i)
            i.setTextAlignment(1, Qt.AlignRight); i.setTextAlignment(2, Qt.AlignRight)

    def fill_log(self, st):
        L = [f"Sequences read from file: {st.reads_in:,} (weighted {st.weighted_in:,})",
             f"Kept after cleaning:      {st.reads_kept:,}  ({st.unique_kept:,} distinct sequences)", ""]
        if st.dropped:
            L.append("Dropped:")
            L += [f"   {n:>8,}  {why}" for why, n in st.dropped.most_common()]
        else:
            L.append("Nothing was dropped.")
        L += ["", "Cleaning steps: upper-case, U->T, remove gaps/digits, IUPAC codes -> N, trim N ends,",
              "drop too short / too long / many-N / low-complexity / homopolymer-artifact sequences, merge identical sequences",
              "(';size=N' counts in headers are honoured)."]
        self.log.setPlainText("\n".join(L))

    def export(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save results", "gutscan_results.csv", "CSV (*.csv)")
        if p: to_csv(self.res, self.st, p); self.msg.setText(f"Saved {p}")


def main():
    app = QApplication(sys.argv); app.setStyleSheet(STYLE)
    w = Main(); w.show(); sys.exit(app.exec())


if __name__ == "__main__":
    main()
