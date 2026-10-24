# Code ↔ Report Traceability Map

> **Status.** This map documents the *earlier* draft of the report. The report has
> since been rewritten (`main.tex`) around the null-model ladder, and the analysis
> now lives in `protein_go_geometry.ipynb`, whose figures are
> `figures/fig_*.png`. The discrepancy notes below remain valid and were addressed
> in the rewrite.

Companion document for the report *Geometric Structure of GO Labels in Protein
Embeddings: Linear Probes, Spectral Concentration, and a Preliminary Analysis of
Superposition* (October 2026).

| Item | Location |
| --- | --- |
| Code | `Protein.ipynb` (91 cells, 11 numbered regions) |
| Data | `protein_dataset/` (see `protein_dataset/README.md`) |
| Figures referenced by the report | `figures/*.png` |
| Small 5-protein fixtures | `swissprot_5_proteins.csv`, `swissprot_5_embeddings.npy`, `protein_embedding.npy` |

Cell numbers below are 0-based positions in the current `Protein.ipynb`. Each entry
also names the region (`region 6.7` etc.) so the reference survives cells being inserted.

---

## 1. Report section → notebook region

| Report section | Notebook region | What the code does there |
| --- | --- | --- |
| §1 Motivation, Q1–Q4 | §6.4, §6.5 (cells 62, 64) | Numerical rank of `W`; rank of random row subsets for M = 2 … 75 |
| §1 Motivation (approximate structure) | §6.1–§6.7 (cells 56–68) | SVD, cumulative energy, effective rank, isotropic nulls |
| §2.1 Data (Table 1) | §2.1–§2.5, §3, §4.4 (cells 21–31, 40) | Filtering, label counts, 75-label selection, `Y` construction, artefact save |
| §2.1 "recorded data files" | `protein_dataset/` + root fixtures | See Section 5.1 — the list in the report is out of date |
| §2.2 Probes (Eq. 1–2), metrics | §4.1, §4.2, §5.2–§5.4 (cells 34, 36, 45, 47, 49) | ESM-2 load, mean pooling, `OneVsRestClassifier`, ROC-AUC / AP, `W` extraction |
| §3.1 Q1: collinearity (Eq. 3–5) | §6.4, §6.5, §7.1 (cells 62, 64, 71) | Rank test + cosine-similarity matrix |
| §3.2 Q2: coplanarity (Eq. 6) | §6.5 (cell 64) | Rank of 3-row subsets |
| §3.3 Q3: M-vector dependence (Eq. 7) | §6.4, §6.5 (cells 62, 64) | Full-rank check and subset ranks |
| §3.4 Q4: spectral concentration (Eq. 8–11) | §6.1–§6.3 (cells 56, 58, 60) | SVD, `E(k)`, `p_i`, effective rank |
| §3.5 Label correlation (Eq. 12, residualization) | §7.1–§7.4 (cells 71, 73, 75, 77) | Co-occurrence, label correlation, Spearman, residualized rank correlation |
| §3.6 Random-direction baseline | §6.6, §6.7 (cells 66, 68) | Row-norm-matched isotropic nulls, 500 and 1,000 draws |
| §4.1 Pairwise results, Fig. 1 | §7.1, §7.2 (cells 71, 73) | Pair table and scatter |
| §4.2 Spectrum results, Figs. 2–3 | §6.1, §6.2 (cells 56, 58) | Spectrum and cumulative-energy prints/plots |
| §4.3 Table 2 (51.96 / 60.91) | §6.7 (cell 68) | 1,000-draw isotropic null |
| §4.3 second batch (57.746 / 49.645) | §8.2, §8.3 (cells 82, 84) | 200 label-shuffle permutations, row-norm-matched variant |
| §4.3 "another label-shuffling experiment" | §8.2 (cell 82, raw variant) | Unmatched label-shuffle results |
| §4.4 Table 3 (subset ranks) | §6.4, §6.5 (cells 62, 64) | `rank(W) = 75` and subset checks |
| §5 Discussion, §7 Conclusion | §10 (cell 90) | Written summary of the same findings |
| §6 experiments 1–6 | not implemented | Proposed future work; none of it exists in the notebook yet |

## 2. Equations → code

| Report item | Code |
| --- | --- |
| Eq. (1) logistic probe `σ(w·x + b)` | cell 45 (`LogisticRegression`), numeric check in cell 53 |
| Eq. (2) stacked `W ∈ R^{75×480}` | cell 49 (`np.vstack([clf.coef_ ...])`) |
| Eq. (3) cosine similarity | cell 60 (`cosine_sim`), cell 71 (`weight_cosine`) |
| Eq. (4) `|c_ij| = 1` | cell 64 (M = 2 rank check); `cosine_sim` available in cell 60 |
| Eq. (5) `rank([w_i; w_j]) < 2` | cell 64, `M=2` row |
| Eq. (6) `rank([w_i; w_j; w_k]) < 3` | cell 64, `M=3` row |
| Eq. (7) `rank(V) < M` | cell 64 (loop over M) |
| Eq. (8) SVD `W = UΣVᵀ` | cells 56, 60, 62 (`np.linalg.svd`) |
| Eq. (9) `E(k)` | cell 56 (`cumulative_energy`), plotted in cell 58 |
| Eq. (10) `p_i = σ_i² / Σσ_j²` | cells 60, 68 |
| Eq. (11) effective rank | cells 60, 68, 80 (`spectral_metrics`) |
| Eq. (12) `C(75,2) = 2775` | cell 71 (`pairs` loop), asserted implicitly by cell 73 (`Number of label pairs: 2775`) |
| §3.5 residualized rank correlation | cell 75 (`rankdata` + `lstsq` + `spearmanr` → 0.4872); cell 77 (same with `pearsonr` → 0.4889) |
| §3.6 row-norm-matched random directions | cell 66 (500 draws), cell 68 (1,000 draws) |
| Eq. (13) `S = ‖Q₁ᵀQ₂‖_F² / k` (future work) | not implemented |

## 3. Numbers in the report → stored output

Every number quoted in the report resolves to one of the cells below. "Stored"
means the value appears verbatim in that cell's saved output.

| Report value | Source cell | Stored output line |
| --- | --- | --- |
| 9,895 protein records | 23, 31 | `Total records: 9895` / `Y shape: (9895, 75)` |
| 480 embedding dim | 36, 40 | `X shape: (9895, 480)` |
| 75 labels | 27, 31 | `Selected labels: 75` |
| 36,188 positives | 31 | `Positive labels: 36188` |
| 3.6572 labels/protein | 31 | `Labels per protein (mean): 3.6572006` |
| 9,786 proteins with ≥1 label | 31 | `Proteins with at least one selected label: 9786` |
| ROC-AUC 0.7841 / AP 0.2226 | 47 | `Mean ROC-AUC: 0.7841453256012005`, `Mean Average Precision: 0.22261515298787826` |
| 2,775 label pairs | 73 | `Number of label pairs: 2775` |
| Spearman ρ = 0.5074, p = 1.701e-181 | 73 | `Spearman correlation: 0.5074`, `P-value: 1.701e-181` |
| residualized ρ = 0.4889 (p = 1.010e-166) | 77 | `Partial Spearman correlation: 0.4889`, `Nominal p-value: 1.010e-166` |
| earlier residualized ρ = 0.4872 | 75 | `Partial rank correlation: 0.4872280648593841` |
| example pair 0.4177 / 0.4057 / 332 | 71 | first row of the "Most similar weight-vector pairs" table |
| σ₁ ≈ 41.96 | 56 | `First 20 singular values: [41.9558 31.6301 …]` |
| 7.63 / 34.38 / 54.68 / 70.19 / 82.13 / 90.91 / 97.08 % | 56 | `Cumulative energy:` block |
| Rank90 = 49 | 56, 58, 68 | `Top 50 directions: 90.9058%`, `90% energy: 49 directions`, `Observed rank90: 49` |
| effective rank 54.54 | 60, 68 | `Effective rank: 54.54` |
| null mean 51.959, interval [51, 52] | 68 | `Null rank90 mean: 51.959`, `[51. 52.]` |
| null mean 60.906, interval [60.665, 61.142] | 68 | `Null effective rank mean: 60.9057…`, `[60.6648… 61.1416…]` |
| numerical rank 75 | 62 | `Numerical rank of W: 75` |
| subset ranks 2…75, 0/100 deficient | 64 | `M= 2: mean rank=2.00, rank-deficient subsets=0/100` … |
| 57.746 / [57.360, 58.158] / p = 0.005 | 82–84 | reproduced, **not stored** (see Section 5.2) |
| 49.645 / [49, 50] / p = 0.3582 | 82–84 | reproduced, **not stored** (see Section 5.2) |

## 4. Figures

| Report figure | File | Produced by | Notes |
| --- | --- | --- | --- |
| Fig. 1 weight similarity vs label correlation | `figures/weight_similarity_vs_label_correlation.png` | cell 73 | 1750×1250 px, 250 dpi |
| Fig. 2 singular-value spectrum | `figures/singular_value_spectrum.png` | cell 56 (first plot) | 2000×1250 px, 250 dpi |
| Fig. 3 cumulative spectral energy | `figures/cumulative_spectral_energy.png` | cell 56 (second plot) | 2000×1250 px, 250 dpi |

The notebook contains two further plots that the report does not use: the
alternative cumulative-energy curve in region 6.2 (cell 58, with the 80/90/95/99%
annotation) and the cosine-similarity heatmap in region 6.3 (cell 60). The figures above
were regenerated from the same `W` and plotting code at 250 dpi so that they are
print-ready; the versions stored inside the notebook are smaller inline renders
(about 690 × 470 px), which would look soft at `0.88\textwidth` in print.

## 5. Discrepancies between report and code

### 5.1 The data-file list in §2.1 does not match the disk

The report lists `swissprot_5_embedings.npy`, `swissprot_5_proteins.csv`,
`selected_go_labels.csv`, `go_cc_label_counts.csv` and `esm2_embedding.npy`, and
states that the misspelling "embedings" is preserved from an existing filename.
There is no file with that spelling; the actual artefacts are:

| Report says | Actual file | Shape |
| --- | --- | --- |
| `esm2_embedding.npy` | `protein_dataset/esm2_embeddings.npy` | 9,895 × 480 |
| `swissprot_5_embedings.npy` | `swissprot_5_embeddings.npy` | 5 × 480 (fixture) |
| — | `swissprot_5_proteins.csv` | 5 × 3 (fixture) |
| `selected_go_labels.csv` | `protein_dataset/selected_go_labels.csv` | 75 × 1 |
| `go_cc_label_counts.csv` | `protein_dataset/go_cc_label_counts.csv` | 1,779 × 2 |
| — | `protein_dataset/go_labels.npy` | 9,895 × 75 |
| — | `protein_dataset/protein_ids.csv` | 9,895 × 1 |
| — | `protein_dataset/uniprot_swissprot_go_cc.csv` | 9,895 × 3 |

The five-protein fixtures come from the earlier exploratory notebook
`Untitled1.ipynb` (cell 4); `protein_embedding.npy` is a single 480-dimensional
vector saved by `Untitled.ipynb` (cell 5). None of them is produced by
`Protein.ipynb`. The sentence about preserving the spelling should be deleted and
the list replaced by the table above.

### 5.2 The 200-permutation numbers have no stored output

Section 4.3 of the report cites a row-norm-matched label-shuffle null with
effective-rank mean 57.746, interval [57.360, 58.158], p = 0.005, and Rank90 mean
49.645, interval [49, 50], p = 0.3582. Those values do not appear anywhere in
`Protein.ipynb`, in `Untitled*.ipynb`, or in the checkpoint copies: cell 82 stores
only the `tqdm` bar and "Permutation experiment complete.", and cell 84 (the
`summarize_null` call) has no output at all.

The experiment itself is cells 80 → 82 → 84. Re-running it reproduces the report's
numbers exactly (see Section 6 below), so the values are correct; they simply have no
saved evidence in the notebook. Either re-run cells 82 and 84 and keep the outputs,
or drop the numbers from the report.

### 5.3 "Row-norm-matched" names three different nulls

| Report wording | Actual null | Cell |
| --- | --- | --- |
| "row-norm-matched random-direction baseline" (Table 2, 51.96 / 60.91) | Gaussian directions, matched row norms, **labels untouched** | 66, 68 |
| "200 permutations … row-norm-matched condition" (57.746) | **labels shuffled per column**, then row norms matched | 82, 84 |
| "another label-shuffling experiment without row-norm matching" | same shuffle, raw weights | 82, 84 |

All three keep row norms, but they randomize different objects (directions vs.
labels vs. both). The report should name them by what is randomized, e.g.
"direction-randomization null", "label-shuffle null (norm-matched)", "label-shuffle
null (raw)".

### 5.4 A fourth null exists but is not in the report

Cell 86 runs a **row-permutation** null (shuffle the rows of `Y_train` as a block,
preserving label co-occurrence, 100 refits, row-norm matched). Its stored output is:

```
Effective rank - row permutation
Observed: 54.540   Null mean: 54.322   95% interval: [53.567, 54.983]   p = 0.7327
Rank90 - row permutation
Observed: 49.000   Null mean: 47.980   95% interval: [48, 48]           p = 1.0000
```

This is the null that most directly addresses the report's own question, because it
preserves the label co-occurrence structure. It says the observed concentration is
indistinguishable from what label structure alone produces — a result the report
currently does not discuss, although it is the strongest evidence in the notebook
against a superposition interpretation (see Section 7).

### 5.5 Smaller points

- Table 2's 95% interval for Rank90 is `[51, 52]`; Rank90 is an integer, so the
  interval is discrete. Fine, but worth a footnote.
- §4.3 computes a 10.5% relative reduction and a 5.7% Rank90 reduction; both are
  correct against the stored values (6.366/60.906 = 10.45%; 3/52 = 5.77%).
- §4.4 states that full row rank implies no dependent subset. Correct in exact
  arithmetic. The tolerance used in cell 62 is
  `max(75, 480) * eps * σ₁ ≈ 480 × 2.22e-16 × 41.96 ≈ 4.5e-12`; the report should
  quote it (it currently only says "the tolerance should be reported").
- The report describes the notebook as using "existing analysis results" rather
  than regenerating data, but `Protein.ipynb` does regenerate everything
  (downloads in §1.5, embeddings in §4.2). Worth aligning the wording.

### 5.6 Q1 could quote a margin instead of only invoking full rank

Two quantities that the report asserts verbally are easy to compute from the same
`W` and are currently absent from both the notebook and the paper:

- The closest pair is `GO:0005789` / `GO:0005783` with `|c| = 0.4177`, i.e. an
  angle of about 65°. For two unit rows, `σ₂/σ₁ = sqrt((1-|c|)/(1+|c|)) = 0.641`
  for that pair, so no pair comes close to `|c| = 1`.
- Among the 1,218 pairs whose label correlation is `≤ 0`, the largest cosine is
  0.192, which supports the report's sentence that some pairs stay positively
  aligned despite non-positive label correlation (the exact value is missing from
  the text).

Adding these two numbers would make the negative answer to Q1 quantitative rather
than an inference from `rank(W) = 75`.

## 6. Re-running the unfinished experiment

Cells 80, 82 and 84 were re-executed from the saved embeddings on 24 October 2026
(CPU, 631 s total for 200 label-shuffle and 100 row-permutation refits). The code,
seed (`42`), iteration count (`200`) and pipeline are unchanged, and the run
reproduces the report's §4.3 values to four decimals:

| Quantity | Report | Re-run |
| --- | --- | --- |
| mean ROC-AUC / mean AP | 0.7841 / 0.2226 | 0.784145 / 0.222615 |
| label shuffle, norm-matched, ER mean | 57.746 | 57.7462 |
| label shuffle, norm-matched, ER 95% interval | [57.360, 58.158] | [57.3599, 58.1589] |
| label shuffle, norm-matched, ER p | 0.005 | 0.0050 |
| label shuffle, norm-matched, Rank90 mean | 49.645 | 49.6450 |
| label shuffle, norm-matched, Rank90 p | 0.3582 | 0.3582 |
| row permutation, ER mean / p | 54.32 / 0.73 | 54.3218 / 0.7327 |
| row permutation, Rank90 mean / p | 47.98 / 1.00 | 47.98 / 1.0 |

A separate check of cells 66 and 68, which only need `W`, reproduces the stored
isotropic baselines as well (Rank90 null mean 51.962 with 500 draws and 51.959 with
1,000 draws; effective-rank null mean 60.906, interval [60.665, 61.142]). So every
number in Table 2 and in region 6.6/6.7 is reproducible from the saved data.

The only problem is bookkeeping: cells 82 and 84 were never saved with their
output, so a reader of the notebook sees no evidence for §4.3. Re-running the two
cells and saving the notebook fixes that.

## 7. What the four nulls actually say

With the values now verified, the report's remark that the results are "not fully
consistent" can be replaced by a much sharper statement. The four nulls differ in
exactly one thing — how much of the label structure they preserve — and, once row
norms are matched, they line up monotonically:

| Null | Randomized | Preserved | Null effective rank | Null Rank90 |
| --- | --- | --- | --- | --- |
| Observed `W` | — | — | **54.540** | **49** |
| Isotropic directions (region 6.7, cell 68) | weight directions | row norms | 60.906 [60.665, 61.142] | 51.96 [51, 52] |
| Label shuffle, raw (region 8.2, cell 82) | labels within column | positives per label | 52.508 [52.037, 53.026] | 44.24 [44, 45] |
| Label shuffle, norm-matched (region 8.2, cell 82) | labels within column | positives + row norms | 57.746 [57.360, 58.159] | 49.645 [49, 50] |
| Row permutation (region 8.4, cell 86) | protein–label pairing | label co-occurrence + row norms | 54.322 [53.567, 54.983] | 47.980 [48, 48] |

Read as a ladder:

1. The observed matrix is more concentrated than random *directions* (54.54 vs
   60.91), which is what §4.3 reports.
2. It is also more concentrated than weights learned from *shuffled labels* once
   the row norms are matched (54.54 vs 57.75, p = 0.005), i.e. it carries structure
   that independent per-label shuffling destroys.
3. But it is **indistinguishable** from weights fitted to the *same co-occurrence
   structure with the protein pairing broken* (54.54 vs 54.32, p = 0.73).

Steps 2 and 3 are consistent: shuffling labels removes their co-occurrence structure
and the null becomes less concentrated than the observation, while keeping the
co-occurrence structure and breaking only the protein pairing reproduces the
observation exactly. Under this reading the concentration tracks the co-occurrence
structure of the labels rather than protein-specific information in the embeddings,
and it is the strongest argument in the workspace against treating the spectrum as
direct evidence of superposition. The report currently omits this null entirely.

A caveat on top of that reading: the raw label-shuffle null (52.51) is *more*
concentrated than the observed matrix, so whether a shuffle null lands above or
below the observation depends on the row-norm matching step. That is exactly why
every null quoted in the report needs its matching rule stated (see Section 5.3).

*Values in this table: isotropic from cell 68 and row permutation from cell 86, both
stored; the two label-shuffle rows come from the re-run in Section 6 because cells
82/84 were saved without output.*

## 8. Suggested additions to the notebook

1. Save `W` once (`np.save("protein_dataset/probe_weights.npy", W)`) so figures and
   nulls can be recomputed without refitting.
2. Persist the output of cell 84 — this is the only result in the report that
   currently has no evidence on disk.
3. Merge `spectral_profile` (region 6.7) and `spectral_metrics` (region 8.1); they
   compute the same quantities under two names, which invites inconsistency.
4. Add a direct Q1/Q2 check (smallest `|c_ij|` distance from 1, smallest
   `σ₂/σ₁` over all pairs and triples) instead of inferring from full rank.
5. Report `σ_min/σ_max` and the rank tolerance next to the rank result.
6. Add the row-permutation null (region 8.4) to the report — it changes the conclusion.
