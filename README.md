# Kruskal-Wallis + Dunn's post-hoc tools

A Streamlit app and a command-line script that compare three or more groups using the Kruskal-Wallis test,
followed by Dunn's post-hoc test with a multiple-comparison correction.
No data is included: the scripts expect your own CSV file.

## Setup

Install Python 3, then in a terminal in this folder:

```
py -m pip install -r requirements.txt
```

## Data layout

One CSV with one row per measurement:

| Material | Specimen | Sa (um) |
|---|---|---|
| Benchmark | 1 | 63 |
| Benchmark | 1 | 30 |
| Material 1 | 1 | 100 |
| … | … | … |

- **Group column** (for example `Material`): the groups being compared.
- **Value column** (for example `Sa (um)`): the measurement.
- **Specimen column** (optional): when a specimen is measured more than once, for example on both matching
  fracture faces, the repeats can be averaged so each specimen counts once.

Other columns are ignored. Upload the file in the app's sidebar, or save it as `data.csv` next to the scripts.

## The app: `kruskal_dunn_app.py`

Double-click `Run Kruskal-Wallis app.bat`, or:

```
py -m streamlit run kruskal_dunn_app.py
```

- **Data:** choose the group, value and specimen columns, and pick which groups to include and in what order.
- **Test:** all pairwise comparisons or each group against a control. Holm, Bonferroni,
  Benjamini-Hochberg or no correction. Adjustable significance level.
- **Results:** H, degrees of freedom, p, effect size ε² = H / (N − 1), and a table of medians, quartiles and
  mean ranks. Dunn's z, raw and adjusted p for each comparison. Results download as CSV, and a text report
  can be copied into notes.
- **Plot:** box plot, points, or both, with significance brackets (* p < 0.05, ** p < 0.01, *** p < 0.001).
  Colours, labels, font size and figure size are adjustable. Figures save to the Desktop as PNG, SVG or PDF.

## The command-line script: `kruskal_dunn.py`

```
py kruskal_dunn.py data.csv
py kruskal_dunn.py data.csv --specimen Specimen --control Benchmark
```

| Option | Meaning |
|---|---|
| `--group`, `--value` | column names (guessed if left out) |
| `--specimen` | average the values that share a group and this column |
| `--control` | also compare every group with this one |
| `--correction` | `Holm` (default), `Bonferroni`, `Benjamini-Hochberg` or `None` |

## Method

- **Kruskal-Wallis:** `scipy.stats.kruskal`, which corrects for ties.
- **Dunn's test:** ranks are pooled across all included groups, and
  z = (R̄ᵢ − R̄ⱼ) / √[(N(N + 1)/12 − Σ(t³ − t)/(12(N − 1))) · (1/nᵢ + 1/nⱼ)], with a two-sided p from the
  normal distribution. Σ(t³ − t) is the tie correction.
- **Corrections:** these apply only to the comparisons that were run. Comparing each group with a control
  therefore gives smaller adjusted p values than comparing all pairs.

The statistics are in `kw_dunn.py`, which needs only numpy, pandas and scipy and can be imported by other scripts.
