"""
Kruskal-Wallis test with Dunn's post-hoc test, from the command line.

    py kruskal_dunn.py data.csv
    py kruskal_dunn.py data.csv --specimen Specimen --control Benchmark

The group and value columns are guessed unless given with --group and --value.
With --specimen, values that share a group and specimen are averaged first.
With --control, each group is compared with the control as well as all pairs.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

import kw_dunn

parser = argparse.ArgumentParser(description="Kruskal-Wallis test with Dunn's post-hoc test.")
parser.add_argument("csv", nargs="?", default=Path(__file__).parent / "data.csv", help="data file (default: data.csv)")
parser.add_argument("--group", help="column holding the group names")
parser.add_argument("--value", help="column holding the measurements")
parser.add_argument("--specimen", help="average values that share a group and this column")
parser.add_argument("--control", help="also compare every group with this one")
parser.add_argument("--correction", default="Holm", choices=list(kw_dunn.CORRECTIONS))
args = parser.parse_args()

df = pd.read_csv(args.csv, encoding="utf-8-sig")
g_guess, v_guess, _ = kw_dunn.guess_columns(df)
group_col, value_col = args.group or g_guess, args.value or v_guess
for col in (group_col, value_col, args.specimen):
    if col and col not in df.columns:
        sys.exit(f"Column '{col}' not found. Columns are: {', '.join(map(str, df.columns))}")

groups = kw_dunn.groups_from_table(df, group_col, value_col, args.specimen)
if args.control and args.control not in groups:
    sys.exit(f"Control group '{args.control}' not found. Groups are: {', '.join(groups)}")

names = list(groups)
sets = []
if args.control:
    sets.append((f"each group vs {args.control}", kw_dunn.comparisons(names, args.control)))
sets.append(("all pairs", kw_dunn.comparisons(names)))

title = f"{value_col} by {group_col}" + (f", averaged per {args.specimen}" if args.specimen else "")
print(kw_dunn.report(groups, sets, args.correction, title))
