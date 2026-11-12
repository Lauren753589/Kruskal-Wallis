"""
Kruskal-Wallis + Dunn's post-hoc app
------------------------------------
Run from a terminal with:
    py -m streamlit run kruskal_dunn_app.py

Expected CSV layout (one row per measurement):
    Group column    : e.g. Material
    Value column    : e.g. Sa (um)
    Specimen column : optional; values that share a group and specimen can be averaged first
"""

import io
import re
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import kw_dunn

HERE = Path(__file__).parent
DEFAULT_DATA = HERE / "data.csv"
DESKTOP = next((p for p in (Path.home() / "Desktop", Path.home() / "OneDrive" / "Desktop") if p.is_dir()), HERE)

st.set_page_config(page_title="Kruskal-Wallis + Dunn", layout="wide")
st.title("Kruskal-Wallis with Dunn's post-hoc test")


# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================
@st.cache_data(max_entries=10)
def read_csv(data):
    return pd.read_csv(io.BytesIO(data), encoding="utf-8-sig")


def unique_path(path):
    """Add ' (1)', ' (2)', ... so an existing file is never overwritten."""
    out, i = path, 1
    while out.exists():
        out = path.with_name(f"{path.stem} ({i}){path.suffix}")
        i += 1
    return out


def bracket_tiers(pairs, pos):
    """Stack significance brackets: narrow ones low, no two overlapping on the same tier."""
    tiers, placed = [], []
    for a, b, label in sorted(pairs, key=lambda t: abs(pos[t[0]] - pos[t[1]])):
        x0, x1 = sorted((pos[a], pos[b]))
        tier = next((i for i, used in enumerate(tiers) if all(x1 < u0 or x0 > u1 for u0, u1 in used)), None)
        if tier is None:
            tiers.append([])
            tier = len(tiers) - 1
        tiers[tier].append((x0, x1))
        placed.append((x0, x1, tier, label))
    return placed, len(tiers)


# =============================================================================
# 2. SIDEBAR: DATA
# =============================================================================
st.sidebar.header("1. Data")
upload = st.sidebar.file_uploader("Data CSV (otherwise data.csv in this folder is used)", type="csv")
if upload is not None:
    raw = read_csv(upload.getvalue())
elif DEFAULT_DATA.exists():
    raw = read_csv(DEFAULT_DATA.read_bytes())
else:
    st.info("Upload a CSV in the sidebar, or save one as `data.csv` next to this script. "
            "Use one row per measurement, for example:")
    st.table(pd.DataFrame({"Material": ["A", "A", "B", "B"], "Specimen": [1, 1, 1, 1],
                           "Value": [63, 30, 100, 115]}))
    st.stop()

cols = list(raw.columns)
g_guess, v_guess, s_guess = kw_dunn.guess_columns(raw)
if v_guess is None:
    st.error("No numeric column found in the CSV.")
    st.stop()

group_col = st.sidebar.selectbox("Group column", cols, index=cols.index(g_guess))
value_options = [c for c in cols if c != group_col]
value_col = st.sidebar.selectbox("Value column", value_options,
                                 index=value_options.index(v_guess) if v_guess in value_options else 0)
spec_options = ["(none: every row is one value)"] + [c for c in cols if c not in (group_col, value_col)]
specimen_col = st.sidebar.selectbox(
    "Average repeats per specimen", spec_options,
    index=spec_options.index(s_guess) if s_guess in spec_options else 0,
    help="Values that share a group and this column are averaged before testing, so each specimen counts once "
         "(for example the two matching fracture faces).")
specimen_col = None if specimen_col == spec_options[0] else specimen_col

all_groups = list(dict.fromkeys(raw[group_col].dropna().astype(str).str.strip()))
chosen = st.sidebar.multiselect("Groups (in plot order)", all_groups, default=all_groups)

# =============================================================================
# 3. SIDEBAR: TEST
# =============================================================================
st.sidebar.header("2. Test")
mode = st.sidebar.segmented_control("Comparisons", ["All pairs", "Each vs control"],
                                    default="All pairs", required=True)
control = st.sidebar.selectbox("Control group", chosen) if mode == "Each vs control" and chosen else None
correction = st.sidebar.selectbox(
    "Multiple-comparison correction", list(kw_dunn.CORRECTIONS),
    help="Holm and Bonferroni control the family-wise error rate (Holm is never less powerful). "
         "Benjamini-Hochberg controls the false discovery rate.")
alpha = st.sidebar.number_input("Significance level (α)", 0.001, 0.2, 0.05, step=0.01, format="%.3f")

# =============================================================================
# 4. SIDEBAR: PLOT STYLE
# =============================================================================
st.sidebar.header("3. Plot style")
plot_type = st.sidebar.segmented_control("Show", ["Box + points", "Box", "Points"],
                                         default="Box + points", required=True)
show_brackets = st.sidebar.toggle("Significance brackets", True,
                                  help="Drawn for comparisons with adjusted p below α.")
start_zero = st.sidebar.toggle("Start y-axis at 0", False)
with st.sidebar.expander("Labels and size"):
    x_title = st.text_input("X-axis title", group_col)
    y_title = st.text_input("Y-axis title", value_col)
    font_size = st.slider("Font size", 8, 32, 16)
    plot_width = st.number_input("Width (px)", 300, 3000, 900, step=50)
    plot_height = st.number_input("Height (px)", 300, 3000, 600, step=50)
with st.sidebar.expander("Group colours"):
    palette = px.colors.qualitative.D3  # hex codes, which st.color_picker requires
    colours = {g: st.color_picker(g, palette[i % len(palette)], key=f"colour_{g}") for i, g in enumerate(chosen)}

# =============================================================================
# 5. RUN THE TESTS
# =============================================================================
groups = kw_dunn.groups_from_table(raw, group_col, value_col, specimen_col, chosen)
if len(groups) < 2:
    st.warning("Pick at least two groups that have numeric values.")
    st.stop()
if control not in groups:
    control = None if mode == "All pairs" else next(iter(groups))

try:
    kw = kw_dunn.kruskal_wallis(groups)
except ValueError as e:  # e.g. every value identical
    st.error(f"The Kruskal-Wallis test could not be run: {e}")
    st.stop()

pairs = kw_dunn.comparisons(list(groups), control)
res = kw_dunn.dunn(groups, pairs, correction)
res["Significance"] = [kw_dunn.stars(p) if p < alpha else "ns" for p in res["p (adjusted)"]]
desc = kw_dunn.describe(groups)
basis = f"averaged per {specimen_col}" if specimen_col else "every row counted as one value"

# =============================================================================
# 6. TABS
# =============================================================================
tab_res, tab_plot, tab_data = st.tabs(["Results", "Plot", "Data"])

with tab_res:
    with st.container(horizontal=True):
        st.metric("Kruskal-Wallis H", f"{kw['H']:.2f}", border=True)
        st.metric("Degrees of freedom", kw["df"], border=True)
        st.metric("p", f"{kw['p']:.3g}", border=True)
        st.metric("Effect size ε²", f"{kw['epsilon2']:.2f}", border=True,
                  help="ε² = H / (N − 1). Roughly 0.01 small, 0.08 medium, 0.26 large.")
        st.metric("N", kw["N"], border=True)
    if kw["p"] < alpha:
        st.success(f"At least one group differs (p = {kw['p']:.3g} < α = {alpha}). Values are {basis}.")
    else:
        st.info(f"No significant difference between groups (p = {kw['p']:.3g} ≥ α = {alpha}). "
                f"Values are {basis}. The post-hoc results below are shown for completeness only.")

    st.subheader("Descriptive statistics")
    st.dataframe(desc, column_config={c: st.column_config.NumberColumn(format="%.4g")
                                      for c in ["Median", "Q1", "Q3", "Min", "Max", "Mean rank"]})

    st.subheader(f"Dunn's post-hoc test, {correction} correction")
    st.dataframe(res, hide_index=True, column_config={
        "Mean rank difference": st.column_config.NumberColumn(format="%.2f"),
        "z": st.column_config.NumberColumn(format="%.2f"),
        "p": st.column_config.NumberColumn(format="%.4g"),
        "p (adjusted)": st.column_config.NumberColumn(format="%.4g"),
    })

    with st.container(horizontal=True):
        st.download_button("Download descriptive statistics CSV", desc.to_csv(), "descriptive_statistics.csv",
                           icon=":material/download:")
        st.download_button("Download Dunn's test CSV", res.to_csv(index=False), "dunn_posthoc.csv",
                           icon=":material/download:")

    with st.expander("Text report (copy into notes)"):
        sets = [(f"each group vs {control}" if control else "all pairs", pairs)]
        st.code(kw_dunn.report(groups, sets, correction, f"{value_col} by {group_col}, {basis}"), language=None)

with tab_plot:
    names = list(groups)
    fig = go.Figure()
    for g in names:
        v = groups[g]
        box_kw = dict(y=v, x=[g] * len(v), name=g, marker_color=colours.get(g), showlegend=False,
                      boxpoints="all" if plot_type != "Box" else False, jitter=0.35, pointpos=0,
                      marker_size=8, line_width=2)
        if plot_type == "Points":  # hide the box but keep the jittered points, then mark the median
            box_kw.update(fillcolor="rgba(0,0,0,0)", line_width=0)
        fig.add_trace(go.Box(**box_kw))
        if plot_type == "Points":
            fig.add_trace(go.Scatter(x=[g], y=[desc.loc[g, "Median"]], mode="markers", showlegend=False,
                                     marker=dict(symbol="line-ew-open", size=40, color="black",
                                                 line_width=2), hoverinfo="skip"))

    lo = min(v.min() for v in groups.values())
    hi = max(v.max() for v in groups.values())
    span = (hi - (0 if start_zero else lo)) or 1
    step = span * 0.09
    top = hi + span * 0.05

    sig = res[res["p (adjusted)"] < alpha]
    if show_brackets and len(sig):
        pos = {g: i for i, g in enumerate(names)}
        placed, n_tiers = bracket_tiers(
            [(r["Group A"], r["Group B"], kw_dunn.stars(r["p (adjusted)"])) for _, r in sig.iterrows()], pos)
        line = dict(color="black", width=1.5)
        for x0, x1, tier, label in placed:
            y = hi + step * (tier + 1)
            fig.add_shape(type="line", x0=x0, x1=x1, y0=y, y1=y, line=line, opacity=1)
            for x in (x0, x1):
                fig.add_shape(type="line", x0=x, x1=x, y0=y - step * 0.3, y1=y, line=line, opacity=1)
            fig.add_annotation(x=(x0 + x1) / 2, y=y, text=label, showarrow=False, yanchor="bottom",
                               font_size=font_size)
        top = hi + step * (n_tiers + 1)

    fig.update_layout(template="simple_white", width=plot_width, height=plot_height,
                      font=dict(size=font_size, color="black"), margin=dict(t=30),
                      xaxis=dict(title=x_title, categoryorder="array", categoryarray=names),
                      yaxis=dict(title=y_title, range=[0 if start_zero else lo - span * 0.05, top]))
    st.plotly_chart(fig, width="content")
    st.caption(f"Brackets: {correction}-adjusted Dunn's test. * p < 0.05, ** p < 0.01, *** p < 0.001. "
               f"Values {basis}.")

    c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
    file_name = c1.text_input("File name", f"{value_col} by {group_col}")
    fmt = c2.selectbox("Format", ["png", "svg", "pdf"])
    if c3.button("Save to Desktop", type="primary", icon=":material/save:"):
        try:
            safe_name = re.sub(r'[<>:"/\\|?*]', "_", file_name).strip() or "kruskal_dunn"  # no illegal characters
            out = unique_path(DESKTOP / f"{safe_name}.{fmt}")
            with st.spinner("Saving…"):
                fig.write_image(out, width=plot_width, height=plot_height, scale=3 if fmt == "png" else 1)
            st.success(f"Saved to {out}")
        except Exception as e:  # e.g. the file is open in another program
            st.error(f"Could not save the figure: {e}")

with tab_data:
    st.subheader("Values used in the test")
    used = pd.DataFrame([(g, x) for g, v in groups.items() for x in v], columns=[group_col, value_col])
    st.dataframe(used, hide_index=True)
    st.download_button("Download values CSV", used.to_csv(index=False), "values_used.csv",
                       icon=":material/download:")
    st.subheader("Uploaded data")
    st.dataframe(raw, hide_index=True)
