#!/usr/bin/env python3
"""
K. pneumoniae Phylogenomics Dashboard
Run:  python3 dashboard/app.py
Open: http://localhost:8050
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from Bio import Phylo
from pathlib import Path
import dash
from dash import dcc, html, Input, Output, dash_table

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).parent.parent / "results"
SUMMARY_F = BASE / "summary" / "summary.tsv"
TREEFILE  = BASE / "iqtree"  / "iqtree.treefile"
SNP_FILE  = BASE / "snippy"  / "pairwise_snp_matrix.tsv"

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(SUMMARY_F, sep="\t", dtype=str)

for c in ["contigs", "total_length_bp", "N50", "L50", "Ns_per_100kbp"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

QUAST_COLS = [c for c in ["contigs","total_length_bp","N50","L50","Ns_per_100kbp"] if c in df.columns]
KLEB_COLS  = [c for c in ["species","ST","Ybt","Clb","Iuc","Iro","RmpADC","rmpA2"] if c in df.columns]
VIR_COLS   = [c for c in ["Ybt","Clb","Iuc","Iro","RmpADC","rmpA2"] if c in df.columns]
META_COLS  = set(["sample_id"] + QUAST_COLS + KLEB_COLS + ["rmpA","rmpC","rmpD"])
AMR_COLS   = [c for c in df.columns if c not in META_COLS]

amr_bin = df[["sample_id"] + AMR_COLS].copy()
for c in AMR_COLS:
    amr_bin[c] = amr_bin[c].notna().astype(int)
amr_present = [c for c in AMR_COLS if amr_bin[c].sum() >= 2]

snp        = pd.read_csv(SNP_FILE, sep="\t", index_col=0)
snp_samps  = snp.index.tolist()
snp_raw    = snp.values.astype(float)
snp_disp   = snp_raw.copy()
np.fill_diagonal(snp_disp, np.nan)

palette  = px.colors.qualitative.Set3 + px.colors.qualitative.Pastel
sts      = sorted(df["ST"].dropna().unique()) if "ST" in df.columns else []
st_color = {st: palette[i % len(palette)] for i, st in enumerate(sts)}
st_color["Unknown"] = "#cccccc"


# ── Tree ──────────────────────────────────────────────────────────────────────
def build_tree_fig():
    tree = Phylo.read(str(TREEFILE), "newick")
    tree.ladderize()
    terminals = tree.get_terminals()
    tip_y = {t.name: i for i, t in enumerate(terminals)}
    xp, yp = {}, {}

    def calc(clade, x=0.0):
        xl = x + (clade.branch_length or 0.0)
        if clade.is_terminal():
            xp[id(clade)] = xl; yp[id(clade)] = tip_y[clade.name]
        else:
            for ch in clade.clades:
                calc(ch, xl)
            ys = [yp[id(c)] for c in clade.clades]
            xp[id(clade)] = xl; yp[id(clade)] = (min(ys)+max(ys))/2

    calc(tree.root)
    lx, ly = [], []

    def branches(clade):
        x = xp[id(clade)]
        if not clade.is_terminal():
            cys = [yp[id(c)] for c in clade.clades]
            lx.extend([x, x, None]); ly.extend([min(cys), max(cys), None])
            for ch in clade.clades:
                lx.extend([x, xp[id(ch)], None])
                ly.extend([yp[id(ch)], yp[id(ch)], None])
                branches(ch)

    branches(tree.root)

    name2st = {}
    if "ST" in df.columns:
        for _, row in df.iterrows():
            name2st[row["sample_id"]] = str(row.get("ST") or "Unknown")

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lx, y=ly, mode="lines",
                             line=dict(color="#444", width=1),
                             hoverinfo="none", showlegend=False))

    groups = {}
    for tip in terminals:
        groups.setdefault(name2st.get(tip.name, "Unknown"), []).append(tip)

    for st, tips in sorted(groups.items()):
        fig.add_trace(go.Scatter(
            x=[xp[id(t)] for t in tips], y=[yp[id(t)] for t in tips],
            mode="markers+text", name=st,
            text=[t.name for t in tips], textposition="middle right",
            textfont=dict(size=8),
            marker=dict(size=7, color=st_color.get(st, "#999")),
            hovertemplate="<b>%{text}</b><extra>"+st+"</extra>"
        ))

    fig.update_layout(
        height=max(700, len(terminals)*16),
        xaxis=dict(title="Branch length", showgrid=True, gridcolor="#eee"),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        legend=dict(title="Sequence Type"),
        margin=dict(l=10, r=230, t=30, b=40),
        plot_bgcolor="white", paper_bgcolor="white"
    )
    return fig


# ── Layout ────────────────────────────────────────────────────────────────────
app = dash.Dash(__name__, title="Kp Phylogenomics")

app.layout = html.Div([
    html.H2("K. pneumoniae Phylogenomics Dashboard",
            style={"textAlign":"center","marginBottom":"6px","color":"#2c3e50"}),
    html.P(f"{len(df)} genomes  ·  Core-genome SNP phylogeny  ·  IQ-TREE 2 GTR+G",
           style={"textAlign":"center","color":"#7f8c8d","marginBottom":"18px"}),
    dcc.Tabs(id="tabs", value="tab-tree", children=[
        dcc.Tab(label="🌳 Phylogeny",      value="tab-tree"),
        dcc.Tab(label="💊 AMR Heatmap",    value="tab-amr"),
        dcc.Tab(label="🧬 Sequence Types", value="tab-st"),
        dcc.Tab(label="📊 Assembly QC",    value="tab-qc"),
        dcc.Tab(label="🦠 Virulence",      value="tab-vir"),
        dcc.Tab(label="🔢 SNP Matrix",     value="tab-snp"),
    ]),
    html.Div(id="content", style={"padding":"20px"})
], style={"maxWidth":"1400px","margin":"auto","fontFamily":"Arial, sans-serif"})


# ── Callback ──────────────────────────────────────────────────────────────────
@app.callback(Output("content","children"), Input("tabs","value"))
def render(tab):

    if tab == "tab-tree":
        return html.Div([
            html.H3("Maximum-likelihood phylogeny (IQ-TREE 2, GTR+G, 1000 bootstraps)"),
            html.P("Tips coloured by MLST sequence type (Kleborate)."),
            dcc.Graph(figure=build_tree_fig())
        ])

    elif tab == "tab-amr":
        mat = amr_bin[["sample_id"]+amr_present].set_index("sample_id")
        mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
        fig = go.Figure(go.Heatmap(
            z=mat.values, x=mat.columns.tolist(), y=mat.index.tolist(),
            colorscale=[[0,"#f8f9fa"],[1,"#c0392b"]],
            showscale=False, xgap=1, ygap=1,
            hovertemplate="<b>%{y}</b><br>%{x}<extra></extra>"
        ))
        fig.update_layout(
            height=max(500, len(mat)*14+150),
            margin=dict(l=160,r=20,t=40,b=130),
            xaxis=dict(tickangle=-65, tickfont=dict(size=8)),
            yaxis=dict(tickfont=dict(size=9)),
            plot_bgcolor="white"
        )
        return html.Div([
            html.H3("AMR gene presence / absence (AMRFinderPlus)"),
            html.P(f"{len(amr_present)} genes present in ≥2 isolates. Red = present."),
            dcc.Graph(figure=fig)
        ])

    elif tab == "tab-st":
        if "ST" not in df.columns:
            return html.P("ST data not available.")
        counts = df["ST"].fillna("Unknown").astype(str).value_counts()
        counts = pd.DataFrame({"ST": counts.index, "count": counts.values})
        counts["label"] = counts["ST"]
        fig = px.bar(counts, x="label", y="count", color="label",
                     color_discrete_sequence=palette, text="count",
                     labels={"label":"Sequence Type","count":"Isolates"})
        fig.update_traces(textposition="outside")
        fig.update_layout(showlegend=False, height=450, plot_bgcolor="white",
                          xaxis=dict(categoryorder="total descending"))
        children = [html.H3("MLST Sequence Type distribution"), dcc.Graph(figure=fig)]
        if "species" in df.columns:
            sp = df["species"].fillna("Unknown").value_counts()
            sp = pd.DataFrame({"species": sp.index, "count": sp.values})
            sp_fig = px.pie(sp, values="count", names="species",
                            color_discrete_sequence=px.colors.qualitative.Pastel)
            sp_fig.update_layout(title="Species distribution")
            children.append(dcc.Graph(figure=sp_fig))
        return html.Div(children)

    elif tab == "tab-qc":
        fig1 = px.scatter(df, x="total_length_bp", y="N50",
                          color="contigs", size="contigs", hover_name="sample_id",
                          color_continuous_scale="RdYlGn_r",
                          title="N50 vs total assembly length",
                          labels={"total_length_bp":"Total length (bp)","contigs":"# contigs"})
        fig1.update_layout(height=450, plot_bgcolor="white")
        fig2 = px.histogram(df, x="contigs", nbins=20,
                            title="Contig count distribution",
                            color_discrete_sequence=["#3498db"],
                            labels={"contigs":"# contigs"})
        fig2.update_layout(height=350, plot_bgcolor="white")
        tbl = df[["sample_id"]+QUAST_COLS].round(0)
        return html.Div([
            html.H3("Assembly QC (QUAST)"),
            dcc.Graph(figure=fig1), dcc.Graph(figure=fig2),
            html.H4("Per-sample metrics"),
            dash_table.DataTable(
                data=tbl.to_dict("records"),
                columns=[{"name":c,"id":c} for c in tbl.columns],
                sort_action="native", page_size=20,
                style_table={"overflowX":"auto"}, style_cell={"fontSize":12},
                style_header={"fontWeight":"bold","backgroundColor":"#ecf0f1"}
            )
        ])

    elif tab == "tab-vir":
        if not VIR_COLS:
            return html.P("No virulence columns found.")
        vir_bin = df[["sample_id"]+VIR_COLS].set_index("sample_id")
        for c in VIR_COLS:
            vir_bin[c] = vir_bin[c].apply(
                lambda x: 0 if pd.isna(x) or str(x).strip() in ("-","None","nan","") else 1
            )
        vir_bin = vir_bin.loc[vir_bin.sum(axis=1).sort_values(ascending=False).index]
        fig = go.Figure(go.Heatmap(
            z=vir_bin.values, x=vir_bin.columns.tolist(), y=vir_bin.index.tolist(),
            colorscale=[[0,"#f0f3f4"],[1,"#8e44ad"]],
            showscale=False, xgap=3, ygap=1,
            hovertemplate="<b>%{y}</b>  %{x}<extra></extra>"
        ))
        fig.update_layout(
            height=max(500, len(vir_bin)*14+150),
            margin=dict(l=160,r=20,t=40,b=60),
            xaxis=dict(tickfont=dict(size=12)),
            yaxis=dict(tickfont=dict(size=9)),
            plot_bgcolor="white"
        )
        return html.Div([
            html.H3("Virulence loci (Kleborate)"),
            html.P("Purple = detected. Ybt=Yersiniabactin | Clb=Colibactin | "
                   "Iuc=Aerobactin | Iro=Salmochelin | RmpADC/rmpA2=hypermucoviscosity"),
            dcc.Graph(figure=fig)
        ])

    elif tab == "tab-snp":
        fig = go.Figure(go.Heatmap(
            z=snp_disp, x=snp_samps, y=snp_samps,
            colorscale="RdYlGn_r", zmid=25,
            hovertemplate="<b>%{y}</b> vs <b>%{x}</b><br>SNPs: %{z:.0f}<extra></extra>",
            colorbar=dict(title="SNP distance", thickness=15)
        ))
        fig.update_layout(
            height=max(600, len(snp_samps)*12+120),
            margin=dict(l=160,r=100,t=40,b=160),
            xaxis=dict(tickangle=-60, tickfont=dict(size=8)),
            yaxis=dict(tickfont=dict(size=8)),
            plot_bgcolor="white"
        )
        pairs = [
            {"Sample 1": snp_samps[i], "Sample 2": snp_samps[j],
             "SNP distance": int(snp_raw[i,j])}
            for i in range(len(snp_samps))
            for j in range(i+1, len(snp_samps))
            if snp_raw[i,j] <= 25
        ]
        children = [
            html.H3("Pairwise SNP distance matrix (core genome)"),
            html.P("Colour scale centred at 25 SNPs. Green = closely related. "
                   "Diagonal masked."),
            dcc.Graph(figure=fig)
        ]
        if pairs:
            close = pd.DataFrame(pairs).sort_values("SNP distance")
            children += [
                html.H4(f"Pairs within outbreak threshold (≤25 SNPs): {len(pairs)}"),
                dash_table.DataTable(
                    data=close.to_dict("records"),
                    columns=[{"name":c,"id":c} for c in close.columns],
                    sort_action="native", page_size=15,
                    style_table={"overflowX":"auto"}, style_cell={"fontSize":12},
                    style_header={"fontWeight":"bold","backgroundColor":"#ecf0f1"}
                )
            ]
        else:
            children.append(html.P("No pairs within ≤25 SNP threshold."))
        return html.Div(children)

    return html.P("Select a tab.")


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
