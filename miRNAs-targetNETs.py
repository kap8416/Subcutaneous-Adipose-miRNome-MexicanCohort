import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import networkx as nx
from pathlib import Path

# ------------------------- Config -------------------------
UP_XLSX = "metascapemiRNAsup.xlsx"
DOWN_XLSX = "metascapemiRNAsdown.xlsx"

# Palettes (extend if tienes más miRNAs)
PALETTE_UP   = ["#D73027", "#FC8D59", "#F46D43", "#FDAE61", "#FEE090", "#E6550D", "#A63603"]
PALETTE_DOWN = ["#4575B4", "#74ADD1", "#3288BD", "#ABD9E9", "#5E4FA2", "#1F78B4", "#084594"]

# Ring radii y colores por categoría de compartición
RINGS = {"2": 3.8, "3": 2.3, "4+": 1.1}
RING_COLORS = {"2": "#78C679", "3": "#FDDC6C", "4+": "#F16913"}  # verde, amarillo, naranja

# Grosor de aristas por categoría
EDGE_LW = {"2": 0.5, "3": 0.9, "4+": 1.4}

# --------------------- Utils de datos ---------------------
def load_clean(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    # elimina filas totalmente vacías y quita espacios
    df = df.dropna(how="all").applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # elimina columnas completamente vacías (por si vinieran)
    df = df.dropna(axis=1, how="all")
    return df

def build_edges(df: pd.DataFrame) -> pd.DataFrame:
    """Convierte columnas (miRNAs) -> filas (edges miRNA-gene) y filtra a genes compartidos (>=2)."""
    edges = []
    for mir in df.columns:
        for gene in df[mir].dropna():
            edges.append((mir, gene))
    if not edges:
        return pd.DataFrame(columns=["miRNA", "Gene", "Category"])

    edges_df = pd.DataFrame(edges, columns=["miRNA", "Gene"])
    gene_counts = edges_df["Gene"].value_counts()

    # filtrar genes compartidos
    shared = gene_counts[gene_counts >= 2]
    edges_df = edges_df[edges_df["Gene"].isin(shared.index)].copy()

    # categorizar por # de miRNAs que comparten el gen
    def cat(g):
        c = gene_counts[g]
        if c == 2: return "2"
        if c == 3: return "3"
        return "4+"
    edges_df["Category"] = edges_df["Gene"].map(cat)
    return edges_df

# --------------------- Plot circular ----------------------
def plot_circular_network(edges_df: pd.DataFrame,
                          title: str,
                          mirna_palette: list,
                          out_png: str) -> None:
    """Dibuja red circular con anillos de compartición y aristas coloreadas por miRNA de origen."""
    fig, ax = plt.subplots(figsize=(11, 11), facecolor="none")
    ax.set_title(title, fontsize=18, fontweight="bold", loc="left")
    ax.axis("off")

    if edges_df.empty:
        # nada compartido
        ax.text(0.5, 0.5, "No shared targets (≥2)", ha="center", va="center", transform=ax.transAxes,
                fontsize=14, fontweight="bold")
        plt.savefig(out_png, dpi=300, transparent=True, bbox_inches="tight")
        plt.close(fig)
        return

    # Grafo y atributos
    G = nx.Graph()
    for _, r in edges_df.iterrows():
        mir, gene, cat = r["miRNA"], r["Gene"], r["Category"]
        # Guardamos el miRNA "origin" para colorear la arista correctamente
        G.add_edge(mir, gene, mirna_origin=mir, category=cat)

    # Paleta por miRNA (en el orden que aparecen)
    miRNAs = edges_df["miRNA"].unique().tolist()
    miRNA_colors = dict(zip(miRNAs, mirna_palette[:len(miRNAs)]))

    # Layout circular: miRNAs afuera, genes por anillos
    pos = {}
    R_outer = 5.2
    ang_mir = np.linspace(0, 2*np.pi, len(miRNAs), endpoint=False)
    for i, mir in enumerate(miRNAs):
        pos[mir] = (np.cos(ang_mir[i]) * R_outer, np.sin(ang_mir[i]) * R_outer)

    # genes por categoría en anillos concéntricos
    for cat in ["2", "3", "4+"]:
        genes_cat = edges_df.loc[edges_df["Category"] == cat, "Gene"].unique()
        if len(genes_cat) == 0:
            continue
        ang = np.linspace(0, 2*np.pi, len(genes_cat), endpoint=False)
        r = RINGS[cat]
        for i, gene in enumerate(genes_cat):
            pos[gene] = (np.cos(ang[i]) * r, np.sin(ang[i]) * r)

    # Dibujar aristas (color del miRNA de origen)
    for u, v, d in G.edges(data=True):
        mir_origin = d.get("mirna_origin")
        cat = d.get("category", "2")
        color = miRNA_colors.get(mir_origin, "#888888")
        x1, y1 = pos[u]; x2, y2 = pos[v]
        ax.plot([x1, x2], [y1, y2], color=color, lw=EDGE_LW.get(cat, 0.6), alpha=0.75)

    # Dibujar nodos
    # miRNAs (afuera)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=miRNAs,
        node_color=[miRNA_colors[m] for m in miRNAs],
        node_size=950, linewidths=1.2, edgecolors="white", ax=ax
    )
    # genes por categoría
    for cat in ["2", "3", "4+"]:
        nodes = edges_df.loc[edges_df["Category"] == cat, "Gene"].unique().tolist()
        if nodes:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=nodes,
                node_color=RING_COLORS[cat],
                node_size=160, linewidths=0.6, edgecolors="white", ax=ax
            )

    # Etiquetas: TODOS los nodos (miRNAs en negritas)
    outline = [pe.withStroke(linewidth=2.0, foreground="white", alpha=0.9)]
    # miRNAs
    for m in miRNAs:
        ax.text(pos[m][0], pos[m][1], m, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color="black", path_effects=outline)
    # genes
    genes_all = edges_df["Gene"].unique().tolist()
    for g in genes_all:
        ax.text(pos[g][0], pos[g][1], g, ha="center", va="center",
                fontsize=7.2, color="black", path_effects=outline)

    # Anillos guía (suaves)
    for r in [RINGS["2"], RINGS["3"], RINGS["4+"], R_outer]:
        circ = plt.Circle((0, 0), r, fill=False, ls="--", lw=0.35, color="lightgray", alpha=0.35)
        ax.add_artist(circ)

    # Leyendas: miRNAs + categorías
    mir_handles = [mpatches.Patch(color=miRNA_colors[m], label=m) for m in miRNAs]
    cat_handles = [
        mpatches.Patch(color=RING_COLORS["2"],  label="Shared by 2 miRNAs"),
        mpatches.Patch(color=RING_COLORS["3"],  label="Shared by 3 miRNAs"),
        mpatches.Patch(color=RING_COLORS["4+"], label="Shared by ≥4 miRNAs"),
    ]
    leg1 = ax.legend(handles=mir_handles, loc="upper center", bbox_to_anchor=(0.5, -0.08),
                     ncol=min(len(mir_handles), 5), frameon=False, fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=cat_handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=3, frameon=False, fontsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)

# ------------------------ Run ----------------------------
def main():
    # UP
    up_df = load_clean(UP_XLSX)
    edges_up = build_edges(up_df)
    plot_circular_network(edges_up, "A) Upregulated miRNAs",
                          PALETTE_UP, "network_up_alllabels.png")
    # DOWN
    down_df = load_clean(DOWN_XLSX)
    edges_down = build_edges(down_df)
    plot_circular_network(edges_down, "B) Downregulated miRNAs",
                          PALETTE_DOWN, "network_down_alllabels.png")
    print("OK -> network_up_alllabels.png, network_down_alllabels.png")

if __name__ == "__main__":
    main()
