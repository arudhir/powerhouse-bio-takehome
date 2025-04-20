# stringutils.py
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt

# Paths
STRING_DIR = Path("~/Desktop/databases/STRING").expanduser()
LINK_FILE = STRING_DIR / "9606.protein.links.detailed.v12.0.txt.gz"
ALIAS_FILE = STRING_DIR / "9606.protein.aliases.v12.0.txt.gz"

def load_string_network(score_threshold=700):
    print(f"📂 Loading STRING interactions from: {LINK_FILE.name} (score ≥ {score_threshold})")
    df = pd.read_csv(LINK_FILE, sep=" ", compression="gzip")
    df = df[df['combined_score'] >= score_threshold]
    G = nx.from_pandas_edgelist(
        df, source="protein1", target="protein2", edge_attr="combined_score"
    )
    return G

def load_alias_map():
    print(f"📂 Loading STRING alias map from: {ALIAS_FILE.name}")
    df = pd.read_csv(ALIAS_FILE, sep="\t", compression="gzip", names=["string_protein_id", "alias", "source"])
    return df[df["source"] == "Ensembl_HGNC"]

def map_genes_to_string_ids(gene_list, alias_df):
    print(f"🔍 Mapping {len(gene_list)} gene symbols to STRING IDs...")
    mapping = alias_df[alias_df['alias'].isin(gene_list)]
    symbol_to_string = dict(zip(mapping['alias'], mapping['string_protein_id']))
    string_to_symbol = dict(zip(mapping['string_protein_id'], mapping['alias']))
    return symbol_to_string, string_to_symbol

def subset_graph(G, gene_list, symbol_to_string, distance=1):
    print(f"🔍 Subsetting graph to genes of interest with distance ≤ {distance}...")
    valid_ids = {symbol_to_string[g] for g in gene_list if g in symbol_to_string}
    neighborhood = set()
    for node in valid_ids:
        if node in G:
            neighborhood |= set(nx.ego_graph(G, node, radius=distance).nodes)
    return G.subgraph(neighborhood).copy()

def annotate_colors(G, pex_map, mt_map):
    print("🎨 Annotating node colors...")
    node_colors = {}
    for sid in G.nodes:
        in_pex = sid in pex_map.values()
        in_mt = sid in mt_map.values()
        if in_pex and in_mt:
            node_colors[sid] = "red"
        elif in_pex:
            node_colors[sid] = "green"
        elif in_mt:
            node_colors[sid] = "blue"
        else:
            node_colors[sid] = "gray"
    nx.set_node_attributes(G, node_colors, name="color")
    return G

def relabel_graph_nodes(G, alias_df):
    print("🔁 Relabeling nodes to gene symbols...")
    alias_df = alias_df[alias_df["source"] == "Ensembl_HGNC"]
    id_to_symbol = dict(zip(alias_df["string_protein_id"], alias_df["alias"]))
    id_to_symbol = {k: v for k, v in id_to_symbol.items() if k in G.nodes}
    return nx.relabel_nodes(G, id_to_symbol)

def draw_colored_graph(G, layout=None, figsize=(12, 12), title=None, node_size=300):
    if layout is None:
        layout = nx.spring_layout(G, seed=42)
    node_colors = [G.nodes[n].get("color", "gray") for n in G.nodes]

    plt.figure(figsize=figsize)
    nx.draw_networkx(
        G,
        pos=layout,
        with_labels=True,
        node_color=node_colors,
        edge_color="lightgray",
        node_size=node_size,
        font_size=8,
        font_weight="bold"
    )
    if title:
        plt.title(title, fontsize=16)
    plt.axis("off")
    plt.tight_layout()
    plt.show()