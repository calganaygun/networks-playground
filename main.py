import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

try:
    import graph_tool.all as gt
except:
    print("graph-tool is not installed, please install it to count the motifs.")
    print("https://graph-tool.skewed.de/static/doc/index.html#installing-graph-tool")


from tqdm import tqdm

def load_datasets(dataset_paths):
    datasets = []
    for path in dataset_paths:
        G = nx.read_edgelist(path)
        print(f"Loaded {path} with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
        print(f"Graph ID for outputs: {G.number_of_nodes() + G.number_of_edges()}")
        datasets.append(G)
        print("--------------------")
    return datasets


def create_random_graph(G, i):
    prefix = G.number_of_nodes() + G.number_of_edges()
    # if file exist load and return
    # if not generate and save
    file_exist = os.path.exists(f"random_graphs/{prefix}/random_graph_{i}.edges")
    if file_exist:
        return nx.read_edgelist(f"random_graphs/{prefix}/random_graph_{i}.edges")
    degree_sequence = [G.degree(n) for n in G.nodes()]
    G_random = nx.random_degree_sequence_graph(degree_sequence, seed=i, tries=10)
    # save as file
    nx.write_edgelist(G_random, f"random_graphs/{prefix}/random_graph_{i}.edges")
    print(f"Generated random graph {i} for {prefix}.")


def generate_random_graphs(G, num_graphs):
    graphs = []
    # genreate random graphs parallely
    for i in range(num_graphs):
        graphs.append(create_random_graph(G, i))
    return graphs


def count_motifs(G):
    # I used graph-tool library to count the motifs
    # convert networkx graph to graph-tool graph
    # beacsue netsci only giving 3 node motifs counts
    g = gt.Graph(directed=False)
    node_map = g.new_vertex_property("int")
    for node in G.nodes():
        node_map[node] = g.add_vertex()
    for edge in G.edges():
        g.add_edge(node_map[edge[0]], node_map[edge[1]])

    # count 3-node and 4-node connected motifs using the k-core decomposition algorithm
    motifs_3 = gt.motifs(g, 3)
    motifs_4 = gt.motifs(g, 4)
    motifs_list = motifs_3[0] + motifs_4[0]
    motifs_counts = motifs_3[1] + motifs_4[1]
    return motifs_list, motifs_counts


def calculate_z_scores(motifs_counts_real, motifs_counts_randoms):
    z_scores = {}
    for i in range(len(motifs_counts_real)):
        random_motif_counts = [
            motifs_counts_random[i] for motifs_counts_random in motifs_counts_randoms
        ]
        real_motif_count = motifs_counts_real[i]
        mean_random = np.mean(random_motif_counts)
        std_random = np.std(random_motif_counts)
        z_score = (real_motif_count - mean_random) / std_random
        z_scores[i] = z_score
    return z_scores


def plot_z_scores(z_scores, graph_id):
    plt.bar(z_scores.keys(), z_scores.values())
    plt.xlabel("Motifs")
    plt.ylabel("Z-score")
    plt.title("Z-score of motifs")
    # show every x ticks
    plt.xticks(np.arange(0, len(z_scores), 1))
    # save plot
    plt.savefig(f"outputs/{graph_id}_z_scores.jpg")
    # clear for next plot
    plt.clf()


def save_motifs_graphs(motifs_list, graph_id):
    for i, n in enumerate(motifs_list):
        gt.graph_draw(n, output=f"outputs/motifs/{graph_id}_{i}.png")
    print(f"Motifs graphs saved as png files to: outputs/{graph_id}_*.png")


def save_motifs_counts_table(motifs_counts_real, motifs_counts_randoms, graph_id):
    df = pd.DataFrame([motifs_counts_real, np.mean(motifs_counts_randoms, axis=0)], columns=range(8))
    df.index = ["Real Graph", "Random Graphs (Mean)"]
    df.to_csv(f"outputs/{graph_id}.csv")
    print(f"Motif counts saved as table to: outputs/{graph_id}.csv")


def plot_motifs_counts(df, graph_id):
    df.plot(kind="bar", title=f"Motifs counts for graph ID {graph_id}")
    plt.xlabel("Motifs")
    plt.ylabel("Count")
    plt.savefig(f"outputs/{graph_id}.jpg")
    plt.clf()
    print(f"Motif counts plot saved as jpg file to: outputs/{graph_id}.jpg")

MOTIF_LISTS = [0] * 8

def main():
    dataset_paths = ["bio-DM-LC/bio-DM-LC.edges", "rt-twitter-copen/rt-twitter-copen.mtx"]
    datasets = load_datasets(dataset_paths)
    for G in datasets:
        graph_id = G.number_of_nodes() + G.number_of_edges()
        print(f"Performing tasks (a) and (b) for graph ID {graph_id}")

        num_graphs = 100
        motifs_counts_randoms = []

        print(f"Counting motifs for real graph...")
        motifs_list, motifs_counts_real = count_motifs(G)
        print(f"Motifs counts for real graph: {motifs_counts_real}\n")

        print(f"Generating {num_graphs} random graphs...")
        random_graphs = generate_random_graphs(G, num_graphs)

        print(f"Counting motifs for {num_graphs} random graphs...")
        for G_random in tqdm(random_graphs):
            motifs_list_random, motifs_count_random = count_motifs(G_random)

            # match motif list with MOTIF_LISTS for reordering the motifs
            # gt.isomorphism method can be used for matching
            motifs_count = [0] * 8
            for i, motif in enumerate(motifs_list_random):
                if MOTIF_LISTS[i] == 0:
                    MOTIF_LISTS[i] = motif
                else:
                    for j, motif_ in enumerate(MOTIF_LISTS):
                        if gt.isomorphism(motif, motif_):
                            motifs_count[j] = motifs_count_random[i]
                            break
            
            motifs_counts_randoms.append(motifs_count)
        motifs_counts_randoms_mean = np.mean(motifs_counts_randoms, axis=0)
        print(f"Motifs counts for random graphs (mean): {motifs_counts_randoms_mean}\n")

        save_motifs_counts_table(motifs_counts_real, motifs_counts_randoms, graph_id)
        df = pd.DataFrame([motifs_counts_real, motifs_counts_randoms_mean], columns=range(8))
        df.index = ["Real Graph", "Random Graphs (Mean)"]
        plot_motifs_counts(df, graph_id)

        save_motifs_graphs(motifs_list, graph_id)

        z_scores = calculate_z_scores(motifs_counts_real, motifs_counts_randoms)
        plot_z_scores(z_scores, graph_id)
        print(f"Z-scores plot saved as jpg file to: outputs/graph_{graph_id}.jpg")
        print("--------------------------------------------------")

if __name__ == "__main__":
    main()