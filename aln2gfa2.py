import sys, time
import networkx as nx
import matplotlib.pyplot as plt

#usage : python aln2gfa2.py aln.sam

start_time = time.time()

def create_graph(fichier):
    f = open(fichier, "r")
    sam = f.readlines()
    futur_graph = {}
    for aln in sam:
        if aln[0] == "@":
            continue;
        aln = aln.split()
        #aln[0] = aln[0].split("_")[0] + "_" + aln[0].split("_")[1]
        #aln[2] = aln[2].split("_")[0] + "_" + aln[2].split("_")[1]
        if aln[0] in futur_graph:
            futur_graph[aln[0]].append([aln[2], aln[3], aln[12].split(":")[2]])
        else:
            futur_graph[aln[0]] = [[aln[2], aln[3], aln[12].split(":")[2]]]
    G = nx.Graph()
    G.add_nodes_from(futur_graph.keys())
    for ref in futur_graph:
        for mate in futur_graph[ref]:
            G.add_edge(ref, mate[0])
    return(G);

def visualiser_graph(graph):
    nx.draw_networkx(graph, with_labels=True, node_size=500)
    plt.show()
    return;

def main():
    G = create_graph(sys.argv[1])
    nx.write_edgelist(G, "edgelist.csv")
    return(0);

main()
print("Computed in {} seconds".format(time.time() - start_time))
