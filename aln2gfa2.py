import sys, time, gfapy, re
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

def verifier_graph(graph):
    liste = list(nx.connected_components(graph))
    f = open("verif_graph.txt", "w")
    for net in liste:
        for comp in list(net):
            f.write(comp + " ")
        f.write("\n")
    f.close()
    return;

def create_gfa2(graph):
    gfa2 = open(sys.argv[1][0:-4] + ".gfa2", "w")
    gfa = gfapy.Gfa(version='gfa2')
    gfa.append("H\tVN:Z:2.0")
    with open(sys.argv[1]) as sam:
        for aln in sam:
            aln = aln.split()
            if aln[0][0:3] == "@SQ":
                gfa.append("S\t" + aln[1].split(":")[1] + "\t" + aln[2].split(":")[1] + "\t*")
                continue
            elif aln[0][0] == "@":
                continue
            #gfa.segment[aln[0]].sequence = aln[9]
            cig = re.split('(\d+)', aln[5])
            D = 0
            I = 0
            taille_clip = 0
            beg2 = int(aln[3]) - 1
            for i in range(0, len(cig)):
                if cig[i] == "S":
                    taille_clip += int(cig[i - 1])
                elif cig[i] == "D":
                    D += int(cig[i - 1])
                elif cig[i] == "I":
                    I += int(cig[i - 1])
            taille_aln = int(len(aln[9])) - taille_clip
            end2 = beg2 + taille_aln - I
            if cig[2] == "S" or cig[2] == "H":
                beg1 = cig[1]
            else:
                beg1 = "0"
            if cig[-1] == "S":
                end1 = str(len(aln[9]) - int(cig[-2]))
            else:
                end1 = str(len(aln[9])) + "$"
            #if aln[1] == "16":
            #    beg1, end1 = end1, beg1
            if cig[-1] == "H":
                end2 -= int(cig[-2])
            reverse = "+"
            if aln[1] == "16":
                reverse = "-"
            #if gfa.segment[aln[2]] == end2:
            #    end2 += "$"
            #cig = gfapy.Alignment(aln[5])
            gfa.append("E\t*\t" + aln[0] + "+" + "\t" + aln[2] + reverse + "\t" + beg1 + "\t" + end1 + "\t" + str(beg2) + "\t" + str(end2) + "\t" + "*")
            #il faut trouver comment integrer le cigar
    for line in gfa.lines:
        gfa2.write(str(line) + "\n")
            #print("end1 calcul tip : {}".format(beg1 + cig.length_on_query()))
    return;

def main():
    G = create_graph(sys.argv[1])
    nx.write_edgelist(G, sys.argv[1][0:-4] + "_edgelist.csv")
    verifier_graph(G);
    #create_gfa2(G)
    return(0);

main()
print("Computed in {} seconds".format(time.time() - start_time))
