import networkx as nx
import argparse
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--graph", help='bundled graph')
parser.add_argument("-l","--length",help="contig length")
parser.add_argument("-o","--output",help="output file")
args = parser.parse_args()
G = nx.Graph()
cpus = multiprocessing.cpu_count()
with open(args.graph,'r') as f:
    for line in f:
        attrs = line.split()
        G.add_edge(attrs[0],attrs[2],mean=float(attrs[4]),stdev=float(attrs[5]),bsize=int(attrs[6]),ori=attrs[1]+attrs[3])


contig_length = {}
with open(args.length,'r') as f:
    for line in f:
        attrs = line.split()
        if attrs[0] in G.nodes():
            contig_length[attrs[0]] = int(attrs[1])

nx.set_node_attributes(G,'length',contig_length)
#print contig_length
repeat_nodes = {}

def get_centrality(subg):
    centralities = nx.betweenness_centrality(subg)
    mean = np.mean(centralities.values())
    stdev = np.std(centralities.values())
    for node in centralities:
        if centralities[node] >= mean + 3*stdev:
            repeat_nodes[node] = centralities[node]

def centrality_wrapper(graph):
    pool =  ThreadPool(cpus)
    for subg in nx.connected_component_subgraphs(graph):
        #get_centrality(subg)
        if len(subg.nodes()) >= 50:
            result = pool.map(get_centrality,[subg])
    pool.close()
    pool.join()
G_copy = G.copy()

ofile = open(args.output,'w')

for i in xrange(3):
    centrality_wrapper(G_copy)
    for node in repeat_nodes:
        if G_copy.has_node(node):
            G_copy.remove_node(node)
	    ofile.write(str(node)+'\t'+str(repeat_nodes[node])+'\n')


 
#for u,v,data in G_copy.edges(data=True):
#    print u +"\t"+data[u][v]['ori'][0]+v+"\t"+data[u][v]['ori'][1]+"\t"+str(data[u][v]["mean"])+"\t"+str(data[u][v]["stdev"])+"\t"+str(data[u][v]["bsize"])
#nx.write_gml(G_copy,args.output)

