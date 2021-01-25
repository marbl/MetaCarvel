import sys
import numpy as np
import networkx as nx

contig_coverage = {}
contig_degree = {}
contig2links = {}
central_nodes = {}
contig_length = {}
#contig coverages
with open(sys.argv[1],'r') as f:
    for line in f:
        attrs = line.split()
        contig_coverage[attrs[0]] = float(attrs[1])

#contig degree, bundled links
G = nx.MultiGraph()
with open(sys.argv[2],'r') as f:
    for line in f:
        attrs = line.split()
        G.add_edge(attrs[0],attrs[2])

for node in G.nodes():
    contig_degree[node] = G.degree(node)

#invalidated links
with open(sys.argv[3],'r') as f:
    for line in f:
        attrs = line.split()
        contig2links[attrs[0]] = int(attrs[1])

#skewed links
skewed_edges = {}
for node in G.nodes():
    s_count = 0
    for neighs in G.neighbors(node):
        if node in contig_coverage and neighs in contig_coverage:
            if contig_coverage[node] >= 2*contig_coverage[neighs]:
                    s_count += 1
    skewed_edges[node] = s_count*1.0/len(list(G.neighbors(node)))

#centralities
centralities = {}
with open(sys.argv[4],'r') as f:
    for line in f:
        attrs = line.split()
        centralities[attrs[0]] = float(attrs[1])

#centralities = nx.betweenness_centrality(G)
#lengths
with open(sys.argv[5],'r') as f:
    for line in f:
        attrs = line.split()
        contig_length[attrs[0]] = int(attrs[1])

repeats = {}
mean = np.mean(list(centralities.values()))
stdev = np.std(list(centralities.values()))

for contig in centralities:
    repeats[contig] = 1

p_coverage = np.percentile(list(contig_coverage.values()),75)
p_invalidated = np.percentile(list(contig2links.values()),75)
p_degree = np.percentile(list(contig_degree.values()),75)
p_skewed = np.percentile(list(skewed_edges.values()),75)

avg_coverage = np.mean(list(contig_coverage.values()))
other_repeats = {}
coverage_outliers = {}
links_outliers = {}
skewed_outliers = {}
degree_outliers = {}

for contig in contig_coverage:
    if contig_coverage[contig] >= p_coverage:
        coverage_outliers[contig] = 1

for contig in skewed_edges:
    if skewed_edges[contig] >= p_skewed:
        skewed_outliers[contig] = 1

for contig in contig2links:
    if contig2links[contig] >= p_invalidated:
        links_outliers[contig] = 1

for contig in contig_degree:
    if contig_degree[contig] >= p_degree:
        degree_outliers[contig] = 1

for contig in coverage_outliers:
    if contig in links_outliers and contig in degree_outliers:
        other_repeats[contig] = 1

repeat_contigs = set()
for key in repeats:
    repeat_contigs.add(key)

for key in other_repeats:
    repeat_contigs.add(key)


with open(sys.argv[2],'r') as f:
    for line in f:
        attrs = line.split()
        dist = float(attrs[4])
        #if contig_coverage[attrs[0]] >= 3.5*avg_coverage or contig_coverage[attrs[2]] >=33.5*avg_coverage:
        #   continue
        if mean != 0 and stdev != 0 and attrs[0] in repeats or attrs[2] in repeats:
            continue
        if attrs[0] in other_repeats or attrs[2] in other_repeats:
            continue
        if dist < 0:
            if abs(dist) >= contig_length[attrs[0]] or abs(dist) >= contig_length[attrs[2]]:
                if abs(dist) >= contig_length[attrs[0]]:
                    repeat_contigs.add(attrs[0])
                if abs(dist) >= contig_length[attrs[2]]:
                    repeat_contigs.add(attrs[2])
                continue
            else:
                print(line.strip())
                continue
        print(line.strip())

ofile = open(sys.argv[6],'w')
#pool =  ThreadPool(cpus)
for each in repeat_contigs:
    ofile.write(each+'\n')
