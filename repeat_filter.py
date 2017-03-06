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
    skewed_edges[node] = s_count

#centralities
centralities = {}
with open(sys.argv[4],'r') as f:
	for line in f:
		attrs = line.split()
		centralities[attrs[0].strip("\"")] = float(attrs[1])

#lengths
with open(sys.argv[5],'r') as f:
	for line in f:
		attrs = line.split()
		contig_length[attrs[0]] = int(attrs[1])

repeats = {}
mean = np.mean(centralities.values())
stdev = np.std(centralities.values())

for contig in centralities:
	if centralities[contig] >= mean + 3*stdev and mean + 3*stdev != 0:
		repeats[contig] = 1

p_coverage = np.percentile(contig_coverage.values(),75)
p_invalidated = np.percentile(contig2links.values(),75)
p_degree = np.percentile(contig_degree.values(),75)
p_skewed = np.percentile(skewed_edges.values(),75)

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
	if contig in skewed_outliers and contig in links_outliers and contig in degree_outliers:
		other_repeats[contig] = 1


with open(sys.argv[2],'r') as f:
	for line in f:
		attrs = line.split()
		dist = float(attrs[4])
		if attrs[0] in repeats or attrs[2] in repeats:
			continue
		if attrs[0] in other_repeats or attrs[2] in other_repeats:
			continue
		if dist < 0:
			if abs(dist) >= contig_length[attrs[0]] or abs(dist) >= contig_length[attrs[2]]:
				continue
			else:
				print line.strip()
		print line.strip()	
