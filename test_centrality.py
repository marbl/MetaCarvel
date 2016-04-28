import networkx as nx
import numpy as np
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--oriented_graph", help="Oriented graph")
	#parser.add_argument("-r", "--repeats", help="repeat_file")
	args = parser.parse_args()

	G = nx.MultiGraph()
	G = nx.read_gml(args.oriented_graph)
	#print len(G.edges())
	nodes = {}
	for node,data in G.nodes_iter(data=True):
		nodes[node] = data['label']

	centrality = nx.betweenness_centrality(G,weight='weight')
	num = centrality.values()
	mean = np.mean(num)
	stdev = np.std(num)
	'''
	for w in sorted(centrality, key=centrality.get, reverse=True):
		if centrality[w] >= mean + 3*stdev:
  			print str(nodes[w])+"\t"+ str(centrality[w]) + "\t1"
  		else:
  			print str(nodes[w])+"\t"+ str(centrality[w]) + "\t0"
	'''
	for w in sorted(centrality, key=centrality.get, reverse=True):
                if centrality[w] >= mean + 3*stdev:
                        print str(nodes[w])
if __name__ == '__main__':
	main()
