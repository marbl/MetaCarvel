import networkx as nx
from collections import deque
from networkx.drawing.nx_agraph import write_dot


'''
This method validates if separation pair given by SPQR tree is valid source-sink pair for bubble
The algorithm is same as the one in Marygold
'''

def test_pair(G,source,sink):

	if G.has_edge(source,sink) or G.has_edge(sink,source):
		return []

	visited = {}
	visited_nodes = {}
	# visited_nodes = set()
	# visited_nodes.add(source)
	visited_nodes[source] = True
	Q = deque()
	at_sink = False
	for edge in G.out_edges(source):
		Q.appendleft(edge)
		visited[edge] = True

	while not len(Q) == 0:
		#print len(Q)
		go_ahead = True
		curr_edge = Q.pop()
		u = curr_edge[0]
		v = curr_edge[1]
		visited_nodes[v] = True
		if v == sink:
			at_sink = True
			continue
		for edge in G.in_edges(v):
			if edge not in visited:
				go_ahead = False
				break

				break

		if go_ahead:
			visited[edge] = True
			for edge in G.out_edges(v):
				if edge not in visited:
					Q.appendleft(edge)
					visited[edge] = True
	
	
	
	if at_sink:
		return visited_nodes.keys()
	else:
		return []


'''
This method finds out all shortest paths between source and sink in subg.
It returns a list of paths with first path being a heaviest path 
'''
def get_all_shortest_paths(subg, source, sink):
	print "source = " + source
	print "sink = " + sink
	paths = []
	'''
	Invert the edge weights and find shortest paths
	'''

	for u,v,data in subg.edges(data=True):
		subg[u][v]['bsize'] = 1.0/data['bsize']

	'''
	Call shortest path now, remove nodes in that except source and sink and then repeat
	'''
	if subg.has_edge(source,sink):
		subg.remove_edge(source,sink)
	'''
	First remove any edge between source and sink
	'''

	while True:
		try:
			s_path = nx.shortest_path(subg,source,sink,weight='bsize')
		except:
			return paths
		#print s_path
		if s_path == []:
			break

		paths.append(s_path)

		'''
		Remove other nodes from the graph
		'''

		for node in s_path:
			if node != source and node != sink:
				subg.remove_node(node)


	return paths


'''
This method takes a graph G as input and removes transitive edges using 
Gene Myer's algorithm from Fragment assembly paper
'''
def transitive_reduction(G):
	
	for u in G.nodes():
		for v in G.nodes():
			for w in G.nodes():
				if G.has_edge(u,w) and G.has_edge(w,v) and G.has_edge(u,v):
					G.remove_edge(u,v)
	return G 

def main():

	G = nx.read_gml("small.gml")
	#G = nx.read_gml("small.gml")
	nx.write_gexf(G,'original.gexf')
	pairmap = {}

	with open('smallseppairs','r') as f:
		for line in f:
			attrs = line.split()
			if attrs[0] <= attrs[1]:
				key = attrs[0] +'$'+ attrs[1]
			else:
				key = attrs[1] +'$'+ attrs[0]
			pairmap[key] = attrs[2:]

	validated = {}
	contig2id = {}
	cnt = 0
	write_dot(G,'graph.dot')

	# for key in pairmap:
	# 	print len(pairmap[key])

	in_bubble = {}
	valid_source_sink = []
	all_bubble_paths = {} #stores all heaviest paths in bubble
	source_and_sinks = {}
	'''
	Here, first validate each source sink pair. To do this, sort them with largest number of nodes in the
	biconnected component.
	'''
	pair_list = sorted(pairmap, key=lambda k: len(pairmap[k]), reverse=True)
	print pair_list
	for key in pair_list:
		contigs = key.split('$')
		if contigs[0] not in validated and contigs[1] not in validated:
			res = test_pair(G,contigs[0],contigs[1])
			if len(res) != 0:
				print 'validated'
				cnt += 1
				valid_source_sink.append((contigs[0],contigs[1]))
				validated[contigs[0]] = 1
				source_and_sinks[contigs[0]] = 1
				source_and_sinks[contigs[1]] = 1
				validated[contigs[1]] = 1
				comp = []
				for contig in res:
					in_bubble[contig] = 1
					validated[contig] = 1
					comp.append(contig)

				subg = G.subgraph(comp)
				
				subg = subg.to_undirected()

				# print subg.nodes()
				# print subg.edges()

				'''
				Find all heaviest paths between source and sink in this bubble and store it for layout
				'''
				subg1 = subg.copy()
				paths = get_all_shortest_paths(subg,contigs[0],contigs[1])
				# print len(paths)
				# print len(paths[0])
				# print len(subg1.nodes())
				all_bubble_paths[key] = paths
				#print "============================"
			# else:
				# print "invalid"
	print cnt

	'''
	Here, find now the new graph by collapsing bubbles
	'''

	G_new = nx.DiGraph()
	colmap = {}
	typemap = {}
	bub_count = 1000
	for each in valid_source_sink:
		#print 'here'
		G_new.add_edge(each[0],str(bub_count))
		G_new[each[0]][str(bub_count)]['bsize'] = 50000
		G_new.add_edge(str(bub_count),each[1])	
		G_new[str(bub_count)][each[1]]['bsize'] = 50000
		typemap[each[0]] = 'source'
		typemap[each[1]]= 'sink'
		typemap[str(bub_count)] = 'bubble'
		colmap[each[0]] = 'red'
		colmap[each[1]] = 'red'
		colmap[str(bub_count)] = 'green'
		bub_count += 1

	for edge in G.edges():
		if edge[0] not in validated and edge[1] not in validated:
			data= G.get_edge_data(edge[0],edge[1])
			G_new.add_edge(edge[0],edge[1],bsize=data['bsize'])
			typemap[edge[0]] = 'contig'
			typemap[edge[1]]= 'contig'
			colmap[edge[0]] = 'black'
			colmap[edge[1]] = 'black'
		else:
			if edge[0] in validated and edge[1] not in validated or edge[1] in validated and edge[0] not in validated:
				
				data= G.get_edge_data(edge[0],edge[1])
				G_new.add_edge(edge[0],edge[1],bsize=data['bsize'])
				if edge[0] in validated:
					typemap[edge[1]] = 'contig'
					colmap[edge[0]] = 'red'
					colmap[edge[1]] = 'black'
				else:
					typemap[edge[0]] = 'contig'
					colmap[edge[0]] = 'black'
					colmap[edge[1]] = 'red'

	for node in source_and_sinks:
		for node1 in source_and_sinks:
			if G.has_edge(node,node1):
				G_new.add_edge(node,node1)
				colmap[node] = 'red'
				colmap[node1] = 'red'
			if G.has_edge(node1,node):
				G_new.add_edge(node1,node)
				colmap[node] = 'red'
				colmap[node1] = 'red'

	for node in G_new.nodes(data=True):
		try:
			node[1]['type'] = typemap[node[0]]
		except:
			node[1]['type'] = 'contig'

	'''
	Output the simplified Graph
	'''
	for node in G_new.nodes(data=True):
		#print node
		m = node[1]
		node[1]['color'] = colmap[node[0]]
	#nx.set_node_attribute(G_new,'color',colmap)
	print len(G_new.nodes())
	print len(G_new.edges())
	nx.write_gexf(G_new,'simplified.gexf')
	write_dot(G_new,'simplified.dot')

	'''
	Make the new simplified graph undirected
	'''
	#G_new = G_new.to_undirected()

	G_reduced = nx.DiGraph()

	for subg in nx.weakly_connected_component_subgraphs(G_new):
		subg = transitive_reduction(subg)
		for u,v in subg.edges():
			G_reduced.add_edge(u,v)

	nx.write_gexf(G_reduced,'t_reduced.gexf')

if __name__ == '__main__':
	main()
