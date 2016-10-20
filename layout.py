import networkx as nx
from collections import deque
#from networkx.drawing.nx_agraph import write_dot
import operator
import argparse


revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N','a':'t','c':'g','g':'c','t':'a',' ':'',}[B] for B in x][::-1])

def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.iteritems():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa

def test_pair(G,source,sink,members):

	# if G.has_edge(source,sink) or G.has_edge(sink,source):
	# 	return []

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
		if v not in members:
			return False
		visited_nodes[v] = True
		if v == sink:
			at_sink = True
			continue
		for edge in G.in_edges(v):
			if edge not in visited:
				go_ahead = False
				break

		if go_ahead:
			visited[edge] = True
			for edge in G.out_edges(v):
				if edge not in visited:
					Q.appendleft(edge)
					visited[edge] = True
	
	
	
	if at_sink:
		return True
	else:
		return False


'''
This method finds out all shortest paths between source and sink in subg.
It returns a list of paths with first path being a heaviest path 
'''
def get_all_shortest_paths(subg, source, sink):
	all_paths = nx.all_simple_paths(subg,source,sink)
	id2path = {}
	id2weight = {}
	id = 1
	for path in all_paths:
		id2path[id] = path
		
		wt = 0
		for i in xrange(0,len(path)-1):
			wt += subg[path[i]][path[i+1]]['bsize']

		id2weight[id] = wt
		id += 1

	sorted_path = sorted(id2weight, key=lambda k: id2weight[k], reverse=True)
	ret = []
	for key in sorted_path:
		ret.append(id2path[key])

	return ret


'''
Instead of finding all shortest paths, remove node on heaviest shortest path and repeat
'''
def get_variants(subg,source,sink):
	#print subg.edges(data=True)
	subg1 = subg.copy()
	for u,v,data in subg1.edges(data=True):
		if data['bsize'] == 0:
			subg1[u][v]['bsize'] = 10
		else:
			subg1[u][v]['bsize'] = 1.0/data['bsize']
	paths = []
	path = nx.shortest_path(subg1,source,sink,weight='bsize')
	paths.append(path)
	# if len(path) == 2:
	# 	return paths
	# for each in path:
	# 	if each != source and each != sink:
	# 		subg1.remove_node(each)

	# while True:
	# 	print 'here'
	# 	try:
	# 		path = nx.shortest_path(subg1,source,sink,weight='bsize')
	# 		paths.append(path)
	# 		for each in path:
	# 			if each != source and each != sink:
	# 				subg1.remove(each)
	# 		if len(path) == 2:
	# 			return paths

	# 	except:
	# 		return paths

	return paths
'''	
This method takes a graph and makes it acyclic by removing lowest cost edge in a cycle
'''

def make_acyclic(G):
	G_copy = G.copy()
	F = []
	original_G = G.copy()
	while not nx.is_directed_acyclic_graph(G_copy):
		#iterate through cycles in G
		
		for cycle in nx.simple_cycles(G_copy):
			min_weight = 100000
			min_u = 0
			min_v = 0
			#Find minimum weight edge in the cycle, weight
			#here is bundle size
			#TODO: start with smallest cycle by sorting
			#print G.edges(data=True)
			for i in xrange(0,len(cycle)-1):
				u = cycle[i]
				v = cycle[i+1]
				if G[u][v]['bsize'] < min_weight:	
					min_weight = G[u][v]['bsize']
					min_u = u
					min_v = v
			if G[cycle[- 1]][cycle[0]]['bsize'] < min_weight:
				min_weight = G[cycle[-1]][cycle[0]]['bsize']
				min_u = cycle[-1]
				min_v = cycle[0]

			#reduce the edge weights by min_weight and remove the edge if its weight is 0
			if min_weight != 100000:
				for i in xrange(0,len(cycle)-1):
					u = cycle[i]
					v = cycle[i+1]
					G[u][v]['bsize'] -= min_weight
				
				G[cycle[-1]][cycle[0]]['bsize'] -= min_weight
				G.remove_edge(min_u,min_v)
				F.append((min_u,min_v,original_G.get_edge_data(min_u,min_v)))
				G_copy = G.copy()
				break

	#Now try adding edges from F to G, TODO do in non-increasing order

		if len(G.edges()) == 0:
			continue
		# if len(G.nodes()) == 0:
		# 	continue
		for edge in F:
			u = edge[0]
			v = edge[1]
			G.add_edge(u,v,edge[2])
			if not nx.is_directed_acyclic_graph(G):
				G.remove_edge(u,v)

	return G


'''
Helper to no_of_paths method
'''
def no_of_paths_helper(subg,source,sink,dp):
	# print "source = " + source
	# print "sink = " + sink
	if source == sink:
		return 1
	if dp[source] != -1:
		return dp[source]
	ret = 0
	for u,v in subg.out_edges(source):
		#print u,v
		ret += no_of_paths_helper(subg,v,sink,dp)
	dp[source] = ret
	return ret

'''
This method takes a DAG as input with source and sink and outputs number of paths
between source and sink
'''
def no_of_paths(subg,source,sink):
	#subg = nx.topological_sort(subg)
	dp = {}
	dp[source] = -1
	dp[sink] = -1
	for node in subg.nodes():
		dp[node] = -1
	return  no_of_paths_helper(subg,source,sink,dp)


'''
This  is main method
'''
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a','--assembly', help='Contig assembly', required=True)
	parser.add_argument('-g','--oriented_graph', help='Oriented Graph of Contigs', required=True)
	parser.add_argument('-s','--seppairs', help='Separation pairs detected in the graph', required=True)
	parser.add_argument('-o','--output', help='Output file for scaffold sequences', required=True)

	args = parser.parse_args()

	G = nx.read_gml(args.oriented_graph)
	#G = nx.read_gml("small.gml")
	nx.write_gexf(G,'original.gexf')
	pairmap = {}
	pair_list = []
	with open(args.seppairs,'r') as f:
		for line in f:
			attrs = line.split()
			if attrs[0] <= attrs[1]:
				key = attrs[0] +'$'+ attrs[1]
			else:
				key = attrs[1] +'$'+ attrs[0]
			pairmap[key] = attrs[2:]
			pair_list.append(key)

	validated = {}
	contig2id = {}
	cnt = 0
	#write_dot(G,'graph.dot')

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
	#pair_list = sorted(pairmap, key=lambda k: len(pairmap[k]), reverse=True)

	# for key in pair_list:
	# 	print pairmap[key]

	comp_to_id = {}
	id_to_comp = {}
	comp_to_pair = {}
	id_to_longest_path = {}
	comp2pairs = {}
	prev_comp = ''
	id = 1
	for key in pair_list:
		comp = pairmap[key]
		if comp[0] == prev_comp:
			continue

		comp_to_id[comp[0]] = str(id)
		comp2pairs[str(id)] = []
		id_to_comp[str(id)] = comp
		comp_to_pair[str(id)] = []
		id_to_longest_path[str(id)] = -1
		id += 1
		prev_comp = comp[0]

	for key in pair_list:
		c = pairmap[key][0]
		comp_id = comp_to_id[c]
		comp_to_pair[comp_id].append(key)

	valid_comps = {}

	for key in pair_list:
		contigs = key.split('$')
		'''
		First find the subgraph of bicomponent. Check if current source sink pair is longer that previously
		validated source sink pair. If yes then only validate current source sink pair.
		'''
		subg = G.subgraph(pairmap[key])
		comp_id = pairmap[key][0]
		comp_id = comp_to_id[comp_id]
		res = test_pair(G,contigs[0],contigs[1],pairmap[key])
		
		if res:
			cnt += 1
			#validated[contigs[0]] = 1
			source_and_sinks[contigs[0]] = 1
			source_and_sinks[contigs[1]] = 1
			#validated[contigs[1]] = 1
			#subg = G.subgraph(comp)
			valid_comps[comp_id] = 1
			

			
	source = {}
	sink = {}
	source_sink_to_comp = {}
	#print len(valid_comps)
	cnt = 0
	bubble_to_graph = {}
	for key in valid_comps:
		pairs = comp_to_pair[key]
		#print "Length of pairs = " + str(len(pairs))
		subg = G.subgraph(id_to_comp[key])
		if not nx.is_directed_acyclic_graph(subg):
			subg = make_acyclic(subg)
		if nx.is_directed_acyclic_graph(subg):

			#print subg.nodes()
			max_path = 0
			max_pair = -1
			#print pairs
			for pair in pairs:
				#print pair
				pair1 = pair.split('$')
				no_paths = no_of_paths(subg,pair1[0],pair1[1])
				if no_paths > max_path:
					max_path = no_paths
					max_pair = pair

			if max_pair != -1:
				# print "max_path = " + str(max_path)
				# print "max_pair = " + str(max_pair)
				# paths = get_variants(subg,max_pair.split('$')[0],max_pair.split('$')[1])
				# print paths
				cnt += 1
				bubble_to_graph[key] = subg
				valid_source_sink.append(max_pair)
				source[max_pair.split('$')[0]] = 1
				sink[max_pair.split('$')[1]] = 1
				source_sink_to_comp[max_pair.split('$')[0]] = key
				source_sink_to_comp[max_pair.split('$')[1]] = key
				for contig in id_to_comp[key]:
					in_bubble[contig] = 1
					validated[contig] = 1
		# else:
		# 	subg = make_acyclic

	#print cnt

	'''
	Here, find now the new graph by collapsing bubbles
	TODO: Preserve node and edge attributes from the original non-collapsed graph
	'''
	#node to info map
	node_info = {}
	for node in G.nodes(data=True):
		node_info[node[0]] = node[1] 
	G_new = nx.DiGraph()

	# print source
	# print sink
	# for each in source:
	# 	print len(G.in_edges(each))

	# for each in sink:
	# 	print len(G.out_edges(each))

	# print source
	# print sink
	for key in valid_comps:
		G_new.add_node(str(key))
	for u,v,data in G.edges(data=True):
		if u not in validated and v not in validated:
			G_new.add_edge(u,v,data)

	for node in G.nodes():
		if node not in source and node not in sink:
			for each in source:
				if G.has_edge(node,each):
					#print 'here'
					data = G.get_edge_data(node,each)
					G_new.add_edge(node,source_sink_to_comp[each],data)
			for each in sink:
				if G.has_edge(each,node):
					#print 'here'
					data = G.get_edge_data(each,node)
					G_new.add_edge(source_sink_to_comp[each],node,data)

	for s in source:
		for t in sink:
			if source_sink_to_comp[s] != source_sink_to_comp[t]:
				if G.has_edge(s,t):
					data = G.get_edge_data(s,t)
					G_new.add_edge(source_sink_to_comp[s],source_sink_to_comp[t],data)
				if G.has_edge(t,s):
					data = G.get_edge_data(t,s)
					G_new.add_edge(source_sink_to_comp[t],source_sink_to_comp[s],data)
	for node in G_new.nodes(data=True):
		if node[0] in node_info:
			info = node_info[node[0]]
			for each in info:
				node[1][each] = info[each]
			node[1]['type'] = 'contig'
		else:
			node[1]['type'] = 'bubble'
			#node[1]['size'] = len(bubble_to_graph[node[0]].nodes())
	'''
	Output the simplified Graph
	'''
	# for node in G_new.nodes(data=True):
	# 	#print node
	# 	m = node[1]
	# 	node[1]['color'] = colmap[node[0]]
	#nx.set_node_attribute(G_new,'color',colmap)
	print len(G_new.nodes())
	print len(G_new.edges())
	#nx.write_gexf(G_new,'simplified.gexf')
	#write_dot(G_new,'simplified.dot')
	nx.write_gml(G_new,'simplified.gml')

	'''
	In this simplified, for each weakly connected component, find out the heaviest linear path. If path
	goes through the bubble, choose the heaviest path in the bubble and continue
	'''

	alternative_contigs = [] #this stores all variants. Tag these as variants while writing to file
	primary_contigs = []
	for subg in nx.weakly_connected_component_subgraphs(G_new):
		# print len(subg.nodes())
		# print 'here'
		#First get all edges
		edges = subg.edges(data=True)
		#sort edges by weights
		sorted_edges = sorted(edges,key = lambda tup: tup[2]['bsize'], reverse=True)
		#print sorted_edges
		#create a new graph
		G_sorted = nx.Graph()
		#add edges to this graph until for is created, this will be undirected graph and it will have
		#'B' and 'E' nodes
		nodes = set()
		for edge in sorted_edges:
			u = edge[0]
			v = edge[1]
			data = edge[2]
			orientation  = data['orientation']
			u = u + '$' + orientation[0]
			v = v + '$' + orientation[1]
			if not G_sorted.has_node(u) and not G_sorted.has_node(v):
				G_sorted.add_edge(u,v,data)
				nodes.add(u.split('$')[0])
				nodes.add(v.split('$')[0])

		#add edges between B and E nodes of same contig
		for node in nodes:
			G_sorted.add_edge(node+'$B',node+'$E')

		#print len(G_sorted.edges())
		#print G_sorted.edges(data=True)
		#now trace out all linear paths in this, each will be a scaffold
		for small_subg in nx.connected_component_subgraphs(G_sorted):
			#print small_subg.edges()
			p = []
			for node in small_subg.nodes():
				if small_subg.degree(node) == 1:
					p.append(node)


			if len(p) == 2:
				path = nx.shortest_path(small_subg,p[0],p[1])

				#if path has a bubble node, insert the contigs on the heaviest path on the bubble
				new_path = []
				new_path_ind = 0
				for i in xrange(1,len(path),2):
					node = path[i].split('$')[0]
					if node not in bubble_to_graph:
						new_path.append(path[i-1])
						new_path.append(path[i])
						new_path_ind += 2
						continue

					bubble_graph = bubble_to_graph[node]
					curr_source = ''
					curr_sink = ''
					for node in bubble_graph.nodes():
						if node in source:
							curr_source = node
						if node in sink:
							curr_sink = node
					try:
						bubble_paths = get_variants(bubble_graph,curr_source,curr_sink)
					except:
						continue
						
					heaviest = bubble_paths[0]
					if len(new_path) > 0:
						# print new_path[new_path_ind-1].split('$')[0]
						# print new_path
						if G.has_edge(new_path[new_path_ind-1].split('$')[0],heaviest[0]):
							continue
						else:
							heaviest.reverse()

					for each in heaviest:
						#print 'appending heaviest'
						orient = G.node[each]['orientation']
						if orient == 'FOW':
							new_path.append(each+'$B')
							new_path.append(each+'$E')
							new_path_ind += 2

						if orient == 'REV':
							new_path.append(each+'$E')
							new_path.append(each+'$B')	
							new_path_ind += 2

					for i in xrange(1,len(bubble_paths)):
						print 'in alternate path'
						alt_path = []
						curr_path = bubble_paths[i]
						for each in curr_path:
							o_node = G.node
							if G.node[each]['orientation'] == 'FOW':
								alt_path.append(each+'$B')
								alt_path.append(each+'$E')

							if G.node[each]['orientation'] == 'REV':
								alt_path.append(each+'$E')
								alt_path.append(each+'$B')

						alternative_contigs.append(alt_path)
				primary_contigs.append(new_path)
				#print new_path

	assembly = open(args.assembly,'r')
	sequences = parse_fasta(assembly.readlines())
	ofile = open(args.output,'w')

	scaffold_id = 1
	for scaffold in primary_contigs:
		scaff_string = ''
		for i in xrange(0,len(scaffold) - 1,2):
			curr = scaffold[i]
			next = scaffold[i+1]
			contig = curr.split('$')[0]
			start = curr.split('$')[1]
			end = next.split('$')[1]
			if start == 'B' and end == 'E':
				scaff_string += sequences[contig]
			else:
				scaff_string += revcompl(sequences[contig])

		chunks = [scaff_string[i:i+80] for i in xrange(0,len(scaff_string),80)]
		ofile.write('>scaffold_'+str(scaffold_id)+'\n')
		for chunk in chunks:
			ofile.write(chunk+'\n')
		scaffold_id += 1


if __name__ == '__main__':
	main()


