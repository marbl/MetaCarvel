import networkx as nx

G = nx.DiGraph()


def main():
	G = nx.read_gml("c++/RF_oriented.gml")
	#G = nx.read_gml("oriented.gml")

	print nx.is_directed_acyclic_graph(G)
	#Doing this for each connected component
	original_G = G.copy()
	#Feedback edge set 
	F = []
	node_info = {}
	for node in G.nodes(data=True):
		node_info[node[0]] = node[1]
	G_copy = G.copy()
	# print len(G.nodes())
	# print len(G.edges())
	#Iterate until G_copy is acyclic
	for G in nx.connected_component_subgraphs(original_G.to_undirected()):
		G_copy = nx.DiGraph()
		#G_temp  = nx.DiGraph()
		for u,v,data in G.edges(data=True):
			if original_G.has_edge(u,v):
				G_copy.add_edge(u,v,attr_dict=data)
				#G_temp.add_edge(u,v,data=data)
			else:
				G_copy.add_edge(v,u,attr_dict=data)
				#G_temp.add_edge(v,u,data=data)

		G_temp = G_copy.copy()
		G = G_copy.copy()	
		if len(G.nodes()) < 2:
			continue

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

		#Now since we have DAG, find topological sortin
		top_sort = nx.topological_sort(G)
		for i in xrange(0,len(top_sort)):
			print node_info[top_sort[i]]['label']
		for i in xrange(0,len(top_sort)-1):
			if G.has_edge(top_sort[i],top_sort[i+1]):
				print G[top_sort[i]][top_sort[i+1]]['mean']
				print G[top_sort[i]][top_sort[i+1]]['stdev']

		print "============================"
if __name__ == '__main__':
	main()