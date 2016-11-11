import networkx as nx
from Bio import SeqIO
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a','--assembly',help="initial assembly in fasta format",required=True)
	parser.add_argument('-e','--edges',help="sorted list of oriented edges",required=True)
	parser.add_argument('-o','--output',help="output file to write scaffolds",required=True)
	G = nx.Graph()
	args= parser.parse_args()
	handle = open(args.assembly, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

	contigs = set()
	seen_nodes = set()
	with open(args.edges,'r') as f:
		for line in f:
			attrs = line.split()
			v1 = attrs[0] + ':' + attrs[1]
			v2 = attrs[2] + ':' + attrs[3]
			mean = float(attrs[4])
			stdev = float(attrs[5])
			bsize = float(attrs[6])

			if v1 not in seen_nodes and v2 not in seen_nodes:
				G.add_edge(v1,v2,mean=mean,stdev=stdev,bsize=bsize)
				seen_nodes.add(v1)
				seen_nodes.add(v2)
				contigs.add(attrs[0])
				contigs.add(attrs[2])

	for contig in contigs:
		G.add_edge(contig+':B',contig+':E')


	records = []
	sum = 0
	scaff_id = 1
	scaffolded = {}
	G_copy = G.copy()
	for g in nx.connected_component_subgraphs(G):
		p = []
		for node in g.nodes():
			if g.degree(node) == 1:
				p.append(node)

		if len(p) == 2:
			path = nx.shortest_path(g,p[0],p[1])
			sum += len(path)/2
			curr_contig = ""
			for i in xrange(0,len(path)-1,2):
				contig = path[i].split(':')[0]
				scaffolded[contig] = 1
				first = path[i].split(':')[1]
				second = path[i+1].split(':')[1]
				if first+second == 'BE':
					curr_contig += str(record_dict[contig].seq)
				else:
					curr_contig += str(record_dict[contig].reverse_complement().seq)

				G_copy.remove_node(path[i])
				G_copy.remove_node(path[i+1])

			id = str("scaffold_"+str(scaff_id))
			records.append(SeqRecord(Seq(str(curr_contig),generic_dna), id = id))
			scaff_id += 1

	for g in nx.connected_component_subgraphs(G_copy):
		p = []
		for node in g.nodes():
			if g.degree(node) == 1:
				p.append(node)

		if len(p) == 2:
			path = nx.shortest_path(g,p[0],p[1])
			sum += len(path)/2
			curr_contig = ""
			for i in xrange(0,len(path)-1,2):
				contig = path[i].split(':')[0]
				scaffolded[contig] = 1
				first = path[i].split(':')[1]
				second = path[i+1].split(':')[1]
				if first+second == 'BE':
					curr_contig += str(record_dict[contig].seq)
				else:
					curr_contig += str(record_dict[contig].reverse_complement().seq)

				for j in xrange(0,100):
					curr_contig += 'N'

			id = str("scaffold_"+str(scaff_id))
			records.append(SeqRecord(Seq(str(curr_contig),generic_dna), id = id))
			scaff_id += 1

	for contig in record_dict:
		if contig not in scaffolded:
			#print contig
			id = str("scaffold_"+str(scaff_id))
			if len(str(record_dict[contig].seq)) >= 10000:
				records.append(SeqRecord(Seq(str(record_dict[contig].seq),generic_dna), id = id))
				scaff_id += 1
	f = open(args.output,'w')
	SeqIO.write(records,f,'fasta')
if __name__ == '__main__':
	main()
