import networkx as nx
import argparse

contig_coverage = {}

def find_plasmids(G,to_write):
    count = 0
    ofile = open(to_write,'w')
    for subg in nx.weakly_connected_component_subgraphs(G):
        is_cyclic = True
        for node in subg.nodes():   
            if subg.in_degree(node) == 1 and subg.out_degree(node) == 1:
                continue
            else:
                is_cyclic = False
                break

        if is_cyclic:
            ofile.write(str(subg.nodes())+'\n')


def find_tandem_repeats(G,to_write):
    sum = 0
    for key in contig_coverage:
        sum += contig_coverage[key]

    avg = sum*1.0/len(contig_coverage)
    print(avg)

    nodes = []
    for node in G.nodes():
        if G.in_degree(node) == 1 and G.out_degree(node) == 1:
            pre = G.predecessors(node)[0]
            suc = G.successors(node)[0]
            cov_pre = contig_coverage[pre]
            cov_suc = contig_coverage[suc]
            attrs = node.split('_')
            curr_cov = contig_coverage[node]
            if curr_cov >= 1.5* cov_pre and curr_cov >= 1.5*cov_suc:
                nodes.append(node)


    ofile = open(to_write,'w')
    for node in nodes:
        ofile.write(str(node)+'\n')


def find_interspersed_repeats(G,to_write):
    sum = 0
    for key in contig_coverage:
        sum += contig_coverage[key]

    avg = sum*1.0/len(contig_coverage)
    print(avg)

    nodes = []
    for node in G.nodes():
        if G.in_degree(node) >= 5 or G.out_degree(node) >= 5:
            if contig_coverage[node] >= 5*avg:
                nodes.append(node)


    ofile = open(to_write,'w')
    for node in nodes:
        ofile.write(str(node)+'\n')

def find_three_bubbles(G,to_write,seppairs):
    explored = {}
    ofile = open(to_write,'w')
    with open(seppairs,'r') as f:
        for line in f:
            attrs = line.split()
            if len(attrs) == 5:
                ofile.write(line)

def find_four_bubbles(G,to_write,seppairs):
    explored = {}
    ofile = open(to_write,'w')
    with open(seppairs,'r') as f:
        for line in f:
            attrs = line.split()
            if len(attrs) == 6:
                ofile.write(line)

def find_complex_bubbles(G,to_write,seppairs):
    ofile = open(to_write,'w')
    with open(seppairs,'r') as f:
        for line in f:
            attrs = line.split()
            if len(attrs) > 6:
                ofile.write(line)


def find_cycles(G,to_write):
    G = G.to_undirected()
    ofile = open(to_write,'w')
    basis = nx.cycle_basis(G)
    for each in basis:
        ofile.write(str(each)+'\n')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--working_dir',help='directory where all files exist',required=True)
    args = parser.parse_args()
    with open(args.working_dir+'/contig_coverage','r') as f:
        for line in f:
            attrs = line.split()
            contig_coverage[attrs[0]] = float(attrs[1])

    G = nx.read_gml(args.working_dir+'/oriented.gml')

    find_plasmids(G,args.working_dir+'/plasmids')
    #find_tandem_repeats(G,args.working_dir+'/tandem_repeats')
    #find_interspersed_repeats(G,args.working_dir+'/interspersed_repeats')
    find_three_bubbles(G,args.working_dir+'/three_bubbles',args.working_dir+'/bubbles.txt')
    find_four_bubbles(G,args.working_dir+'/four_bubbles',args.working_dir+'/bubbles.txt')
    find_complex_bubbles(G,args.working_dir+'/complex_bubbles',args.working_dir+'/bubbles.txt')
    find_cycles(G,args.working_dir+'/cycles')



if __name__ == '__main__':
    main()





