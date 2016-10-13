#include <iostream>
#include <set>
#include <map>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "cmdline/cmdline.h"

#include <ogdf/basic/Graph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/decomposition/BCTree.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/decomposition/StaticSPQRTree.h>
#include <ogdf/decomposition/Skeleton.h>

using namespace std;
using namespace ogdf;

unordered_map<node,string> id2contig;
unordered_map<string,node> revid2contig;
unordered_map<int,string> intid2contig;
class Link
{
public:
    int id;
    string contig_a;
    string contig_a_orientation;
    string contig_b;
    string contig_b_orientation;
    double mean;
    double stdev;
    int bundle_size;
	Link() {};
    Link(int id, string contig_a, string contig_a_orientation, string contig_b, string contig_b_orientation, double mean, double stdev);
	Link(int id, string contig_a, string contig_a_orientation, string contig_b, string contig_b_orientation, double mean, double stdev, int bundle_size);
	double getmean();
	double getstdev();
	string getlinkorientation();
    string getcontigs();
    string getfirstcontig();
    string getsecondcontig();
    string getfirstorietation();
    string getsecondorientation();
    int get_bundle_size();
    int getid();
};  

Link :: Link(int id, string contig_a, string contig_a_orientation, string contig_b, string contig_b_orientation, double mean, double stdev, int bundle_size)
{
	this->id = id;
	this->contig_a = contig_a;
	this->contig_b = contig_b;
	this->contig_a_orientation = contig_a_orientation;
	this->contig_b_orientation = contig_b_orientation;
	this->mean = mean;
	this->stdev = stdev;
    this->bundle_size = bundle_size;
}

Link :: Link(int id, string contig_a, string contig_a_orientation, string contig_b, string contig_b_orientation, double mean, double stdev)
{
    this->id = id;
    this->contig_a = contig_a;
    this->contig_b = contig_b;
    this->contig_a_orientation = contig_a_orientation;
    this->contig_b_orientation = contig_b_orientation;
    this->mean = mean;
    this->stdev = stdev;
}

string Link :: getfirstcontig()
{
    return this->contig_a;
}

string Link :: getsecondcontig()
{
    return this->contig_b;
}

string Link :: getfirstorietation()
{
    return this->contig_a_orientation;
}

string Link :: getsecondorientation()
{
    return this->contig_b_orientation;

}

int Link :: get_bundle_size()
{
    return this->bundle_size;
}

double Link :: getmean()
{
	return this->mean;
}

double Link :: getstdev()
{
	return this->stdev;
}

string Link :: getlinkorientation()
{
	return this->contig_a_orientation + this->contig_b_orientation;
}

string Link :: getcontigs()
{
    return contig_a +"$"+contig_b;
}

int Link :: getid()
{
    return this->id;
}

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

class Bicomponent {
	private:
	  std::set<int> memberNodes;
    
    std::set<int>::iterator nextIter(std::set<int>::iterator iter) {
      return ++iter;
    }
    
  public:
    Bicomponent(std::set<int> mn) {
      memberNodes = mn;
    } //constructor
};


int searchList (SList <edge> S,edge e) 
{
		int x = 0;
		
		for(SListIterator <edge> i = S.begin(); i.valid(); ++i, ++x) 
		{
				if(*i == e) 
				{
					return x;
				}
		}
		return -1;
}

string getTypeString(node &n, StaticSPQRTree &s) {
	std::string res = "unkown";	
	int type = s.typeOf(n);
	switch (type) {
		case 0:
			res = "S";
			break;
		case 1:
			res = "P";
			break;
		case 2:
			res = "R";
			break;
	}
	return res;
}

void write_dot(Graph G, map<int,int> sk2origin, string file,Skeleton &sk)
{
	ofstream of(getCharExpr(file));
	of << "digraph {"<<endl;
	edge e;
	forall_edges(e,G)
	{
		if(!sk.isVirtual(e))
		{
			int source = sk2origin[e->source()->index()];
			int target = sk2origin[e->target()->index()];
			of <<"\t"<<source<<"->"<<target<<endl;
		}
	}
	of<<"}";
}

vector<pair<int,int> > findTwoVertexCuts(Bicomponent &bicomp, Skeleton &sk, unordered_map<int,int> sk2orig, std::string type) 
{
	const Graph &G = sk.getGraph();
	int virtualCount;
	vector<pair<int,int> > pairs;
	edge e;
	
	node n1;
	const int nrNodes = G.numberOfNodes();
	//cout<<"Number of nodes = "<<nrNodes<<endl;
	int allnodes[nrNodes];
	int count = 0;
	
	forall_nodes(n1, G) {
		allnodes[count] = sk2orig[n1->index()];
		count++;
	}
	//cout<<"Done"<<endl;
	if (type == "R") {
		//cout<<"R"<<endl;
		//A virtual edge in an R node represents a two vertex cut
		forall_edges(e,G) {
			if (sk.isVirtual(e))
				pairs.push_back(make_pair(sk2orig[e->source()->index()], sk2orig[e->target()->index()]));
		} //forall edges
	}//if
	else if (type == "P") {
		//Node associated with p-nodes with two or more virtual edges are 2-vertex cuts
		//cout<<"P"<<endl;
		virtualCount = 0;
		forall_edges(e,G) {
			if (sk.isVirtual(e)) {
				virtualCount++;
				if (virtualCount > 1) {
					pairs.push_back(make_pair(sk2orig[e->source()->index()], sk2orig[e->target()->index()]));
					break;
				}//if
			}//if
		}//forall_edges
	}//else if
	else if (type == "S") 
	{
		//cout<<"S"<<endl;
		// A virtual edge in an S node represents a 2-vertex cuts
		forall_edges(e,G) {
			if (sk.isVirtual(e)) 
				pairs.push_back(make_pair(sk2orig[e->source()->index()], sk2orig[e->target()->index()]));
			
		} //forall edges
		
		// All non-adjacent nodes in an S-node are cut-vertices
		for (int i = 0; i < nrNodes-1; i++)
				for(int j = i+1; j < nrNodes; j++)
						if (find(pairs.begin(),pairs.end(),make_pair(allnodes[i], allnodes[j])) == pairs.end() && find(pairs.begin(),pairs.end(),make_pair(allnodes[j], allnodes[i])) == pairs.end())
							pairs.push_back(make_pair(allnodes[i], allnodes[j]));
	}//else if
	//cout<<pairs.size()<<endl;
	return pairs;
} //getTwoVertexCuts

std::set<int> getBiComponent(GraphCopy *GC, BCTree *p_bct, node bcTreeNode) 
{
	node n;
	edge e;
	std::set<int> memberNodes; // Members of the N-node
	
	const Graph &auxGraph = p_bct->auxiliaryGraph();
	//GraphCopy GC(auxGraph); 					                     //copy of original
	SList <edge> componentEdges = p_bct->hEdges(bcTreeNode); //edges in component bcTreeNode
	forall_edges (e, auxGraph) {							 						   //Check if edge belongs to component
		//cerr << "Testing edge " << e << endl;
		if (searchList(componentEdges,e) == -1) {			         //If not, delete edge from copy 
			GC->delEdge(GC->copy(e));
		}
	}
	forall_nodes(n, auxGraph)
	{								               //Delete nodes without edges
		if (!GC->copy(n)->degree()) 
		{
			//cerr << "Deleting node: " << n->index() << endl;
			GC->delNode(GC->copy(n));
		} //if
		else 
		{
		  int index = p_bct->original(n)->index();
		  memberNodes.insert(index);
		} //else
	}// forall_nodes
	return memberNodes;
}

node original(node &n, BCTree &bc, const GraphCopy &GC, Skeleton &sk)
{
	node np;
	np = bc.original(GC.original(sk.original(n)));
	return np;
}

int main(int argc, char* argv[])
{	
	cmdline ::parser pr;
    pr.add<string>("oriented_graph",'l',"list of oriented links",true,"");
    pr.add<string>("output",'o',"output file tow write sep pairs",true,"");
    pr.parse_check(argc,argv);
    Graph G;
    ifstream linkfile(getCharExpr(pr.get<string>("oriented_graph")));
    ofstream ofile(getCharExpr(pr.get<string>("output")));
    string line;
    //unordered_map<int, Link> linkmap;
   
    unordered_map<string,node> revid2contig;
    int contig_id = 1, linkid = 0;
    while(getline(linkfile,line))
    {
    	//cout<<line<<endl;	
    	string a,b,c,d;
    	double e,f;
        int g;
    	istringstream iss(line);
    	if(!(iss >> a >> b >> c >> d >> e >> f >> g))
    		break;
    	//Link l(linkid,a,b,c,d,e,f,g);
        //Link l(linkid,a,b,c,d,e,f);
        node first = 0, second = 0;
       	if(revid2contig.find(a) == revid2contig.end())
       	{
       		
       		first = G.newNode(contig_id);
       		id2contig[first] = a;
       		intid2contig[contig_id] = a;
       		revid2contig[a] = first;
       		contig_id++;
       	}
       	if(revid2contig.find(c) == revid2contig.end())
       	{
       		second = G.newNode(contig_id);
       		intid2contig[contig_id] = c;
       		revid2contig[c] = second;
       		id2contig[second] = c;
       		contig_id++;
       	}
       	// cout<<first<<"\t"<<second<<endl;
       	// G.newEdge((node)first,(node)second);
       	// cout<<"edge added"<<endl;
        //contigs2bundle[a+c] = g;
    }
    linkfile.close();
    //cout<<"Nodes: "<<G.numberOfNodes()<<endl;
    ifstream linkfile1(getCharExpr(pr.get<string>("oriented_graph")));
    while(getline(linkfile1,line))
    {
    	//cout<<line<<endl;	
    	string a,b,c,d;
    	double e,f;
        int g;
    	istringstream iss(line);
    	if(!(iss >> a >> b >> c >> d >> e >> f >> g))
    		break;
    	//Link l(linkid,a,b,c,d,e,f,g);
        //Link l(linkid,a,b,c,d,e,f);
        node first = revid2contig[a];
        node second = revid2contig[c];
       	cout<<first<<"\t"<<second<<endl;
       	edge x = G.newEdge(node(first),node(second));
       	//cout<<"edge added"<<endl;
        //contigs2bundle[a+c] = g;
    }

	cerr<<"Nodes: "<<G.numberOfNodes()<<endl;
	cerr<<"Edges: "<<G.numberOfEdges()<<endl;
	// GraphAttributes GA(G, GraphAttributes::nodeId);
	// bool ok = GraphIO::readGML(GA,G,"test_graph/oriented.gml");
	//since this is giving an error, lets just read tsv file and construct graph ourself

	// GraphIO::writeDOT(G,"tmp/original.dot");
	// cout<<ok<<endl;
	// if(ok)
	// {
	// 	cout<<"Graph loaded correctly!"<<endl;
	// }
	// else
	// {
	// 	cout<<"Graph loaded incorrectly!"<<endl;
	// }
	
	
	//decompose into connected components
	int nrCC = 0;
	NodeArray<int> node2cc(G);
	nrCC = connectedComponents(G, node2cc);
	cerr<<"Number of connected components = "<<nrCC<<endl;

	node startNodes[nrCC];
	int index = 0;
	node n;
	forall_nodes(n, G) 
	{
		if (node2cc[n] == index)
		{
			startNodes[index] = n;
			index++;
		}
		if (index == nrCC)
			break;
	}	
	set<int> memberNodes;
	unordered_map<int,int> sk2orig; // node mapping
	//Building BC tree for each component
	Graph G_new;
	int new_node_index = 1;
	map<int,vector<node> > nodemapping;
	for(int j = 0;j < nrCC; j++)
	{
		BCTree bc(G,startNodes[j]);
		BCTree *p_bct = &bc;
		cerr<<"Number of Biconnected Components = "<<bc.numberOfBComps()<<endl;

		if(bc.numberOfBComps() == 0)
		{
			continue;
			//do some special processing here
		}
		//Now, for each Biconnected Component, build SPQR tree 
		//Connected Components in auxgraph are the biconnected components of original graph


		const Graph &auxgraph = p_bct->auxiliaryGraph();
		cerr<<"graph made"<<endl;
		node bcTreeNode;
		forall_nodes(bcTreeNode,bc.bcTree())
		{

			if(bc.typeOfBNode(bcTreeNode) == 0)
			{
				GraphCopy GC(p_bct->auxiliaryGraph());
				memberNodes = getBiComponent(&GC,p_bct,bcTreeNode);
				cerr<<memberNodes.size()<<endl;
				Bicomponent bicomp(memberNodes);
				//cer<<"membernodes found"<<endl;
				//Now Generate SPQR tree for this component

				bool biconnected = isBiconnected(GC);
		        int  nrEdges     = GC.numberOfEdges();
		        bool loopfree    = isLoopFree(GC);
		        if(!biconnected || nrEdges <= 2 || !loopfree) 
		        {
		        	continue;
                    cerr << "Graph is not a valid input for SPQR-tree decomposition!" << endl;
                    cerr << "Reason(s):" << endl;
                    if (!biconnected)
                            cerr << "-> Graph is not biconnected" << endl;
                    if (nrEdges <= 2)
                            cerr << "-> Graph has "<< nrEdges << " edge(s). Should be more than 2." << endl;
                    if (!loopfree)
                            cerr << "-> Graph is not loop free" << endl;
           	
		        }

				StaticSPQRTree spqr(GC);
				//cout<<"SPQR generated"<<endl;
				const Graph &T = spqr.tree();
				//cout<<"SPQR tree made"<<endl;
				GraphIO::writeDOT(T,"tmp/spqr.dot");
				// cout<<"S nodes: "<<spqr.numberOfSNodes()<<endl;
				// cout<<"P nodes: "<<spqr.numberOfPNodes()<<endl;
				// cout<<"R nodes: "<<spqr.numberOfRNodes()<<endl;
				int c = 0;
				GraphCopy GCopy(T);
				node n,Nn,cn;
				forall_nodes(n, T) 
				{
					const Graph &Gn = spqr.skeleton(n).getGraph(); // Print the skeleton of a tree node to dis

					// Generate hash table: sk2orig[Skeleton node] = Original node 
					forall_nodes(Nn, Gn) 
					{
						cn = original(Nn,bc,GC,spqr.skeleton(n)); //Node in original graph G
						sk2orig[Nn->index()] = cn->index();
					}
									
					string type = getTypeString(n, spqr);	
					//Get 2-vertex cuts
					vector<pair<int, int> > pairs = findTwoVertexCuts(bicomp,spqr.skeleton(n) , sk2orig, type);
					for(int i = 0;i < pairs.size();i++)
					{
						ofile<<intid2contig[pairs[i].first]<<"\t"<<intid2contig[pairs[i].second];
						for(set<int> :: iterator it = memberNodes.begin(); it != memberNodes.end();++it)
						{
							ofile<<"\t"<<intid2contig[*it];
						}
						ofile<<endl;
					}
				}
				// forall_nodes(n,T)
				// {
					
				// 	int degree = n->degree();
				// 	const Graph &Gn = spqr.skeleton(n).getGraph();
				// 	//only if current node is leaf node, add it in new graph
					
					
				// 	vector<node> nodes_in_skel;
				// 	node Nn;
				// 	forall_nodes(Nn,Gn)
				// 	{	
				// 		node cn = original(Nn,bc,GC,spqr.skeleton(n));
				// 		sk2origin[Nn->index()] = cn->index();
				// 		if(degree == 1)
				// 			nodes_in_skel.push_back(cn);
				// 	}
				// 	nodemapping[new_node_index] = nodes_in_skel;
				// 	write_dot(Gn,sk2origin,"tmp/comp"+to_string(c)+".dot",spqr.skeleton(n));
				// 	if(degree == 1)
				// 	{
				// 		G_new.newNode(new_node_index);
				// 		new_node_index++;
				// 	}
				// 	c++;
				// 	//retrieve skeleton of nodes in SPQR tree
				// 	/*
				// 	const Graph &Gn = spqr.skeleton(n).getGraph();
				// 	node Nn;
				// 	forall_nodes(Nn,Gn)
				// 	{	
				// 		node cn = original(Nn,bc,GC,spqr.skeleton(n));
				// 		sk2origin[Nn->index()] = cn->index();
				// 	}
				// 	//write_dot(Gn,sk2origin,"tmp/comp"+to_string(c)+".dot",spqr.skeleton(n));
				// 	//find mapping of these nodes to original nodes in graph meaning : sk2origin[skeleton_node] = original_node
				// 	//c++;	
				// 	/*
				// 	map<int,int> :: iterator it;
				// 	for(it = sk2origin.begin();it != sk2origin.end();++it)
				// 	{
				// 		cout<<it->first<<"\t"<<it->second<<endl;
				// 	}
				// 	cout<<"=========="<<endl;
				// 	c = "2";
				// 	*/
				// }

			}
		}	
	}
	//add edges in this new graph based on original graph
	return 0;
}
