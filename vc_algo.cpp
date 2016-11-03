#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <string>
#include <numeric>
#include <map>
#include <algorithm>
#include <stdlib.h>

#include "cmdline/cmdline.h"

#include <NetworKit/graph/Graph.h>
#include <NetworKit/io/GMLGraphReader.h>
#include <NetworKit/centrality/ApproxBetweenness.h>
#include <NetworKit/centrality/Betweenness.h>
#include <NetworKit/centrality/ApproxBetweenness2.h>


using namespace std;
using namespace NetworKit;

char* getCharExpr(string s)
{
	char *a=new char[s.size()+1];
	a[s.size()]=0;
	memcpy(a,s.c_str(),s.size());
	return a;
}

int main(int argc, char* argv[]) {
	
	cmdline::parser p;
	p.add<string>("oriented_graph",'g',"give graph of oriented contigs",true,"");
	p.add<string>("repeats_file",'r',"file where to write repeated nodes",true,"");
	p.add<double>("epsilon",'f',"sampling parameter",false,0.2);
	//clock_t start = clock();
	p.parse_check(argc,argv);
	GMLGraphReader reader;
	Graph g = reader.read(getCharExpr(p.get<string>("oriented_graph")));

	cerr << "Nodes: " << g.nodes().size() << endl;
	cerr << "Edges: " << g.edges().size() << endl;
	vector<pair<unsigned long,unsigned long>> all_edges = g.edges();
	double delta = 0.1;
	double epsilon = p.get<double>("epsilon");

	map<int,string> id2contigs;

	ifstream cfile;
	cfile.open(getCharExpr(p.get<string>("oriented_graph")));
	string line;
	bool isidpresent = false;
	int id;
	string contig;
	while(getline(cfile,line))
	{
		istringstream iss(line);

		if(line.find("id") != string::npos)
		{
			istringstream iss(line);
			string a;
			int b;
			iss >> a >> b;
			id = b;
			isidpresent = true;
		}
		if(line.find("label") != string::npos)
		{
			istringstream iss(line);
			string a,b;
			iss>>a>>b;
			//cout<<id<<"\t"<<b<<endl;;
			id2contigs[id] = b;
			isidpresent = false;
		}
		if(line.find("edge") != string::npos)
			break;
	}

	
	ApproxBetweenness centrality(g, delta, epsilon,0);
	//Closeness centrality(g,false,true);
	centrality.run();
	ofstream myfile;
	
	myfile.open(getCharExpr(p.get<string>("repeats_file")));
	vector<pair<node, double>> nodeRanks(g.nodes().size());
	nodeRanks = centrality.ranking();
	vector<double> scores = centrality.scores();
	
	vector<double> scores1 = scores;
	sort(scores1.begin(),scores1.end());
	double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
	double mean = sum / scores.size();
	int first_quartile = scores1.size()/4;
	int median = scores1.size()/2;
	int third_quartile = first_quartile + median;
	double range = scores1.at(third_quartile) - scores1.at(first_quartile);
	double filter = scores1.at(third_quartile) + 1.5*range;

	double sq_sum = std::inner_product(scores.begin(), scores.end(), scores.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / scores.size() - mean * mean);
	cerr<<"Mean: "<<mean<<endl;
	cerr<<"Std Dev: "<<stdev<<endl;
	cerr<<"bound: "<<mean + 3*stdev<<endl;
	cerr<<"filter: "<<filter<<endl;
	double bound = mean + 3*stdev;
	//double bound = 1e-5;
	
	cerr<<"done with algorithm"<<endl;

	
	for (int i = 0; i < g.nodes().size(); i++) {
		int isrepeat = 0;
		if(nodeRanks.at(i).second >= bound) 			 	
		//if(nodeRanks.at(i).second >= bound)
		{
			isrepeat = 1;
		}
		string String = static_cast<ostringstream*>( &(ostringstream() << nodeRanks.at(i).first) )->str();
		if(isrepeat == 1)
			myfile <<id2contigs[nodeRanks.at(i).first]<<"\t"<<nodeRanks.at(i).second<<endl;
	}
	
	return 0;
}

