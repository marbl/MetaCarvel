#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <queue>

#include "cmdline/cmdline.h"

using namespace std;



const int FOW = 1, REV = 2, NIL = 0;

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

bool pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {
  return firstElem.second < secondElem.second;

}

struct SortLinkByBundle
{
    bool operator()(const Link& lhs, const Link& rhs) const
    {
        return lhs.bundle_size > rhs.bundle_size;
    }
};


struct SortLinkByNeighborSize
{
    map<string,int> contig2length;
    SortLinkByNeighborSize(map<string,int>& contig2length1)
    {
        contig2length = contig2length1;
    }

    bool operator()(const Link& lhs, const Link& rhs)
    {
        return contig2length[lhs.contig_b] > contig2length[rhs.contig_b];
    }
};

struct SortLinkByDegree
{
    map<string,int> contig2degree;
    SortLinkByDegree(map<string,int>& contig2degree1)

    {
        contig2degree = contig2degree1;
    }

    bool operator()(const Link& lhs, const Link& rhs)
    {
        return contig2degree[lhs.contig_b] > contig2degree[rhs.contig_b];
    }
};

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

//add bundle size

class Node
{
public:
    string contig;
    int length;
    int degree;
    Node () {}
    Node(string contig, int length);
    Node(string contig, int length, int degree);
};

Node :: Node(string contig, int length)
{
    this->contig = contig;
    this->length = length;
}

Node :: Node(string contig, int length, int degree)
{
    this->contig = contig;
    this->length = length;
    this->degree = degree;
}

//comparators for node class

struct MoreThanByLength
{
    bool operator()(const Node& lhs, const Node& rhs) const
    {
        return lhs.length > rhs.length;
    }
};

struct MoreThanByDegree
{
    bool operator()(const Node& lhs, const Node& rhs) const
    {
        return lhs.degree > rhs.degree;
    }
};


map<string, vector<Link> > adjacency;
map<string, vector<Link> > revadjacency;
map<string, int> ctg2orient;
map<int, bool> invalidlinks;
map<string, int> contig2length;
map<string, int> contigs2bundle;
map<string, int> contig2degree;
ofstream invalidfile("hmp/invalid_counts");

int findorientation(string node_to_orient)
{
    cout<<"finding orientation for node "<<node_to_orient<<endl;
    int curr_fow = 0, curr_rev = 0;
    if(adjacency.find(node_to_orient) != adjacency.end())
    {
        vector<Link> neighbors = adjacency[node_to_orient];
        for(int i = 0;i < neighbors.size();i++)
        {
            Link link = neighbors[i];
            if(invalidlinks.find(link.getid()) == invalidlinks.end())
            {
                string neigh = link.getsecondcontig();
                if(ctg2orient.find(neigh) != ctg2orient.end() && ctg2orient[neigh] != NIL)
                {
                    int orientation = ctg2orient[neigh];
                    if(orientation == FOW)
                    {
                        if(link.getlinkorientation() == "EB")
                            curr_fow += link.get_bundle_size();
                            //curr_fow++;

                        if(link.getlinkorientation() == "BB")
                            curr_rev += link.get_bundle_size();
                            //curr_rev++;
                    }
                    if(orientation == REV)
                    {
                        if(link.getlinkorientation() == "EE")
                            curr_fow += link.get_bundle_size();
                            //curr_fow++;

                        if(link.getlinkorientation() == "BE")
                            curr_rev += link.get_bundle_size();
                            //curr_rev++;
                    }
                }
            }
        }
    }
    if(revadjacency.find(node_to_orient) != adjacency.end())
    {
        //retrieve adjacency list
        vector<Link> neighbors = revadjacency[node_to_orient];
        //check if any of the neighbors is oriented, if yes then use that to orient current node
        for(int i = 0;i < neighbors.size();i++)
        {
            Link link = neighbors[i];
            if(invalidlinks.find(link.getid()) == invalidlinks.end())
            {
                string neigh = link.getfirstcontig();
                if(ctg2orient.find(neigh) != ctg2orient.end() && ctg2orient[neigh] != NIL)
                {
                    int orientation = ctg2orient[neigh];
                    //cerr<<"Using node "<<neigh<<" with orientation "<<orientation<<endl;
                    if(orientation == FOW)
                    {
                        //cerr<<link.getlinkorientation()<<endl;
                        if(link.getlinkorientation() == "EB")
                            curr_fow += link.get_bundle_size();
                            //curr_fow++;

                        if(link.getlinkorientation() == "EE")
                            curr_rev += link.get_bundle_size();
                            //curr_rev++;
                    }

                    if(orientation == REV)
                    {
                        if(link.getlinkorientation() == "BB")
                            curr_fow += link.get_bundle_size();
                            //curr_fow++;
                        
                        if(link.getlinkorientation() == "BE")
                            curr_rev += link.get_bundle_size();
                            //curr_rev++;
                        
                    }
                }
            }
        }
        //cerr<<"currfow = "<<curr_fow<<endl;
        //cerr<<"currev  = "<<curr_rev<<endl; 
        //cout<<"FOW = "<<curr_fow<<"   REV = "<<curr_rev<<endl;

        if(curr_fow >= curr_rev)
        {
            //cerr<<node_to_orient<<" orientation = fwd"<<endl;
            return FOW;
        }
        else
        {
            //cerr<<node_to_orient<<" orientation = rev"<<endl;
            return REV;
        }
    }
    return NIL;
}

void invalidatelinks(string v,int orientation)
{ 
    int count = 0;
    cerr<<"invalidating..."<<v<<endl;
    if(adjacency.find(v) != adjacency.end())
    {
        vector<Link> neighbors = adjacency[v];
        for(int i = 0;i < neighbors.size();i++)
        {
            Link link = neighbors[i];
            string target = link.getsecondcontig();
            //cerr<<"calculating for "<<target<<endl;
            if(ctg2orient[target] != NIL)
            {
                int neighorientation = ctg2orient[target];
                if(orientation == FOW && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "EB")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BE")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == FOW && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "EE")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "BB")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
            }
        }
    }
    
    if(revadjacency.find(v) != revadjacency.end())
    {
        vector<Link> neighbors = revadjacency[v];
        for(int i = 0;i < neighbors.size();i++)
        {
            Link link = neighbors[i];
            string target = link.getfirstcontig();
            //cerr<<"calculating for "<<target<<endl;
            if(ctg2orient[target] != NIL)
            {
                int neighorientation = ctg2orient[target];
                if(orientation == FOW && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "EB")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cout<<link.getlinkorientation()<<"\t"<<"BE"<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BE")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "EE")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cout<<orientation<<"\t"<<neighorientation<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == FOW && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BB")
                    {
                        invalidlinks[link.getid()] = true;
                        count += link.get_bundle_size();
                        //cout<<orientation<<"\t"<<neighorientation<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
            }
        }
    }
    invalidfile<<v<<"\t"<<count<<endl;
}

int get_degree(string start)
{
    int degree = 0;
    if(adjacency.find(start) != adjacency.end())
    {
        vector<Link> neighbors = adjacency[start];
        degree += neighbors.size();
    }
    if(revadjacency.find(start) != revadjacency.end())
    {
        vector<Link> neighbors = revadjacency[start];
        degree += neighbors.size();
    }
    return degree;
}

void bfs(string start, string strategy)
{
  //Priority Queue based BFS using length as priority
    if(strategy == "length")
    {
        std :: priority_queue<Node,vector<Node>, MoreThanByLength> Q;
        Node n(start,contig2length[start]);
        Q.push(n);
        while(!Q.empty())
        {
            Node n = Q.top();
            Q.pop();
            //cout<<Q.size()<<endl;
            string u = n.contig;
            //sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByBundle());
            sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByNeighborSize(contig2length));
            //sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByDegree(contig2degree));
            for(int i = 0;i < int(adjacency[u].size());i++)
            {
                Link l = adjacency[u][i];
                string v = l.getsecondcontig();
                if(ctg2orient[v] == NIL)
                {
                    int orientation = findorientation(v);
                    ctg2orient[v] = orientation;
                    invalidatelinks(v,orientation);
                    Node n(v,contig2length[v]);
                    Q.push(n);
                }
                
                else
                {
                    invalidatelinks(v,ctg2orient[v]);
                }
            }
        }
    }
    //priority based BFS using degree as priority
    if(strategy == "degree")
    {
        std :: priority_queue<Node,vector<Node>, MoreThanByDegree> Q;
        Node n(start,contig2length[start],get_degree(start));
        Q.push(n);
        while(!Q.empty())
        {
            Node n = Q.top();
            Q.pop();
            string u = n.contig;
            //sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByBundle());
            sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByDegree(contig2degree));
            //sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByNeighborSize(contig2length));
            for(int i = 0;i < int(adjacency[u].size());i++)
            {
                Link l = adjacency[u][i];
                string v = l.getsecondcontig();
                if(ctg2orient[v] == NIL)
                {
                    int orientation = findorientation(v);
                    ctg2orient[v] = orientation;
                    invalidatelinks(v,orientation);
                    Node n(v,contig2length[v],get_degree(v));
                    Q.push(n);
                }
                
                else
                {
                    invalidatelinks(v,ctg2orient[v]);
                }
            }
        }
    }

    //Choose node by bundle size
    if(strategy == "bsize")
    {   
        queue<string> Q;
        Q.push(start);
        while(!Q.empty())
        {
            string u = Q.front();
            Q.pop();
            sort(adjacency[u].begin(),adjacency[u].end(),SortLinkByBundle());
            for(int i = 0;i < int(adjacency[u].size());i++)
            {
                Link l = adjacency[u][i];
                string v = l.getsecondcontig();
                if(ctg2orient[v] == NIL)
                {
                    int orientation = findorientation(v);
                    ctg2orient[v] = orientation;
                    invalidatelinks(v,orientation);
                    Q.push(v);
                }
                
                else
                {
                    invalidatelinks(v,ctg2orient[v]);
                }
            }
        }
    }
}


void get_contig_length(string file)
{
    ifstream lenfile(getCharExpr(file));
    string line;
    while(getline(lenfile,line))
    {
        istringstream iss(line);
        string contig;
        int len;
        iss >> contig >> len;
        contig2length[contig] = len;
    }
}

string get_unoriented_node_by_length()
{
    //cout<<"inside method"<<endl;
    map<string, int> :: iterator it;
    int max_len = -1;
    string max_contig = "";
    for(it = ctg2orient.begin(); it != ctg2orient.end();++it)
    {
        //cout<<it->first<<"\t"<<it->second<<endl;
        if(it->second == NIL)
        {
            if(contig2length[it->first] > max_len)
            {
                max_len = contig2length[it->first];
                max_contig = it->first;
            }
        }
    }
    if(max_contig != "")
        return max_contig;
    return "done";
}

string get_unoriented_node_by_degree()
{
    map<string, int> :: iterator it;
    int max_degree = -1;
    string max_contig = "";
    for(it = ctg2orient.begin(); it != ctg2orient.end();++it)
    {
        if(it->second == NIL)
        {
            if(get_degree(it->first) > max_degree)
            {
                max_degree = contig2length[it->first];
                max_contig = it->first;
            }
        }
    }
    if(max_contig != "")
        return max_contig;
    return "done";
}

int main(int argc, char* argv[])
{
	
    cmdline ::parser pr;
    pr.add<string>("bundled_graph",'l',"list of bundled links",true,"");
    pr.add<string>("contig_length",'c',"contig lengths",true,"");
    pr.add("length",'\0',"sort contigs by size");
    pr.add("bsize",'\0',"sort contigs by bundle size");
    pr.add("degree",'\0',"sort contigs by degree");
    pr.add<string>("output",'o',"output file",true,"");
    pr.parse_check(argc,argv);
    map<string,double> contig2coverage;
    get_contig_length(pr.get<string>("contig_length"));
    string line;
    /*
    ifstream covfile("contig_coverage");
    while(getline(covfile,line))
    {
        string contig;
        double cov;
        istringstream iss(line);
        if(!(iss >> contig >> cov))
            break;
        contig2coverage[contig] = cov;
    }*/
    ifstream linkfile(getCharExpr(pr.get<string>("bundled_graph")));
    ofstream ofile(getCharExpr(pr.get<string>("output")));
    int linkid = 0;
    map<int, Link> linkmap;
    while(getline(linkfile,line))
    {
    	string a,b,c,d;
    	double e,f;
        int g;
    	istringstream iss(line);
    	if(!(iss >> a >> b >> c >> d >> e >> f >> g))
    		break;
    	Link l(linkid,a,b,c,d,e,f,g);
        //Link l(linkid,a,b,c,d,e,f);
        ctg2orient[a] = NIL;
        ctg2orient[c] = NIL;
    	linkmap[linkid] = l;
        //contigs2bundle[a+c] = g;
    	if(adjacency.find(a) == adjacency.end())
    	{
    		vector<Link> neighbors;
    		neighbors.push_back(l);
    		adjacency[a] = neighbors;
    	}
    	else
    	{
    		vector<Link> neighbors = adjacency[a];
    		neighbors.push_back(l);
    		adjacency[a] = neighbors;

    	}
        if(revadjacency.find(c) == adjacency.end())
        {
            vector<Link> neighbors;
            neighbors.push_back(l);
            revadjacency[c] = neighbors;
        }
        else
        {
            vector<Link> neighbors = revadjacency[c];
            neighbors.push_back(l);
            revadjacency[c] = neighbors;
        }
    	linkid++;
    }
    for(map<string,int> :: iterator it = contig2length.begin();it != contig2length.end(); ++it)
    {
        contig2degree[it->first] = get_degree(it->first);
    }
    //assign orientation to any node
    int maxlength = -1;
    string maxnode = "";
    if(pr.exist("degree"))
    {
        for(map<string, vector<Link> > :: iterator it = adjacency.begin(); it != adjacency.end();++it)
        {
            vector<Link> neighs = it->second;
            if(int(neighs.size()) > maxlength)
            {
                maxlength = neighs.size();
                maxnode = it->first;
            }
        }
    }
    else
    {
        for(map<string,int> ::iterator it = contig2length.begin(); it != contig2length.end();++it)
        {
            string contig = it->first;
            int length = it->second;
            if(length > maxlength)
            {
                maxlength = length;
                maxnode = contig;
            }
        }
    }
    string strategy;
    if(pr.exist("degree"))
    {
        strategy = "degree";
    }
    if(pr.exist("bsize"))
    {
        strategy = "bsize";
    }
    if(pr.exist("length"))
    {
        strategy = "length";
    }
    ctg2orient[maxnode] = FOW;
    invalidatelinks(maxnode,FOW);
    bfs(maxnode,strategy);
    string nd="";
    if(strategy == "bsize" || strategy == "length")
    {
        nd =get_unoriented_node_by_length();
        //cout<<nd<<endl;
    }
    else
    {
        nd =get_unoriented_node_by_degree();
    }
    while(nd != "done")
    {        
        //cout<<nd<<endl;
        ctg2orient[nd] = FOW;
        //cout<<nd<<"\t"<<ctg2orient[nd]<<endl;
        bfs(nd,strategy);
        if(strategy == "bsize" || strategy == "length")
        {
            //cout<<nd<<endl;
            nd =get_unoriented_node_by_length();
        }
        else
        {
            nd =get_unoriented_node_by_degree();
        }
    }

    int nodecounter = 1;
    map<string, int> contig2node;
    ofile << "graph ["<<endl;
    ofile << "  directed 1"<<endl;
    map<string, int> :: iterator x;
    map<string, int> actual_repeats;
    /*
    ifstream repfile("actual_repeats");
    while(getline(repfile,line))
    {
        string contig;
        double cov;
        istringstream iss(line);
        if(!(iss >> contig >> cov))
            break;
        actual_repeats[contig] = cov;
       // cout<<contig<<endl;
    }
	*/
    for(x = ctg2orient.begin(); x != ctg2orient.end(); ++x)
    {
    	string o = (x->second == 1)?"FOW":"REV";
    	string contig = x->first;
    	ofile<< "  node ["<<endl;
    	ofile<< "   id "<<nodecounter<<endl;
    	ofile<< "   label \"" <<contig<<"\""<<endl;
    	ofile<< "   orientation \""<<o<<"\""<<endl;
        ofile<< "   length \""<<contig2length[contig]<<"\""<<endl;
        string ans = "";
    	ofile<< "  ]"<<endl;
    	contig2node[contig] = nodecounter;
    	nodecounter++; 
    }
    //cerr<<"Here";
    for(map<int, Link> :: iterator it = linkmap.begin(); it!= linkmap.end();++it)
    {
        int id = it->first;
        //cerr<<"Here";
        if(invalidlinks.find(id) == invalidlinks.end())
        {
            //cout<<"HEre1"<<endl;
            Link link = it->second; 
            //cout<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
        	ofile<<"  edge ["<<endl;
        	ofile<<"   source "<<contig2node[link.getfirstcontig()]<<endl;
        	ofile<<"   target "<<contig2node[link.getsecondcontig()]<<endl;
        	ofile<<"   orientation "<<link.getlinkorientation()<<endl;
		/*
            string x = link.getfirstcontig() +"$"+link.getsecondcontig();
            if (edge2cov.find(x) == edge2cov.end())
            {
                ofile<<"   label "<<"NIL"<<endl;
            }
            else
            {
                ofile<<"   label \""<<edge2cov[x]<<"\""<<endl;
            }
            */
        	ofile<<"   mean \""<<link.getmean()<<"\""<<endl;
        	ofile<<"   stdev "<<link.getstdev()<<endl;
        	ofile<<"  ]"<<endl;
        }
    }
    ofile<<"]"<<endl;
    return 0;
}
