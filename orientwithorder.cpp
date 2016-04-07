#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_map>

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

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

//add bundle size



map<string, vector<Link> > adjacency;
map<string, vector<Link> > revadjacency;
map<string, int> ctg2orient;
map<int, bool> invalidlinks;
map<string, int> contig2length;
map<string, int> contigs2bundle;
map<string, int> contig2degree;
unordered_map<int,bool> assigned;
unordered_map<string,int> contigs2links;
unordered_map<string,int> contig2id;
unordered_map<int,string> id2contig;

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
        cout<<"FOW = "<<curr_fow<<"   REV = "<<curr_rev<<endl;

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
    //cerr<<"invalidating..."<<v<<endl;
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
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BE")
                    {
                        invalidlinks[link.getid()] = true;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == FOW && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "EE")
                    {
                        invalidlinks[link.getid()] = true;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "BB")
                    {
                        invalidlinks[link.getid()] = true;
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
                        //cout<<link.getlinkorientation()<<"\t"<<"BE"<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BE")
                    {
                        invalidlinks[link.getid()] = true;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == REV && neighorientation == FOW)
                {
                    if(link.getlinkorientation() != "EE")
                    {
                        invalidlinks[link.getid()] = true;
                        //cout<<orientation<<"\t"<<neighorientation<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
                if(orientation == FOW && neighorientation == REV)
                {
                    if(link.getlinkorientation() != "BB")
                    {
                        invalidlinks[link.getid()] = true;
                        //cout<<orientation<<"\t"<<neighorientation<<endl;
                        //cerr<<link.getfirstcontig()<<"\t"<<link.getfirstorietation()<<"\t"<<link.getsecondcontig()<<"\t"<<link.getsecondorientation()<<endl;
                    }
                }
            }
        }
    }
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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

/*helper to find component*/

set<int> bfs(int start,unordered_map<int,vector<int> > undirected_graph)
{   
    set<int> component;
    component.insert(start);
    assigned[start] = true;
    queue<int> Q;
    Q.push(start);
    while(!Q.empty())
    {
        int n = Q.front();
        //cout<<Q.size()<<endl;
        Q.pop();
        for(int i = 0;i < undirected_graph[n].size();i++)
        {
            int neighbor = undirected_graph[n][i];
            if(component.find(neighbor) == component.end())
            {
                component.insert(neighbor);
                Q.push(neighbor);
                assigned[neighbor] = true;
            }
        }
    }
    return component;
}

/*
    Gives connected components as a vector of sets using adjacency list
*/
vector<set<int> > getConnectedComponends()
{
    vector<set<int> > components;
    map<string, vector<Link> > :: iterator it;
    unordered_map<int,vector<int> > undirected_graph;

    for(it = adjacency.begin(); it != adjacency.end(); ++it)
    {
        vector<Link> neighbors = it->second;
        if(undirected_graph.find(contig2id[it->first]) == undirected_graph.end())
        {
            vector<int> neighs;
            undirected_graph[contig2id[it->first]] = neighs;
        }
        for(vector<Link> :: iterator it1 = neighbors.begin(); it1 != neighbors.end();++it1)
        {
            Link l = *it1;
            undirected_graph[contig2id[it->first]].push_back(contig2id[l.getsecondcontig()]);
            if(undirected_graph.find(contig2id[l.getsecondcontig()]) == undirected_graph.end())
            {
                vector<int> neighs;
                undirected_graph[contig2id[l.getsecondcontig()]] = neighs;
            }
            undirected_graph[contig2id[l.getsecondcontig()]].push_back(contig2id[it->first]);
        }
    }
    unordered_map<int, vector<int> > :: iterator it2;
    for(it2 = undirected_graph.begin(); it2 != undirected_graph.end(); ++it2)
    {
        if(assigned.find(it2->first) == assigned.end())
        {
            components.push_back(bfs(it2->first,undirected_graph));
        }
    }
    return components;
}

/*
Function that takes a component and returns a weighted graph based on edge weight function.
*/
unordered_map<int,vector<pair<int,int> > > getWeightedGraph(set<int> component)
{
    //cerr<<"in the method"<<endl;
    unordered_map<int,vector<pair<int,int> > > weighted_graph;
    map<string,int> forward_scores;
    map<string,int> reverse_scores;  
    //cerr<<component.size()<<endl;     
    for(set<int> :: iterator it = component.begin(); it != component.end(); ++it)
    {
        int nd = *it;
        string curr_node = id2contig[nd];
        //cerr<<curr_node<<endl;
        if(adjacency.find(curr_node) != adjacency.end())
        {
            vector<Link> adj = adjacency[curr_node];
            for(int i = 0; i < adj.size();i++)
            {
                Link link = adj[i];
                string adj_contig = link.getsecondcontig();
                if(link.getlinkorientation() == "BE" || link.getlinkorientation() == "EB")
                {
                    string key;
                    if(curr_node <= adj_contig)
                    {
                        key = curr_node +"$"+adj_contig;
                    }
                    else
                    {
                        key = adj_contig + "$" + curr_node;
                    }

                    if(forward_scores.find(key) == forward_scores.end())
                    {
                        forward_scores[key] = 0;
                    }
                    forward_scores[key] += link.get_bundle_size();
                }
                if(link.getlinkorientation() == "BB" || link.getlinkorientation() == "EE")
                {
                    string key;
                    if(curr_node <= adj_contig)
                    {
                        key = curr_node +"$"+adj_contig;
                    }
                    else
                    {
                        key = adj_contig + "$" + curr_node;
                    }

                    if(reverse_scores.find(key) == reverse_scores.end())
                    {
                        reverse_scores[key] = 0;
                    }
                    reverse_scores[key] += link.get_bundle_size();
                }
            }
        }
    }
    //cerr<<"done"<<endl;
    //cerr<<"size of fow scores "<<forward_scores.size()<<endl;
    //cerr<<"size of rev scores "<<reverse_scores.size()<<endl;
    //finds edge weights for just forward orientations and forward+reverse orientations
    for(map<string,int> :: iterator it = forward_scores.begin(); it != forward_scores.end();++it)
    {
        string key = it->first;
        vector<string> contigs = split(key,'$');
        //cerr<<contigs[0]<<"\t"<<contigs[1]<<endl;
        int fow_weight = it->second;
        int edge_weight;
        int node0 = contig2id[contigs[0]];
        int node1 = contig2id[contigs[1]];
        if(weighted_graph.find(node0) == weighted_graph.end())
        {
            vector<pair<int,int> > neighbors;
            weighted_graph[node0] = neighbors;
        }
        if(weighted_graph.find(node1) == weighted_graph.end())
        {
            vector<pair<int,int> > neighbors;
            weighted_graph[node1] = neighbors;
        }
        if(reverse_scores.find(key) != reverse_scores.end())
        {
            int rev_weight = reverse_scores[key];
            edge_weight = abs(fow_weight - rev_weight);      
        }
        else
        {
            edge_weight = fow_weight;
        }
        //cerr<<"Edge weight = "<<edge_weight<<endl;
        pair<int,int> first = make_pair(node0,edge_weight);
        pair<int,int> second = make_pair(node1,edge_weight);
        weighted_graph[node1].push_back(first);
        weighted_graph[node0].push_back(second);
        //cerr<<"first done"<<endl;
    }
    //find edge weights for just reverse orientations
    for(map<string,int> :: iterator it = reverse_scores.begin(); it != reverse_scores.end();++it)
    {
        string key = it->first;
        vector<string> contigs = split(key,'$');
        int node0 = contig2id[contigs[0]];
        int node1 = contig2id[contigs[1]];
        if(forward_scores.find(key) == forward_scores.end())
        {
            int edge_weight = reverse_scores[key];
            if(weighted_graph.find(node0) == weighted_graph.end())
            {
                vector<pair<int,int> > neighbors;
                weighted_graph[node0] = neighbors;
            }
            if(weighted_graph.find(node1) == weighted_graph.end())
            {
                vector<pair<int,int> > neighbors;
                weighted_graph[node1] = neighbors;
            }
            pair<int,int> first = make_pair(node0,edge_weight);
            pair<int,int> second = make_pair(node1,edge_weight);
            weighted_graph[node1].push_back(first);
            weighted_graph[node0].push_back(second);
        }
    }
    return weighted_graph;
}


unordered_map<int,vector<int> > Prim_MST(unordered_map<int,vector<pair<int,int> > > graph)
{
    int start = graph.begin()->first;
    unordered_map<int,int> parent;
    unordered_map<int,int> distance;
    set<pair<int,int > > Q;
    for(unordered_map<int,vector<pair<int,int> > > :: iterator it = graph.begin(); it != graph.end();++it)
    {
        distance[it->first] = 1<<28;
        Q.insert(make_pair(distance[it->first],it->first));
    }
    Q.erase(Q.find(make_pair(distance[start],start)));
    distance[start] = 0;
    Q.insert(make_pair(distance[start],start));
    parent[start] = -1;
    while(!Q.empty())
    {
        pair<int,int> top = *Q.begin();
        //cout<<Q.size()<<endl;
        Q.erase(Q.begin());
        int u = top.second;
        vector<pair<int,int> >  neighbors = graph[u];
        for(int i = 0;i < neighbors.size();i++)
        {
            int v = neighbors[i].first;
            int weight = neighbors[i].second;
            //cerr<<"v = "<<v<<"; weight = "<<weight<<"; distance[v] = "<<distance[v]<<endl;
            if(distance[v] > weight && Q.find(make_pair(distance[v],v)) != Q.end())
            {
                parent[v] = u;
                Q.erase(Q.find(make_pair(distance[v],v)));
                distance[v] = weight;
                Q.insert(make_pair(weight,v));
            }
        }
    }
    unordered_map<int,vector<int> > ret;
    for(unordered_map<int,int> :: iterator it = parent.begin();  it != parent.end(); ++it)
    {
        //cout<<"here"<<endl;
        int u = it->first;
        int v = it->second;
        if(v == -1)//start node
            continue;
        if(ret.find(u) == ret.end())
        {
            vector<int> neighbors;
            ret[u] = neighbors;
        }
        if(ret.find(v) == ret.end())
        {
            vector<int> neighbors;
            ret[v] = neighbors;
        }
        if(v < u)
            ret[v].push_back(u);
        else
            ret[u].push_back(v);
    }
    return ret;
}

unordered_map<int,vector<int> > getMaximumSpanningTree(unordered_map<int,vector<pair<int,int> > > graph)
{
    //first negate all edge weights
    int num_nodes = graph.size();
    for(unordered_map<int,vector<pair<int,int> > >:: iterator it = graph.begin(); it != graph.end(); ++it)
    {
        int node = it->first;
        vector<pair<int,int> > actual_edges = it->second;
        for(int i = 0;i < actual_edges.size();++i)
        {
            pair<int,int> edge = actual_edges[i];
            edge.second = -1*edge.second;
            actual_edges[i] = edge;
        }
        graph[node] = actual_edges;
    }
    //cerr<<"Graph negated"<<endl;
    return Prim_MST(graph);
}


pair<int,int> first_dfs(unordered_map<int,vector<int> > graph, int root, unordered_map<int,pair<int,int> > &values)
{
    //check if root is leaf
    if(graph[root].size() == 0)
    {
        pair<int,int> val = make_pair(0,1);
        values[root] = val;
        return val;
    }
    //find (sum of path and subtree size for each node)
    int root_distance = 0;
    int root_treesize = 0;
    for(int i = 0; i < graph[root].size();i++)
    {
        int neighbor = graph[root][i];
        pair<int,int> new_values = first_dfs(graph,neighbor,values);
        root_distance += (new_values.first + new_values.second);
        root_treesize += (new_values.second);
    }
    //increment treesize by 1 to count for root
    root_treesize++;
    values[root] = make_pair(root_distance,root_treesize);
    return values[root];
}

void second_dfs(unordered_map<int,vector<int> > graph, int root, unordered_map<int,pair<int,int> > &values,unordered_map<int,int> &total_distances)
{
    int total_nodes = values.size();
    for(int i = 0; i < graph[root].size();i++)
    {
        int curr_node = graph[root][i];
        int new_distance = total_distances[root] -  values[curr_node].first - values[curr_node].second+ 1 + values[curr_node].first + total_nodes - values[curr_node].second - 1;
        total_distances[curr_node] = new_distance;
        second_dfs(graph,curr_node,values,total_distances);
    }
}

void start_orientation(unordered_map<int,vector<int> > graph, int root)
{
    //cout<<id2contig[root]<<endl;
    for(int i = 0;i < graph[root].size(); i++)
    {
        int v = graph[root][i];
        if(ctg2orient[id2contig[v]] == NIL)
        {
            int orientation = findorientation(id2contig[v]);
            ctg2orient[id2contig[v]] = orientation;
            invalidatelinks(id2contig[v],orientation);
            start_orientation(graph,v);
        }
    }
}

int main(int argc, char* argv[])
{
	
    cmdline ::parser pr;
    pr.add<string>("bundled_graph",'l',"list of bundled links",true,"");
    pr.add<string>("contig_length",'c',"contig lengths",true,"");
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
    int contigid = 0;
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
        if(contig2id.find(a) == contig2id.end())
        {
            contig2id[a] = contigid;
            id2contig[contigid] = a;
            contigid++;
        }
        if(contig2id.find(c) == contig2id.end())
        {
            contig2id[c] = contigid;
            id2contig[contigid] = c;
            contigid++;
        }

        contigs2links[a+c] = linkid;
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
    //first, get all connected components in the graph
    vector<set<int> > components = getConnectedComponends();

    //now for each connected component, transform it to a weighted graph, where weight is defined in Myer's paper
    for(int i = 0; i < components.size(); i++)
    {
        set<int> component = components[i];
        unordered_map<int,vector<pair<int,int> > > weighted_graph = getWeightedGraph(component);
        unordered_map<int,vector<int> > mst = getMaximumSpanningTree(weighted_graph);
        unordered_map<int,pair<int,int> > values;
        first_dfs(mst,mst.begin()->first,values);
        unordered_map<int,int> distances;
        second_dfs(mst,mst.begin()->first,values,distances);
        int maxnode,maxdist = -1;
        for(unordered_map<int,int> :: iterator it = distances.begin(); it != distances.end();++it)
        {
            if(it->second > maxdist)
            {
                maxdist = it->second;
                maxnode = it->first;
            }
        }
        ctg2orient[id2contig[maxnode]] = FOW;
        unordered_map<int,vector<int> > mst_all;
        for(unordered_map<int,vector<int> > :: iterator it = mst.begin(); it != mst.end(); ++it)
        {
            if(mst_all.find(it->first) == mst_all.end())
            {
                vector<int> neigh;
                mst_all[it->first] = neigh;
            }
            vector<int> neigh = it->second;
            for(int i = 0; i < neigh.size();i++)
            {
                int v = neigh[i];
                mst_all[it->first].push_back(v);
                if(mst_all.find(v) == mst_all.end())
                {
                    vector<int> neigh;
                    mst_all[v] = neigh;
                }
                mst_all[v].push_back(it->first);
            }
        }
        start_orientation(mst_all,maxnode);
    }
        //cout<<"done"<<endl;
        /*
        for(unordered_map<int,vector<pair<int,int> > > :: iterator it = weighted_graph.begin(); it != weighted_graph.end();++it)
        {
            cout<<id2contig[it->first]<<endl;
            vector<pair<int, int> > neighbors = it->second;
            for(int i = 0;i < neighbors.size();i++)
            {
                pair<int,int> neigh = neighbors[i];
                cout<<"\t("<<id2contig[neigh.first]<<","<<neigh.second<<")"<<endl;
            }
        }
        cout<<"==========="<<endl;
        
        unordered_map<int,vector<int> > mst = getMaximumSpanningTree(weighted_graph);
        cerr<<"done"<<endl;
        for(unordered_map<int,vector<int> > :: iterator it = mst.begin(); it != mst.end(); ++it)
        {
            cout<<id2contig[it->first]<<"(";
            vector<int> neighs = it->second;
            for(int i = 0; i < neighs.size();i++)
            {
                cout<<"\t"<<id2contig[neighs[i]];
            }
            cout<<")"<<endl;
        }
    }*/
    /*
    //To print each connected component
    for(int i = 0;i < components.size();i++)
    {
        set<string> component = components[i];
        for(set<string> :: iterator it = component.begin(); it != component.end(); ++it)
        {
            cout<<*it<<"\t";
        }
        cout<<endl;
    }*/
    
    int nodecounter = 1;
    map<string, int> contig2node;
    ofile << "graph ["<<endl;
    ofile << "  directed 1"<<endl;
    map<string, int> :: iterator x;
    map<string, int> actual_repeats;
    
    for(x = ctg2orient.begin(); x != ctg2orient.end(); ++x)
    {
    	string o = (x->second == 1)?"FOW":"REV";
    	string contig = x->first;
    	ofile<< "  node ["<<endl;
    	ofile<< "   id "<<nodecounter<<endl;
    	ofile<< "   label \"" <<contig<<"\""<<endl;
    	ofile<< "   orientation \""<<o<<"\""<<endl;
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
        	ofile<<"   mean \""<<link.getmean()<<"\""<<endl;
        	ofile<<"   stdev "<<link.getstdev()<<endl;
        	ofile<<"  ]"<<endl;
        }
    }
    ofile<<"]"<<endl;
    return 0;
}
