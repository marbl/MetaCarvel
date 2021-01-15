#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>

#include "cmdline/cmdline.h"

using namespace std;

class Link
{
private:
	int id;
	string contig_a;
	string contig_a_orientation;
	string contig_b;
	string contig_b_orientation;
    int bundle_size;
	double mean;
	double stdev;
public:
	Link() {};
	Link(int id, string contig_a, string contig_a_orientation, string contig_b, string contig_b_orientation, double mean, double stdev);
	double getmean();
	double getstdev();
	string getlinkorientation();
    string getcontigs();
    string getfirstcontig();
    string getsecondcontig();
    string getfirstorietation();
    string getsecondorientation();
    void set_bundle_size(int size);
    int get_bundle_size();
    int getid();
};

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

void Link :: set_bundle_size(int size)
{
    this->bundle_size = size;
}

int Link :: get_bundle_size()
{
    return this->bundle_size;
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


int main(int argc, char* argv[])
{
    cmdline ::parser pr;
    pr.add<string>("contigs",'l',"contig links",true,"");
    pr.add<string>("output",'o',"output file",true,"");
    pr.add<string>("bgraph",'b',"bundled graph in gml format",true,"");
    pr.add<int>("cutoff",'c',"number of mate pairs to support an edge",false,3);
    pr.parse_check(argc,argv);

    ifstream linkfile(getCharExpr(pr.get<string>("contigs")));
    ofstream ofile(getCharExpr(pr.get<string>("output")));
    ofstream g(getCharExpr(pr.get<string>("bgraph")));

    string line;
    int linkid = 1;
    int cutoff = pr.get<int>("cutoff");
    map<int, Link> linkmap;

    while(getline(linkfile,line))
    {
    	string a,b,c,d;
    	double e,f;
    	istringstream iss(line);
    	if(!(iss >> a >> b >> c >> d >> e >> f))
    		break;
    	Link l(linkid,a,b,c,d,e,f);
    	linkmap[linkid] = l;
    	linkid++;
    }
    
    //Store links for a pair of contigs and orientation. For each possible pair, there can be 4 orientations

    map<int, Link> :: iterator it;
    //first map: contig-> second_map, second_map: orientation->links
    map<string, map<string, vector<Link> > > contig_to_links;
    for(it = linkmap.begin(); it!= linkmap.end(); ++it)
    {
        Link link = it->second;
        string contigs = link.getcontigs();
        string orientations = link.getlinkorientation();
        if(contig_to_links.find(contigs) == contig_to_links.end())
        {
            map<string, vector<Link> > orientation_to_links;
            if(orientation_to_links.find(orientations) == orientation_to_links.end())
            {
                vector<Link> links;
                links.push_back(link);
                orientation_to_links[orientations] = links;
            }
            contig_to_links[contigs] = orientation_to_links;
        }
        else
        {
            map<string, vector<Link> > orientation_to_links = contig_to_links[contigs];
            if(orientation_to_links.find(orientations) == orientation_to_links.end())
            {
                vector<Link> links;
                links.push_back(link);
                orientation_to_links[orientations] = links;
            }
            else
            {
                vector<Link> links = orientation_to_links[orientations];
                links.push_back(link);
                orientation_to_links[orientations] = links;
            }
            contig_to_links[contigs] = orientation_to_links;
        }
    }
    //cerr<<"Links loaded and stored by orientation"<<endl;

    
    //code to test if link are properly loaded
    for(map<string, map<string, vector<Link> > > :: iterator it = contig_to_links.begin(); it != contig_to_links.end(); ++it)
    {
        //cout<<it->first<<endl;
        map<string, vector<Link> >  olinks = it->second;
        for(map<string, vector<Link> > :: iterator it1 = olinks.begin(); it1 != olinks.end(); ++it1)
        {
            //cout<<"\t"<<it1->first<<endl;
            vector<Link> l = it1->second;
            for(int i = 0; i < l.size();i++)
            {   
                Link link = l[i];
                //cout<<"\t\t"<<link.getcontigs()<<endl;
            }
        }
    }
    

    //For each pair of contig, for each possible orientation apply maximal clique algorithm and compress links to 1

    vector<Link> bundled_links;
    map<string, map<string, vector<Link> > > :: iterator linkit;
    for(linkit = contig_to_links.begin(); linkit != contig_to_links.end();++linkit)
    {
        map<string, vector<Link> > orientation_to_links = linkit->second;
        map<string, vector<Link> > :: iterator it1;
        for(it1 = orientation_to_links.begin(); it1 != orientation_to_links.end();it1++)
        {
            string orientation = it1->first;
            vector<Link> links = it1->second;
            //Apply clique algorithm only if number of link with same orientation is more than cutoff
            if(links.size() > cutoff)
            {
                vector< pair<int,double> > begins;
                vector< pair<int,double> > ends;

                for(int i = 0;i < links.size();i++)
                {
                    Link link = links[i];
                    double mean = link.getmean();
                    double stdev = link.getstdev();
                    begins.push_back(make_pair(link.getid(), mean - 3*stdev));
                    ends.push_back(make_pair(link.getid(),mean + 3* stdev));
                }

                //sort begins and ends in increasing order
                sort(begins.begin(),begins.end(),pairCompare);
                sort(ends.begin(),ends.end(),pairCompare);
                int start_index = 0;
                int end_index = 0;
                int curr_clique = 0, best_clique = 0;
                double best_coord = -100000;
                vector<Link> clique_links;
                double begin_left, begin_right, end_left, end_right;
                int printed = 0;
                while(start_index < begins.size() && end_index < ends.size())
                {
                    if(start_index < begins.size() - 1 && begins[start_index].second <= ends[end_index].second)
                    {
                        int linkno = begins[start_index].first;
                        Link curlink = linkmap[linkno];
                        begin_left = curlink.getmean() - 3*curlink.getstdev();
                        begin_right = curlink.getmean() + 3*curlink.getstdev();
                        //cout<<begin_right<<endl;
                        curr_clique++;
                        //cout<<"here"<<endl;
                        if (curr_clique > best_clique)
                        {
                            best_clique = curr_clique;
                            clique_links.clear();
                            best_coord = begin_left;
                        }
                        start_index++;
                    }
                    else
                    {   
                        if((end_index < ends.size()) && ((start_index == begins.size() - 1 || (begins[start_index].second > ends[end_index].second))))
                        {
                            //cout<<"here"<<endl;
                            int linkno = ends[end_index].first;
                            Link curlink = linkmap[linkno];
                            end_left = curlink.getmean() - 3*curlink.getstdev();
                            end_right = curlink.getmean() + 3*curlink.getstdev();

                            if(end_left <= best_coord && end_right >= best_coord)
                            {
                                clique_links.push_back(linkmap[ends[end_index].first]);
                                //cout<<"adding link"<<endl;
                            }
                            curr_clique--;
                            end_index++;
                        }

                    }
                }
                //cerr<<"Clique Done"<<endl;
                //cerr<<best_clique<<endl;
                if(clique_links.size() != 0)
                {
                    double min_range = clique_links[0].getmean() - 3*clique_links[0].getstdev();
                    double max_range = clique_links[0].getmean() + 3*clique_links[0].getstdev();
                    for(int i = 1; i < clique_links.size();i++)
                    {
                        Link link = clique_links[i];
                        double mean = link.getmean();
                        double stdev = link.getstdev();
                        if(mean - 3*stdev > min_range)
                            min_range = mean - 3*stdev;
                        if(mean + 3*stdev < max_range)
                            max_range = mean + 3*stdev;
                    }

                    //write code to log invalid links

                    double newmean, newsd, p = 0,q = 0;
                    for(int i = 0;i < clique_links.size();i++)
                    {
                        Link link = clique_links[i];
                        double tmp = link.getstdev();
                        if(tmp == 0)
                            tmp = 1;
                        tmp  = tmp*tmp;
                        p += link.getmean()*1.0/tmp;
                        q += 1.0/tmp;
                    }
                    newmean = p/q;
                    newsd = 1/sqrt(q);
                    //cout<<clique_links.size()<<endl;
                    Link templink = clique_links[0];
                    Link newlink(0, templink.getfirstcontig(), templink.getfirstorietation(), templink.getsecondcontig(), templink.getsecondorientation(),
                        newmean, newsd);
                    //cout<<clique_links.size()<<endl;
                    newlink.set_bundle_size(clique_links.size());
                    bundled_links.push_back(newlink);
                }
            }
            else
            {   
            	links[0].set_bundle_size(1);
                bundled_links.push_back(links[0]);
                
            }
        }

    }
    int nodeid = 1;
    map<string,int> contig2node;
    for(int i = 0;i < bundled_links.size();i++)
    {
        Link l = bundled_links[i];
        string contiga = l.getfirstcontig();
        string contigb = l.getsecondcontig();
        if(contig2node.find(contiga) == contig2node.end())
        {
            contig2node[contiga] = nodeid;
            nodeid++;
        }
        if(contig2node.find(contigb) == contig2node.end())
        {
            contig2node[contigb] = nodeid;
            nodeid++;
        }
    }

    g <<"graph ["<<endl;
    g <<" directed 1"<<endl;
    for(map<string,int> :: iterator it = contig2node.begin(); it != contig2node.end(); ++it)
    {
        g<<" node ["<<endl;
        g<<"  id "<<it->second<<endl;
        g<<"  label \""<<it->first<<"\""<<endl;
        g<<" ]"<<endl;
    }
    for(int i = 0;i < bundled_links.size();i++)
    {
        Link l = bundled_links[i];
        if (l.get_bundle_size() >= cutoff)
        {
            g<<" edge ["<<endl;
            g<<"  source "<<contig2node[l.getfirstcontig()]<<endl;
            g<<"  target "<<contig2node[l.getsecondcontig()]<<endl;
            g<<"  mean "<<l.getmean()<<endl;
            g<<"  stdev "<<l.getstdev()<<endl;
            g<<"  bsize "<<l.get_bundle_size()<<endl;
            g<<" ]"<<endl;
        }
    }
    g<<"]";
    for(int i = 0;i < bundled_links.size();i++)
    {
        Link l = bundled_links[i];
        if (l.get_bundle_size() >= cutoff)
            ofile<<l.getfirstcontig()<<"\t"<<l.getfirstorietation()<<"\t"<<l.getsecondcontig()<<"\t"<<l.getsecondorientation()<<"\t"<<l.getmean()<<"\t"<<l.getstdev()<<"\t"<<l.get_bundle_size()<<endl;
    }
    //write code to dump to gml file
    return 0;
}
