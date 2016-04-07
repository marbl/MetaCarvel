#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <queue>
#include <numeric>

#include "cmdline/cmdline.h"

using namespace std;

class BedRecord
{
public:
	string contig;
	int start;
	int end;
	char strand;//+ forward - reverse
	BedRecord () {}
	BedRecord(string contig, int start, int end, char strand);

};

BedRecord :: BedRecord(string contig, int start, int end, char strand)
{
	this->contig = contig;
	this->start = start;
	this->end = end;
	this->strand = strand;
}

//change readnames to ids
map<string,BedRecord> first_in_pair;
map<string, BedRecord> second_in_pair;

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

void parse_bed(string path)
{
	ifstream bedfile(getCharExpr(path));
	string line;
	while(getline(bedfile,line))
	{	
		string contig, read;
		char strand;
		int start,end,flag;
		istringstream iss(line);
		iss >> contig >> start >> end >> read >> flag >> strand;
		BedRecord rec(contig,start,end,strand);
		if(read[read.length() -1 ] == '1')
		{
			first_in_pair[read.substr(0,read.length()-2)] = rec;
		}
		else
		{
			second_in_pair[read.substr(0,read.length()-2)] = rec;
		}

	}
}


class LibRecord
{
public:
	string lib_id;
	string read_1;
	string read_2;
	string format;
	double mean;
	double stdev;
	double maximum;
	double minimum;
	string orientation;
	LibRecord() {}
	LibRecord(string lib_id, string read_1, string read_2, string format, double mean, double stdev,double maximum, double minimum, string orientation);
	
};

LibRecord :: LibRecord(string lib_id, string read_1, string read_2, string format, double mean, double stdev,double maximum, double minimum, string orientation)
{
	this->lib_id = lib_id;
	this->read_1 = read_1;
	this->read_2 = read_2;
	this->format = format;
	this->mean = mean;
	this->stdev = stdev;
	this->maximum = maximum;
	this->minimum = minimum;
	this->orientation = orientation;
}

map<string, int> contig2length;
map<string, int> contig2bases;


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



map<string,string> getFastqSequences(string file)
{
	map<string, string> ret;
	ifstream fastqfile(getCharExpr(file));
	string line,seqname,seq;
	bool prevlineseqname = false;
	while(getline(fastqfile,line))
	{
		if(line[0] == '@')
		{
			seqname = line.substr(1);
			int space = int(seqname.find(" "));
			seqname = seqname.substr(0,space);
			prevlineseqname = true;
			continue;
		}
		if(prevlineseqname == true)
		{
			seq = line;
			ret[seqname] = seq;
			prevlineseqname = false;
		}
	}
	return ret;
}

int get_insert_size(int start1, int end1, int start2, int end2)
{
	if(start1 <= start2)
	{
		return end2 - start1 + 1;		
	}
	else
	{
		return end1 - start2 + 1;
	}
}

double estimate_distance(double mean, int start1, int end1, int start2, int end2, int ctg1_length, int ctg2_length, string orientation)
{
	int read1_length = end1 - start1 + 1;
	int read2_length = end2 - start2 + 1;
	int offset1,offset2;
	
	if(orientation == "EB")
	{
		offset1 = ctg1_length - end1;
		offset2 = start2;
	}
	//Need to work out BB and EE properly with reasoning
	if(orientation == "BB")
	{
		offset1 = start1;
		offset2 = start2;
	}
	if(orientation == "EE")
	{
		offset1 = ctg1_length - end1;
		offset2 = ctg2_length - end2;
	}
	if(orientation == "BE")
	{
		offset1 = start1;
		offset2 = ctg2_length - end2;
	}

	return mean - read1_length - read2_length - offset2 - offset1;
}



int main(int argc, char* argv[])
{
	cmdline ::parser pr;
    pr.add<string>("lib_info",'l',"file containing information about library",true,"");
    pr.add<string>("alignment_info",'a',"alignment of read to assembled contigs in bed format",true,"");
    pr.add<string>("contig_file",'d',"file containing length of contigs",true,"");
    pr.add<string>("output",'o',"output file",true,"");
    pr.add<string>("coverage",'c',"contig coverage",true,"");
    pr.parse_check(argc,argv);

    //ifstream linkfile(getCharExpr(pr.get<string>("lib_info")));
    //getFastqSequences(getCharExpr(pr.get<string>("lib_info")));
	//vector<SAMRecord> alignments = parseSAM(pr.get<string>("lib_info"));
	get_contig_length(pr.get<string>("contig_file"));
	ifstream libfile(getCharExpr(pr.get<string>("lib_info")));
	vector<LibRecord> libraries;
	string line;
	while(getline(libfile,line))
	{
		istringstream iss(line);
		string a,b,c,d,e;
		double mean, stdev, minimum, maximum;
		iss >> a >> b >> c >> d >> mean >> stdev >> minimum >> maximum >> e;
		LibRecord record(a,b,c,d,mean,stdev,minimum,maximum,e);
		libraries.push_back(record);
	}
	parse_bed(pr.get<string>("alignment_info"));
	vector<int> insert_sizes;
	//iterate through all records and estimate library size
	
	map<string,BedRecord> :: iterator it;
	for(it = first_in_pair.begin(); it != first_in_pair.end();++it)
	{
		string read = it->first;
		BedRecord first = it->second;
		if(second_in_pair.find(read) != second_in_pair.end())
		{
			BedRecord second = second_in_pair[read];
			if(first.contig == second.contig)
			{
				int insert_size = get_insert_size(first.start, first.end, second.start, second.end);
				//cout<<insert_size<<endl;
				insert_sizes.push_back(insert_size);
			}
		}
	}
	

	double sum = std::accumulate(insert_sizes.begin(), insert_sizes.end(), 0.0);
	double mean = sum / insert_sizes.size();

	double sq_sum = std::inner_product(insert_sizes.begin(), insert_sizes.end(), insert_sizes.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / insert_sizes.size() - mean * mean);
	cerr<<"Mean = "<<mean<<endl;
	cerr<<"Stdev = "<<stdev<<endl;
	//calculate links between contigs based on mate pair information, iterate through maps of mate pairs and find links
	ofstream ofile(getCharExpr(pr.get<string>("output")));
	for(it = first_in_pair.begin(); it != first_in_pair.end(); ++it)
	{
		BedRecord first = it->second;
		string firstcontigend, secondcontigend;
		if(second_in_pair.find(it->first) != second_in_pair.end())
		{
			BedRecord second = second_in_pair[it->first];
			if(first.contig != second.contig)
			{
				if(first.strand == '+' && second.strand == '+')
				{
					firstcontigend = "E";
					secondcontigend = "E";
				}
				if(first.strand == '+' && second.strand == '-')
				{
					firstcontigend = "E";
					secondcontigend = "B";
				}
				if(first.strand == '-' && second.strand == '+')
				{
					firstcontigend = "B";
					secondcontigend = "E";
				}
				if(first.strand == '-' && second.strand == '-')
				{
					firstcontigend = "B";
					secondcontigend = "B";
				}
				double dist = estimate_distance(mean,first.start,first.end,second.start,second.end,contig2length[first.contig],contig2length[second.contig],firstcontigend+secondcontigend);
					
				ofile << first.contig<<"\t"<<firstcontigend<<"\t"<<second.contig<<"\t"<<secondcontigend<<"\t"<<dist<<"\t"<<stdev<<endl;

			}
		}
	}

	ofstream cfile(getCharExpr(pr.get<string>("coverage")));
	map<string,int> :: iterator it2;
	for(it2 = contig2bases.begin(); it2 != contig2bases.end();++it2)
	{
		string contig = it2->first;
		int contiglen = contig2length[contig];
		double coverage = it2->second*1.0/contiglen;
		cfile<<contig<<"\t"<<coverage<<endl;
	}
	return 0;
}
