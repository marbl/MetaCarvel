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

class SAMRecord
{
public:
	int id;
	string read_name;
	string contig_name;
	int alignment_flag;
	int alignment_length;
	int insert_size;
	int alignment_start;
	SAMRecord () {};
	SAMRecord(int id, string read_name, string contig_name,int alignment_flag, int alignment_length, int insert_size, int alignment_start);
	bool isFirst_in_Pair();
	bool isSecond_in_Pair();
	bool is_Fow_Strand();
	bool is_Rev_Strand();
};

SAMRecord :: SAMRecord (int id, string read_name, string contig_name,int alignment_flag, int alignment_length, int insert_size, int alignment_start)
{
	this->id = id;
	this-> read_name = read_name;
	this->contig_name = contig_name;
	this->alignment_flag = alignment_flag;
	this->alignment_length = alignment_length;
	this->insert_size = insert_size;
	this->alignment_start = alignment_start;
}

bool SAMRecord :: isFirst_in_Pair()
{
	return this->alignment_flag & 64;
}

bool SAMRecord :: isSecond_in_Pair()
{
	return this->alignment_flag & 128;
}

bool SAMRecord :: is_Fow_Strand()
{
	return !(this->alignment_flag & 16);
}

bool SAMRecord :: is_Rev_Strand()
{
	return this->alignment_flag & 16;
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
map<string, SAMRecord> first_read_to_record;
map<string, SAMRecord> second_read_to_record;


char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
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

vector<SAMRecord> parseSAM(string file)
{
	vector<SAMRecord> ret;
	ifstream samfile(getCharExpr(file));
	int rec_id = 1;
	string line;
	while(getline(samfile,line))
	{
		string first_two = line.substr(0,3);
		//cout<<first_two<<endl;
		if(first_two[0] == '@')
		{
			if(first_two == "@SQ")
			{
				istringstream iss(line);
				string contigname, contiglength, buf;
				iss >> buf >> contigname >> contiglength;
				int splitindex = contigname.find(":");
				contigname = contigname.substr(splitindex+1);
				splitindex = contiglength.find(":");
				contiglength = contiglength.substr(splitindex+1);
				int len = atoi(getCharExpr(contiglength));
				contig2length[contigname] = len;
				//cout<<contigname<<"\t"<<len<<endl;
			}
		}
		else
		{
			string read_name,flag, contigname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual;
			istringstream iss(line);
			iss >> read_name >> flag >> contigname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
			int alignment_start = atoi(getCharExpr(pos));
			int alignment_length = int(seq.length());
			int insert_size = abs(atoi(getCharExpr(tlen)));
			int alignment_flag = atoi(getCharExpr(flag));
			if(contigname != "*")
			{
				//cout<<rec_id<<"\t"<<read_name<<"\t"<<contigname<<"\t"<<flag<<"\t"<<alignment_length<<"\t"<<insert_size<<"\t"<<alignment_start<<endl;
				rec_id += 1;
				if(contig2bases.find(contigname) != contig2bases.end())
				{
					contig2bases[contigname] += alignment_length;
				}
				else
				{
					contig2bases[contigname] = alignment_length;
				}
				SAMRecord record(rec_id,read_name,contigname,alignment_flag,alignment_length,insert_size,alignment_start);
				if(record.isFirst_in_Pair())
				{
					first_read_to_record[read_name] = record;
				}
				if(record.isSecond_in_Pair())
				{
					second_read_to_record[read_name] = record;
				}
				ret.push_back(record);
			}

		}
	}
	return ret;
}


int main(int argc, char* argv[])
{
	cmdline ::parser pr;
    pr.add<string>("lib_info",'l',"file containing information about library",true,"");
    pr.add<string>("alignment_info",'a',"alignment of read to assembled contigs",true,"");
    pr.add<string>("output",'o',"output file",true,"");
    pr.add<string>("coverage",'c',"contig coverage",true,"");
    pr.parse_check(argc,argv);

    //ifstream linkfile(getCharExpr(pr.get<string>("lib_info")));
    //getFastqSequences(getCharExpr(pr.get<string>("lib_info")));
	//vector<SAMRecord> alignments = parseSAM(pr.get<string>("lib_info"));
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
	vector<int> insert_sizes;
	vector<SAMRecord> alignments = parseSAM(pr.get<string>("alignment_info"));
	//iterate through all records and estimate library size
	
	for(int i = 0;i < alignments.size();i++)
	{
		SAMRecord cur_record = alignments[i];
		string cur_read = cur_record.read_name;
		//check if the read is first in pair and its mate is aligned too
		if(cur_record.isFirst_in_Pair() && first_read_to_record.find(cur_read) != first_read_to_record.end() && second_read_to_record.find(cur_read) !=second_read_to_record.end())
		{
			SAMRecord firstrecord = first_read_to_record[cur_read];
			SAMRecord secondrecord = second_read_to_record[cur_read];
			if(firstrecord.contig_name == secondrecord.contig_name)
			{
				//ask mihai about from which side we should take coordinates
				insert_sizes.push_back(abs(secondrecord.insert_size));
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

	map<string, SAMRecord> :: iterator it;
	ofstream ofile(getCharExpr(pr.get<string>("output")));
	for(it = first_read_to_record.begin(); it != first_read_to_record.end(); ++it)
	{
		string first_read = it->first;
		if(second_read_to_record.find(first_read) != second_read_to_record.end())
		{
			SAMRecord secondrecord = second_read_to_record[first_read];
			SAMRecord firstrecord = it->second;
			//check if contig for both is different
			if(firstrecord.contig_name != secondrecord.contig_name)
			{
				char firstcontigend, secondcontigend;

				bool is_first_fow = firstrecord.is_Fow_Strand();
				bool is_second_fow = secondrecord.is_Fow_Strand();
				if(is_first_fow && is_second_fow)
				{
					firstcontigend = 'E';
					secondcontigend = 'E';
				}
				if(is_first_fow && !is_second_fow)
				{
					firstcontigend = 'E';
					secondcontigend = 'B';
				}
				if(!is_first_fow && is_second_fow)
				{
					firstcontigend = 'B';
					secondcontigend = 'E';
				}
				if(!is_first_fow && !is_second_fow)
				{
					firstcontigend = 'B';
					secondcontigend = 'B';
				}
				double dist = mean - contig2length[firstrecord.contig_name] + firstrecord.alignment_start - secondrecord.alignment_start - secondrecord.alignment_length;
				
				ofile << firstrecord.contig_name<<"\t"<<firstcontigend<<"\t"<<secondrecord.contig_name<<"\t"<<secondcontigend<<"\t"<<dist<<"\t"<<stdev<<endl;
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
