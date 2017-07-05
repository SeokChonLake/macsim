#include "assert_macros.h"
#include "RD.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

#define MAX_RDD_SIZE 256
#define SET_SIZE 2048


class RD_sampler{
	public:
		Addr tag;
		int set;
		int count;
		int appl_id;
};

vector<RD_sampler> reuse;

long long int RDD[10][MAX_RDD_SIZE]={0};
int access_count[SET_SIZE]={0};
int APPL_NUM;
RD_cache_c::RD_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
	bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
	int num_tiles, int interleave_factor, macsim_c* simBase)
: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
	cache_by_pass, core_id, cache_type_info, enable_partition,
	num_tiles, interleave_factor, simBase){
	for(int i=0;i<SET_SIZE;i++)
		access_count[i]=0;

	APPL_NUM = *KNOB(KNOB_NUM_SIM_CORES);
		
}
RD_cache_c::~RD_cache_c()
{
	ofstream File;
	int i;
	File.open("RDD.csv",ios_base::out);
	for(i=0;i<MAX_RDD_SIZE;i++){
		File << i+1;

		for(int j=0;j<APPL_NUM;j++)
			File <<","<< RDD[j][i];
		
		File << endl;

	}
	File.close();




	while(!reuse.empty())
		reuse.pop_back();
}


void RD_cache_c::update_cache_on_access(Addr tag,int set,int appl_id){
	int i,j;
	RD_sampler tmp;
	int find;
//	cout<<"<"<<set<<">";
	
	access_count[set]++;
	find = -1;
	for(i=0;i<reuse.size();i++){
		tmp = reuse[i];
		if(tmp.tag==tag && tmp.set == set&& tmp.appl_id==appl_id)
		{
			find = i;	
			break;
		}
	}
	//not found
	if(find==-1){ // first access
		tmp.tag=tag;
		tmp.set=set;
		tmp.count=access_count[set];
		tmp.appl_id=appl_id;
		reuse.push_back(tmp);
	}
	//found
	else{  // update all
		/*
		for(i = 0; i < reuse.size();i++){
			if(tmp.count>reuse[i].count)
				if(reuse[i].set==set)
					reuse[i].count++;
		}
		*/
		
		if(access_count[set]-reuse[find].count-1>=MAX_RDD_SIZE)
			RDD[appl_id][MAX_RDD_SIZE-1]++;
		else
			RDD[appl_id][access_count[set]-reuse[find].count-1]++;
		reuse[find].count=access_count[set];
		
	}


}

