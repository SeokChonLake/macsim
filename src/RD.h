/************************************************************
 * File		: RD.h
 * Author	:
 * Data		: 
 * SVN		:
 * Description	:
 * ***********************************************************/

#ifndef RD_H
#define RD_H

#include <string>
#include <vector>
#include "cache.h"
/*
class RD_sampler{
	public:
		Addr tag;
		int set;
		int count;
};

vector<RD_sampler> reuse;  // save access for calculate reuse distance
int RDD[MAX_RDD_SIZE]={0};  // memorize reuse distance distribution initialize all 0 
int access_count=0;
*/
class RD_cache_c: public cache_c
{
	public:


		RD_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
			bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
			int num_tiles, int interleave_factor, macsim_c * simBase);
		virtual ~RD_cache_c();
		void update_cache_on_access(Addr tag, int set, int appl_id);
	
	
/*	
		class RD_sampler{
			public:
				Addr tag;
				int set;
				int count;
		};

		vector<RD_sampler> reuse;  // save access for calculate reuse distance
		int RDD[MAX_RDD_SIZE]={0};  // memorize reuse distance distribution initialize all 0 
		int access_count=0;
*/

	protected:
		int	m_num_cores;
		int isProtected;
		int PD;
};

#endif
