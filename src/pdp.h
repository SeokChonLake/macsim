/*
Copyright (c) <2012>, <Georgia Institute of Technology> All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted 
provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of 
conditions and the following disclaimer in the documentation and/or other materials provided 
with the distribution.

Neither the name of the <Georgia Institue of Technology> nor the names of its contributors 
may be used to endorse or promote products derived from this software without specific prior 
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
*/


/**********************************************************************************************
 * File         : pdp.h
 * Description  : pdp cache structure (based on cache_lib at scarab) 
 *********************************************************************************************/

#ifndef PDP_H
#define PDP_H


#include <string>

#include "macsim.h"
#include "global_types.h" 
#include "global_defs.h"
#include "cache.h"

#include <vector>
#include <algorithm>

#define CPU_CORE 0
#define IP_CORE 1

struct rd_record{
        Addr tag;
        int set;
        Counter count;
        int appl_id;

        rd_record() {}
        rd_record(Addr _tag, int _set, Counter _count, int _appl_id) {
            tag = _tag; set = _set; count = _count; appl_id = _appl_id;
        }
};

///////////////////////////////////////////////////////////////////////////////////////////////
/// \brief PDP cache library class
///////////////////////////////////////////////////////////////////////////////////////////////
class pdp_cache_c: public cache_c 
{
    public:
        pdp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num, 
                    bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
                    int num_tiles, int interleave_factor, macsim_c* simBase); 

        virtual ~pdp_cache_c();

        void update_cache_on_access(Addr tag, int set, int appl_id);
        void update_line_on_hit(cache_entry_c* line, int set, int appl_id);
		void update_cache_on_miss(int set, int appl_id);
//        void update_set_on_replacement(Addr tag, int appl_id, int set_id, bool gpuline);
        cache_entry_c* find_replacement_line(int set, int appl_id);
        void initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, 
                                   int appl_id, bool gpuline, int set_id, bool skip);
        void update_cache_policy(Counter m_cycle);
	void cache_set_promotion(int set);
    /* Additional Member for PDP */
    public:
        int m_num_cores;                      /**< number of master IP */
		int m_num_set;						  /**< number of set of cache */
        int m_max_rd;                         /**< maximum reuse distance */
        int *m_pd;                            /**< per-master IP Protecting Distance */                 
        vector<rd_record> m_rd_sampler;       /**< Reuse distance sampler */
        Counter **m_rd_counter;               /**< per-master IP RD Counter Array */
		Counter **m_real_counter;	
		Counter **m_eviction_counter;
        Counter *m_access_counter;            /**< per-set Access Counter */
        bool *m_ip_table;                     /**< type master IP */
        int *m_allocation;                    /**< number of allocated ways */ 
        int *m_occupied;                      /**< number of occupied ways */
		int *pre_hit_rate;
		bool real_rd_flag; /*for real rd in cache now */
    /* Additional Method for PDP */
    public:
		void calculate_way(void);
		void calculate_pd(void);
        bool bypass_monitor(Addr addr, int m_appl_id);
        void update_reuse_distance(Addr tag, int set, int appl_id);
        cache_entry_c* is_there_a_room(int set, int appl_id);
};

#endif // PDP_H
