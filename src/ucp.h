/************************************************************
 * File		: ucp.h
 * Author	:
 * Data		: 
 * SVN		:
 * Description	:
 * ***********************************************************/

#ifndef UCP_H
#define UCP_H

#include <string>

#include "cache.h"

typedef enum
{
	UMON_LOCAL	= 0,
	UMON_GLOBAL	= 1,
	UMON_DSS	= 2
} UMON_Type;

class ucp_cache_c: public cache_c
{
	public:
		ucp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
				bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
				int num_tiles, int interleave_factor, macsim_c* simBase,
				UMON_Type umonType);

		virtual ~ucp_cache_c();

		void update_cache_on_access(Addr tag, int set, int appl_id);
		//		void update_line_on_hit(cache_entry_c* line, int set, int appl_id);
		//		void update_cache_on_miss(int set_id, int appl_id);
		//		void update_set_on_replacement(Addr tag, int appl_id, int set_id, bool gpuline);
		cache_entry_c* find_replacement_line(int set, int appl_id);
		void ATD_initialize_cache_line(cache_entry_c *ins_line, Addr tag, 
				Addr addr, int set_id, bool skip);

		double get_max_mu(int appl_id, int alloc, int balance, int *blk_reqs);
		double get_mu_value(int appl_id, int a, int b);
		void calculate_ways(void);	
		void update_cache_policy(Counter m_cycles);

	protected:
		int			        m_num_cores;		/**< Number of cores */
		cache_set_c***			m_ATD_set;	    	/**< Auxiliary tag directories */
		uns32*          	    	m_num_access;   	/**< Number of accesses per core */
		uns32**				m_num_misses;		/**< Number of misses per core */
		uns32**			    	m_num_hits;	    	/**< Number of hits per core */
		uns32**			    	m_utility;	    	/**< Utility */
		UMON_Type		    	m_umon_type;		/**< Type of UMON */
		int*			    	m_allocations;		/**< Allocated ways per core */
		int*			    	m_occupied;	    	/**< Occupied ways per core */
		int			        m_num_sampling;  	/**< Number of sampling sets */
};

#endif
