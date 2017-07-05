/************************************************************
 * File		: rrip.h
 * Author	:
 * Data		: 
 * SVN		:
 * Description	:
 * ***********************************************************/

#ifndef RRIP_H
#define RRIP_H

#include <string>

#include "cache.h"

typedef enum
{
	SDM_LOCAL	= 0,
	SDM_GLOBAL	= 1,
	SDM_DSS	= 2
} SDM_Type;

typedef enum
{
    NONE        = 0,
    FOLLOWER    = 1,
    SRRIP       = 2,
    BRRIP       = 3
} Policy_Type;

class rrip_cache_c: public cache_c
{
	public:
		rrip_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
				bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
				int num_tiles, int interleave_factor, macsim_c* simBase,
				SDM_Type sdmType);

		virtual ~rrip_cache_c();

		void update_cache_on_miss(Addr addr, int set,  int appl_id);
		void update_line_on_hit(cache_entry_c* line, int set, int appl_id);
//		void update_cache_on_access(Addr tag, int set, int appl_id);
		cache_entry_c* find_replacement_line(int set, int appl_id);
		void update_cache_policy(Counter m_cycle);
		void initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, int appl_id, bool gpuline, int set_id, bool skip);

		Policy_Type determine_policy(void);
		void psel_change(int inc_or_dec);

	protected:
		int                 m_rrpv_len;           /**< Length of rrpv */
		uns32               m_psel[4];         /**< Policy selector */
		int			        m_num_cores;	/**< Number of cores */
		Policy_Type*        m_pstate;       /**< Policy state per set */
		SDM_Type            m_SDM_type;     /**< Type of Set Dueling Monitor */
		int                 m_num_leader_sets; /**< Number of learder sets */
		double              m_epsilon;      /**< Epsilon value for BRRIP */

};

#endif
