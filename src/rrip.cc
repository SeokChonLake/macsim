#include "assert_macros.h"
#include "rrip.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

#define LONG_RRPV       ((1 << m_rrpv_len) - 2)
#define DISTANT_RRPV    ((1 << m_rrpv_len) - 1)

#define GET_PROB (static_cast <double> (rand()) / static_cast <float> (RAND_MAX))

rrip_cache_c::rrip_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
		bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
		int num_tiles, int interleave_factor, macsim_c* simBase,
		SDM_Type sdmType)
: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
		cache_by_pass, core_id, cache_type_info, enable_partition,
		num_tiles, interleave_factor, simBase)
{
	m_num_cores=*KNOB(KNOB_NUM_SIM_CORES);

	/* random seed */
	srand (static_cast <unsigned> (time(0)));

	// Sampling Type
	m_SDM_type = sdmType;
	m_pstate = new Policy_Type[num_set];
	memset(m_pstate, 0, sizeof(Policy_Type) * num_set);
	switch(m_SDM_type)
	{
		case SDM_LOCAL:
			break;
		case SDM_GLOBAL:
			m_num_leader_sets = num_set;
			for (int ii = 0; ii < num_set; ii++)
				if (ii % 2)
					m_pstate[ii] = BRRIP;
				else
					break;
		case SDM_DSS:
			m_num_leader_sets = 8;
			for (int ii = 0; ii < m_num_leader_sets; ++ii)
			{
				m_pstate[ii*2] = SRRIP;
				m_pstate[ii*2 + 1] = BRRIP;
			}
			for (int ii = 0; ii < num_set; ++ii)
				if(m_pstate[ii] == NONE)
					m_pstate[ii] = FOLLOWER;
			break;
	}

	memset(m_psel, 0, sizeof(uns32) * 4);
	m_rrpv_len = 5;
	m_epsilon = 1/32; 

}

rrip_cache_c::~rrip_cache_c()// : ~cache_c()
{
	delete[] m_pstate;
}

void rrip_cache_c::update_cache_on_miss(Addr addr, int set, int appl_id)
{
	switch(m_pstate[set]) {
		case NONE:
			break;
		case FOLLOWER:
			break;
		case SRRIP:
			m_psel[SRRIP]++;
			break;
		case BRRIP:
			m_psel[BRRIP]++;
			break;
		default:
			exit(-1);
	}
}

void rrip_cache_c::initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, int appl_id, bool gpuline, int set_id, bool skip)
{
	ins_line->m_valid           = true;
	ins_line->m_tag             = tag;
	ins_line->m_base            = (addr & ~m_offset_mask);
	ins_line->m_access_counter  = 0;
	ins_line->m_pref            = false;
	ins_line->m_skip            = skip;

	ins_line->m_appl_id         = appl_id;
	ins_line->m_gpuline         = gpuline;

	if (ins_line->m_gpuline) {
		++m_num_gpu_line;
		++m_set[set_id]->m_num_gpu_line;
	}
	else {
		++m_num_cpu_line;
		++m_set[set_id]->m_num_cpu_line;
	}
	// Select insertion policy
	switch(m_pstate[set_id])
	{
		case NONE:
		case FOLLOWER:
			if (m_psel[SRRIP] > m_psel[BRRIP])
				goto LABEL_BRRIP;
		case SRRIP:
			ins_line->m_last_access_time = LONG_RRPV;
			break;
		case BRRIP:
LABEL_BRRIP:
			if (GET_PROB <= m_epsilon)
				ins_line->m_last_access_time = LONG_RRPV;
			else
				ins_line->m_last_access_time = DISTANT_RRPV;
			break;
	}
}

// Extract victim cache line from Long RRPV position
cache_entry_c* rrip_cache_c::find_replacement_line(int set, int appl_id)
{
	int RRIP_tail_ind = 0;
	// If free entry found, return it
	for (int ii = 0; ii < m_assoc; ii++)
	{
		cache_entry_c* line = &(m_set[set]->m_entry[ii]);
		if (!line->m_valid) {
			return line;
		}
	}

	while(1) {
		RRIP_tail_ind = -1;
		// if the long re-reference interval prediction value entry found, return it
		for(int ii = 0; ii < m_assoc; ii++)
		{
			cache_entry_c* line = &(m_set[set]->m_entry[ii]);
			// use m_last_access_time as a rrpv
			if(line->m_last_access_time >= DISTANT_RRPV)
			{
				// TODO: Need to break a tie 
				RRIP_tail_ind = ii;
				return &(m_set[set]->m_entry[ii]);
			}
		}
		if (RRIP_tail_ind == -1) // Can't find long RRPV, do aging for all entry
		{
			for (int ii = 0; ii < m_assoc; ii++)
			{
				cache_entry_c* line = &(m_set[set]->m_entry[ii]);
				line->m_last_access_time++;
			}
		}
		else{
		    	// long RRPV entry found, return it
			return &(m_set[set]->m_entry[RRIP_tail_ind]);
		}
	}
}

void rrip_cache_c::update_line_on_hit(cache_entry_c* line, int set, int appl_id)
{
	line->m_last_access_time = 0;
}

void rrip_cache_c::update_cache_policy(Counter m_cycle)
{
	static uns32 check = 1; 
	uns32 check_point = m_cycle / 500000;

	if (check_point > check)
	{
		m_psel[BRRIP] = 0;
		m_psel[SRRIP] = 0;

		check++;
	}

}
