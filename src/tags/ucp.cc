#include "assert_macros.h"
#include "tags/ucp.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

ucp_cache_c::ucp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
			 bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
			 int num_tiles, int interleave_factor, macsim_c* simBase)
	: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
		  cache_by_pass, core_id, cache_type_info, enable_partition,
		  num_tiles, interleave_factor, simBase)
{
}

ucp_cache_c::~ucp_cache_c()
{
}

void update_cache_on_access(Addr tag, int set, int appl_id);
void update_line_on_hit(cache_entry_c* line, int set, int appl_id);
void update_cache_on_miss(int set_id, int appl_id);
void update_set_on_replacement(Addr tag, int appl_id, int set_id, bool gpuline);
cache_entry_c* find_replacement_line(int set, int appl_id);


double get_max_mu(int core_id, int alloc) /* Look-ahead algorithm */
{
	int max_mu = 0;
	
	for (int ii = 1; ii <= m_assoc; ii++)
	{
		int mu = get_mu_value(core_id, alloc, alloc+ii);
		if(mu > max_mu) max_mu = mu;
	}
	
	return max_mu;
}

double get_mu_value(int core_id, int a, int b)
{
	return (double)((m_num_misses[core_id][a] - m_num_misses[core_id][b]) / (b - a));
}

void set_partition()
{
	int balance = m_assoc;
	int* max_mu = new int[m_num_cores];
	int* block_reqs = new int[m_num_cores];

	for (int ii = 0; ii < m_assoc; ii++)
		m_allocations[ii] = 0;

	while(balance) {
		for (int ii = 0; ii < m_num_cores; ii++) /* get max marginal utiltiy */
		{
			int alloc = m_allocations[ii];
			max_mu[ii] = get_max_mu(ii, alloc, balance);
			block_reqs[ii]

	}
}
