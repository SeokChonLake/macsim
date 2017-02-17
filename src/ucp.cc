#include "assert_macros.h"
#include "ucp.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"

#include <algorithm>

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

ucp_cache_c::ucp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
		bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
		int num_tiles, int interleave_factor, macsim_c* simBase,
		UMON_Type umonType)
: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
		cache_by_pass, core_id, cache_type_info, enable_partition,
		num_tiles, interleave_factor, simBase)
{
	m_num_cores = *KNOB(KNOB_NUM_SIM_CORES);

	// Allocating memory for Auxiliary Tag Directory for every other cores
	m_ATD_set = new cache_set_c**[m_num_cores];
	for (int ii = 0; ii < m_num_cores; ++ii) {
		m_ATD_set[ii] = new cache_set_c*[m_num_sets];
		for (int jj = 0; jj < m_num_sets; ++jj)
		{
			m_ATD_set[ii][jj] = new cache_set_c(m_assoc);
			for (int kk = 0; kk < assoc; ++kk) 
			{
				m_ATD_set[ii][jj]->m_entry[kk].m_valid		= false;
				m_ATD_set[ii][jj]->m_entry[kk].m_access_counter = false;
				m_ATD_set[ii][jj]->m_entry[kk].m_data		= INIT_CACHE_DATA_VALUE;
			}
		}
	}

	// Sampling Type
	m_umon_type = umonType;
	switch(m_umon_type)
	{
		case UMON_LOCAL:
			m_num_sampling = num_set;
			break;
		case UMON_GLOBAL:
			m_num_sampling = 32;
			break;
		case UMON_DSS:
			break;
	}

	// Allocating other informations
	m_allocations = new int[m_num_cores];
	m_occupied = new int[m_num_cores];

	m_num_access = new uns32[m_num_cores];
	m_num_misses = new uns32*[m_num_cores];
	m_num_hits = new uns32*[m_num_cores];
	m_utility = new uns32*[m_num_cores];
	for(int ii = 0; ii < m_num_cores; ii++) {
		m_num_misses[ii] = new uns32[m_assoc];
		m_num_hits[ii] = new uns32[m_assoc];
		m_utility[ii] = new uns32[m_assoc];
	}

	// Initializing 
	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		m_allocations[ii] = m_assoc / m_num_cores;
		cout << m_allocations[ii] << " ";
		m_occupied[ii] = 0;
		m_num_access[ii] = 0;

		for (int jj = 0; jj < m_assoc; ++jj)
		{
			m_num_misses[ii][jj] = 0;
			m_num_hits[ii][jj] = 0;
			m_utility[ii][jj] = 0;
		}
	}
	cout << endl;
}

ucp_cache_c::~ucp_cache_c()// : ~cache_c()
{
	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		for (int jj = 0; jj < m_num_sets; ++jj)
			delete m_ATD_set[ii][jj];
		delete[] m_ATD_set[ii];
	}
	delete[] m_ATD_set;

	delete[] m_allocations;
	delete[] m_occupied;

	for (int ii = 0; ii < m_num_cores; ++ii) 
	{
		delete[] m_num_misses[ii];
		delete[] m_num_hits[ii];
		delete[] m_utility[ii];
	}
	delete[] m_num_misses;
	delete[] m_num_hits;
	delete[] m_utility;
}

static bool comp_func(uns32 i, uns32 j)
{
	return (i > j);
}

void ucp_cache_c::update_cache_on_access(Addr addr, int set, int appl_id)
{
	cache_set_c** ATD_set = m_ATD_set[appl_id];

	if (m_cache_by_pass)
		return ;

	Addr tag;
	int i = 0;
	int lru_ind = 0;
	Counter lru_time = MAX_INT;

	find_tag_and_set(addr, &tag, &set);

	// Increment the access counter
	m_num_access[appl_id]++;

	// Walk through the ATD set
	for (int ii = 0; ii < m_assoc; ++ii) {
		// For each line in based on associativity
		cache_entry_c * line = &(ATD_set[set]->m_entry[ii]);

		// Check for matching tag and validity
		if (line->m_valid && line->m_tag == tag) {
			// If hit, increment the hit counter
			m_num_hits[appl_id][ii]++;
			return ;
		}
	}

	// If miss incurred, find victim and replace the entry based on the policy
	while (i < m_assoc) {
		cache_entry_c* line = &(ATD_set[set]->m_entry[i]);
		// If free entry found, return it
		if (!line->m_valid) {
			lru_ind = i;
			break;
		}

		// Check if this is the LRU entry encountered
		if (line->m_last_access_time < lru_time) {
			lru_ind    = i;
			lru_time   = line->m_last_access_time;
		}
		++i;
	}

	// Initialize the cache line
	ATD_initialize_cache_line(&(ATD_set[set]->m_entry[lru_ind]), tag, addr, set, false);
}


/*
   void update_line_on_hit(cache_entry_c* line, int set, int appl_id)
   {
   }

   void update_cache_on_miss(int set_id, int appl_id)
   {
   }
   void update_set_on_replacement(Addr tag, int appl_id, int set_id, bool gpuline)
   {
   }
   */


// Modification to replacement engine for way-partitioing
// On a cache miss, the replacement engine counts the number of cache blocks
// that belong to the miss-causing applications in the set. If this number is less
// than the number of blocks allocated to the application, then the LRU block
// among all the blocks that do not belong to the application is evicted.
// Otherwise, the LRU block among all the blocks of the miss-causing application is evicted.
cache_entry_c* ucp_cache_c::find_replacement_line(int set, int appl_id)
{
	int lru_ind = -1;
	Counter lru_time = MAX_INT;

	// Initialize the number of blocks belong to each core
	for (int ii = 0; ii < m_num_cores; ++ii)
		m_occupied[ii] = 0;

	// Count the number of occupied blocks per core
	for (int ii = 0; ii < m_assoc; ++ii)
		m_occupied[m_set[set]->m_entry[ii].m_appl_id]++;

	/* Find the invalid cache line or LRU cache line */
	for (int ii = 0; ii < m_assoc; ++ii)
	{
		cache_entry_c* line = &(m_set[set]->m_entry[ii]);
		// If free entry found, return it
		if (!line->m_valid)
		{
			lru_ind = ii;
			break;
		}

		// Check if this is the LRU entry encountered
		if (line->m_last_access_time < lru_time) {
			// find victim entry from other over occupied applications
			if (m_occupied[appl_id] < m_allocations[appl_id])
			{
				if (appl_id == line->m_appl_id) 
					//victim cannot be from miss-causing apps, 
					continue;
				else if (m_occupied[line->m_appl_id] < m_allocations[line->m_appl_id]) 
					// victim cannot be from other less occupied apps's entry
					continue;
			}
			else	// find victim entry from miss_causing applications 
			{
				if (appl_id != line->m_appl_id) 
					// find victim only from miss-causing apps
					continue;	
			}
			lru_ind  = ii;
			lru_time = line->m_last_access_time;
		}
	}
	
	if (lru_ind == -1)
	{
		assert(lru_ind != -1);
	}

	return &(m_set[set]->m_entry[lru_ind]);
}

void ucp_cache_c::ATD_initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, int set_id, bool skip)
{
	ins_line->m_valid            = true;
	ins_line->m_tag              = tag;
	ins_line->m_base             = (addr & ~m_offset_mask);
	ins_line->m_access_counter   = 0;
	ins_line->m_last_access_time = CYCLE;
	ins_line->m_pref             = false;
	ins_line->m_skip             = skip;
}

double ucp_cache_c::get_max_mu(int appl_id, int alloc, int balance, int *blk_reqs)
{
	double max_mu = 0;

	for(int ii = 1; ii <= balance; ii++) {
		int mu = get_mu_value(appl_id, alloc, alloc+ii);
		if(mu > max_mu)
		{
			max_mu = mu; 
			*blk_reqs = ii;
		}
	}
	return max_mu;
}

double ucp_cache_c::get_mu_value(int appl_id, int a, int b)
{
	double mu = m_num_misses[appl_id][a] - m_num_misses[appl_id][b];
	return mu/(b - a);
}

void ucp_cache_c::calculate_ways(void)
{
	int balance = m_assoc;
	int* max_mu = new int[m_num_cores];
	int* blk_reqs = new int[m_num_cores];

	/* Set up the access, missess, hits tables */
	for (int ii = 0; ii < m_num_cores; ++ii)
		std::sort(m_num_hits[ii], m_num_hits[ii] + m_assoc, std::greater<uns32>());
	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		int hit_sum = 0;
		for (int jj = 0; jj < m_assoc; ++jj)
		{
			hit_sum += m_num_hits[ii][jj];
			m_num_misses[ii][jj] = m_num_access[ii] - hit_sum; 
		}
	}

	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		cout << "cores " << ii << endl;
		cout << ".access: " << m_num_access[ii] << endl;
		cout << ".hits: ";
		for (int jj = 0; jj < m_assoc; ++jj)
			cout << " " << m_num_hits[ii][jj];
		cout << endl;
		uns32 temp = 0;
		cout << ".misses: ";
		for (int jj = 0; jj < m_assoc; ++jj) {
			temp += m_num_hits[ii][jj];
			cout << " " << m_num_access[ii] - temp;
		}
		cout << endl;
	}

	memset(m_allocations, 0, m_num_cores * sizeof(int));
	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		m_allocations[ii] = 1;
		balance--;
	}

	while (balance) {

		for (int ii = 0; ii < m_num_cores; ++ii) {
			int alloc = m_allocations[ii];
			max_mu[ii] = get_max_mu(ii, alloc, balance, &blk_reqs[ii]);
		}
		int winner = 0;
		for (int ii = 1; ii < m_num_cores; ++ii)
			if (max_mu[ii] > max_mu[winner])
				winner = ii;
		m_allocations[winner] += blk_reqs[winner];
		balance -= blk_reqs[winner];
	}

	for (int ii = 0; ii < m_num_cores; ++ii)
	{
		memset(m_num_hits[ii], 0, m_assoc * sizeof(uns32));
		memset(m_num_misses[ii], 0, m_assoc * sizeof(uns32));
	}
	memset(m_num_access, 0, m_num_cores * sizeof(uns32));

	delete[] max_mu;
	delete[] blk_reqs;
}

void ucp_cache_c::update_cache_policy(Counter m_cycles)
{
	static uns32 check = 1;

	uns32 check_point = m_cycles / 500000;
	if (check_point >= check)
	{
		cout << "Caclulate ways @ " << m_cycles << endl;
		
		for(int ii = 0; ii < m_num_cores; ++ii)
			cout << m_allocations[ii] << ",";
		cout << endl;
		
		calculate_ways();
	
		for(int ii = 0; ii < m_num_cores; ++ii)
			cout << m_allocations[ii] << ",";
		cout << endl;

		check++;
	}

}
