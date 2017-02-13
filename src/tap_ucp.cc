#include "assert_macros.h"
#include "tap_ucp.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"
#include "retire.h"
#include <algorithm>

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

tap_ucp_cache_c::tap_ucp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
			 bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
			 int num_tiles, int interleave_factor, macsim_c* simBase,
			 UMON_Type umonType)
	: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
		  cache_by_pass, core_id, cache_type_info, enable_partition,
		  num_tiles, interleave_factor, simBase)
{
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
		m_allocations[ii] = m_num_sets / m_num_cores;
		m_occupied[ii] = 0;
        m_num_access[ii] = 0;

		for (int jj = 0; jj < m_assoc; ++jj)
		{
			m_num_misses[ii][jj] = 0;
			m_num_hits[ii][jj] = 0;
			m_utility[ii][jj] = 0;
		}
	}

    // For TAP
    m_core_type = new CORE_Type[m_num_cores];
    m_core_policy_type = new Core_Policy_Type[m_num_cores];
    int leader_count = 2;
    for (int ii = 0; ii < m_num_cores; ii++)
    {
        m_core_policy_type[ii] = LRU;
        string type_name = m_simBase->m_core_pointers[ii]->get_core_type();
        if (type_name.compare("x86") == 0)
            m_core_type[ii] = CPU;
        else if (type_name.compare("ptx") == 0)
        {
            if (leader_count)
            {
                m_core_type[ii] = GPGPU_LEADER;
                if (leader_count == 2) 
                    m_core_policy_type[ii] = LIP;
                m_leaders[--leader_count] = ii;
            }
            else
            {
                m_core_type[ii] = GPGPU_FOLLOWER;
                m_core_policy_type[ii] = LIP; // maybe?
            }
        }
        //TODO: ACCELERATOR
        
        UCP_Mask = 0;
        memset(m_inst_count, 0, 2 * sizeof(Counter));
        memset(m_cycle_count, 0, 2 * sizeof(Counter));
    }
}

tap_ucp_cache_c::~tap_ucp_cache_c()// : ~cache_c()
{
    for (int ii = 0; ii < m_num_cores; ++ii)
    {
        for (int jj = 0; jj < m_num_sets; ++jj)
            delete[] m_ATD_set[ii][jj];
        delete[] m_ATD_set[ii];
    }

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

    delete[] m_core_type;
    delete[] m_core_policy_type;
}

static bool comp_func(uns32 i, uns32 j)
{
    return (i > j);
}

void tap_ucp_cache_c::update_cache_on_access(Addr addr, int set, int appl_id)
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
cache_entry_c* tap_ucp_cache_c::find_replacement_line(int set, int appl_id)
{
	int lru_ind = -1;
	Counter lru_time = MAX_INT;

	// Initialize the number of blocks belong to each core
	for (int ii = 0; ii < m_assoc; ++ii)
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
				if (appl_id == line->m_appl_id) //victim cannot be from miss-causing apps, 
					continue;
				else if (m_occupied[line->m_appl_id] < m_allocations[line->m_appl_id]) // victim cannot be from other less occupied apps's entry
					continue;
			}
			else	// find victim entry from miss_causing applications 
			{
				if (appl_id != line->m_appl_id) // find victim only from miss-causing apps
				       continue;	
			}
			lru_ind  = ii;
			lru_time = line->m_last_access_time;
		}
	}

	assert(lru_ind != -1);

	return &(m_set[set]->m_entry[lru_ind]);
}

void tap_ucp_cache_c::ATD_initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, int set_id, bool skip)
{
	ins_line->m_valid            = true;
	ins_line->m_tag              = tag;
	ins_line->m_base             = (addr & ~m_offset_mask);
	ins_line->m_access_counter   = 0;
	ins_line->m_last_access_time = CYCLE;
	ins_line->m_pref             = false;
	ins_line->m_skip             = skip;
}

void tap_ucp_cache_c::initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr, int appl_id, bool gpuline, int set_id, bool skip)
{
	ins_line->m_valid            = true;
	ins_line->m_tag              = tag;
	ins_line->m_base             = (addr & ~m_offset_mask);
	ins_line->m_access_counter   = 0;
	ins_line->m_last_access_time = CYCLE;
	ins_line->m_pref             = false;
	ins_line->m_skip             = skip;

    ins_line->m_appl_id          = appl_id;
    ins_line->m_gpuline          = gpuline;

    if(gpuline) {
        ++m_num_gpu_line;
        ++m_set[set_id]->m_num_gpu_line;
    }
    else
    {
        ++m_num_cpu_line;
        ++m_set[set_id]->m_num_cpu_line;
    }

    // Core sampling:: first gpu core uses policy 1(LRU), second gpu core uses policy 2 (LIP)
    if(m_core_policy_type[appl_id] == LIP)
        ins_line->m_last_access_time = 1;
    
}
double tap_ucp_cache_c::get_max_mu(int appl_id, int alloc, int balance, int *blk_reqs)
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

double tap_ucp_cache_c::get_mu_value(int appl_id, int a, int b)
{
    double mu = m_num_misses[appl_id][a] - m_num_misses[appl_id][b];
    return mu/(b - a);
}

void tap_ucp_cache_c::calculate_ways(void)
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

    uns32 sum_access = 0;
    uns32 max_cpu_access = 0;
    int XSRATIO = 1;
    for (int ii = 0; ii < m_num_cores; ++ii)
    {
        if (m_core_type[ii] == GPGPU_LEADER || m_core_type[ii] == GPGPU_FOLLOWER)
            sum_access += m_num_access[ii];
        else {
            if (max_cpu_access < m_num_access[ii])
                max_cpu_access = m_num_access[ii];
        }
    }

    if (max_cpu_access * 10 < sum_access)
        XSRATIO = sum_access / max_cpu_access;


    memset(m_allocations, 0, m_num_cores * sizeof(int));

    //TAP-UCP begin
    if (XSRATIO > 1)
        for (int ii = 0; ii < m_num_cores; ++ii)
        {
            if (m_core_type[ii] == GPGPU_LEADER || m_core_type[ii] == GPGPU_FOLLOWER)
                m_num_access[ii] /= XSRATIO;
        }
    //TAP-UCP end

    while (balance) {
        
        for (int ii = 0; ii < m_num_cores; ++ii) {
            //TAP-UCP begin
            if(m_core_type[ii] == GPGPU_LEADER || m_core_type[ii] == GPGPU_FOLLOWER)
                if(UCP_Mask)
                    continue;
            //TAP-UCP end
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

void tap_ucp_cache_c::update_cache_policy(Counter m_cycles)
{
    static uns32 check = 0;
    
    double leaders_cpi[2];

    uns32 check_point = m_cycles / 5000000;
    if (check_point >= check)
    {
        // calculate current interval's cpi
        for (int i = 0; i < 2; i++) {
            m_inst_count[i] = (m_simBase->m_core_pointers[m_leaders[i]]->get_retire()->get_total_insts_retired()
                    - m_inst_count[i]);
            m_cycle_count[i] = (m_simBase->m_core_pointers[m_leaders[i]]->get_cycle_count() 
                    - m_cycle_count[i]);
            leaders_cpi[i] = (double)(m_cycle_count[i] / m_inst_count[i]);
        }

        double cpi_delta = leaders_cpi[0] - leaders_cpi[1];
        if (cpi_delta < 0) 
            cpi_delta *= -1; 

        if(cpi_delta * 20 < leaders_cpi[0] || cpi_delta * 20 < leaders_cpi[1])
            UCP_Mask = 1; // caching gpgpu application is not effective
        else
            UCP_Mask = 0;

        if (UCP_Mask) // caching is not effective, following cache-against policy
        {
            for (int ii = 0; ii < m_num_cores; ii++)
                if(m_core_type[ii] == GPGPU_FOLLOWER)
                    m_core_policy_type[ii] = LIP;
        }
        else    // caching is effective, following cache-frinedly policy
        {
            for (int ii = 0; ii < m_num_cores; ii++)
                if(m_core_type[ii] == GPGPU_FOLLOWER)
                    m_core_policy_type[ii] = LRU;
        }

        cout << "Caclulate ways @ " << m_cycles << endl;
        calculate_ways();
        check++;
    }

}
