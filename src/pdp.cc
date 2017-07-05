#include "assert_macros.h"
#include "pdp.h"
#include "utils.h"

#include "debug_macros.h"

#include "all_knobs.h"

#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define DEBUG(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_CACHE_LIB, ## args)
#define DEBUG_MEM(args...) _DEBUG(*m_simBase->m_knobs->KNOB_DEBUG_MEM_TRACE, ## args)

#define ENABLE_CAPTURE_RDD 1


pdp_cache_c::pdp_cache_c(string name, int num_set, int assoc, int line_size, int data_size, int bank_num,
		bool cache_by_pass, int core_id, Cache_Type cache_type_info, bool enable_partition,
		int num_tiles, int interleave_factor, macsim_c* simBase
		)
: cache_c(name, num_set, assoc, line_size, data_size, bank_num,
		cache_by_pass, core_id, cache_type_info, enable_partition,
		num_tiles, interleave_factor, simBase)
{
    m_num_cores = *KNOB(KNOB_NUM_SIM_CORES);
    m_max_rd = 1024;
	m_num_set = num_set;
    m_pd = new int[m_num_cores];
    m_access_counter = new Counter[num_set];
    m_rd_counter = new Counter*[m_num_cores];
    m_ip_table = new bool[m_num_cores];
    m_allocation = new int[m_num_cores];
    m_occupied = new int[m_num_cores];

	pre_hit_rate = new int[m_num_cores];

	real_rd_flag= false;
    /* Initialization */
    for (int ii = 0; ii < m_num_cores; ii++) {
        // initially set the PD as occupied ways for each master IP

        m_rd_counter[ii] = new Counter[m_max_rd];
		

        for (int jj = 0; jj < m_max_rd; jj++) {
            m_rd_counter[ii][jj] = 0;
        }
        m_allocation[ii] = assoc / m_num_cores;

        /* FIX ME */
        cin >> m_ip_table[ii];
		if(m_ip_table[ii]==0)//cpu
			m_pd[ii]=50;
		else
			m_pd[ii]=m_max_rd;

		cout << m_ip_table[ii] <<endl;
    }
	/* for test */

	for (int ii=0;ii<m_num_cores;ii++){
		cin >> m_allocation[ii];
		cout <<"Way" << ii <<":" << m_allocation[ii]<<endl;
	}
	
		
	
	cout << "ip_table done"<<endl;
	for (int ii = 0; ii<m_num_cores; ii++)
	    cout << "way"<<ii<<"="<<m_allocation[ii]<<endl;
	for (int ii = 0; ii < num_set; ii++) {
        m_access_counter[ii] = 0;
    }
}

pdp_cache_c::~pdp_cache_c()
{
#if ENABLE_CAPTURE_RDD
    ofstream fp;
    fp.open("eviction.csv", ios_base::out);
	fp << "access";
	for (int ii = 0; ii < m_num_cores; ii++){
		fp<<",reuse distance for all"<<ii;
	}
	fp << endl;



    for (int ii = 0; ii < m_max_rd; ii++) {
        fp << ii+1;
        for (int jj = 0; jj < m_num_cores; jj++) {
            fp << "," << m_rd_counter[jj][ii];
        }
        fp << endl;
    }
    fp.close();
#endif

    delete[] m_pd;
    delete[] m_access_counter;
    delete[] m_ip_table;
    delete[] m_allocation;
    delete[] m_occupied;

    for (int ii = 0; ii < m_num_cores; ii++)
        delete[] m_rd_counter[ii];
    delete[] m_rd_counter;
}

/* walk through the sampler and find history
   if found, update the RD counter and re-marking access stamp
   if not found, push the new line to sampler */
void pdp_cache_c::update_reuse_distance(Addr tag, int set, int appl_id) {
    
    m_access_counter[set]++;

    // find inserted line in sampler
    bool found = false;
    int index = -1;
    rd_record temp;
    for (int ii = 0; ii < m_rd_sampler.size(); ii++) {
        temp = m_rd_sampler[ii];
        if (temp.tag == tag && temp.set == set && temp.appl_id == appl_id) {
            found = true; index = ii; break;
        }
    }

    // if not found
    if (!found) {
        rd_record temp = rd_record(tag, set, m_access_counter[set], appl_id);
        m_rd_sampler.push_back(temp);
    }
    // if found
    else {
        Counter rd = m_access_counter[set] - m_rd_sampler[index].count - 1;    // calculate reuse distance


        if (rd >= m_max_rd)
            m_rd_counter[appl_id][m_max_rd - 1]++;
        else
            m_rd_counter[appl_id][rd]++;
        m_rd_sampler[index].count = m_access_counter[set];
        
    }
}

/* Call update reuse distance
 */
void pdp_cache_c::update_cache_on_access(Addr tag, int set, int appl_id) {
    update_reuse_distance(tag, set, appl_id);
}

void pdp_cache_c::initialize_cache_line(cache_entry_c *ins_line, Addr tag, Addr addr,
        int appl_id, bool gpuline, int set_id, bool skip)
{
    ins_line->m_valid               = true;
    ins_line->m_tag                 = tag;
    ins_line->m_base                = (addr & ~m_offset_mask);
    ins_line->m_access_counter      = 0;
    ins_line->m_last_access_time    = CYCLE;
    ins_line->m_pref                = false;
    ins_line->m_skip                = skip;
    ins_line->m_appl_id             = appl_id;
    ins_line->m_gpuline             = gpuline;

    // if the inserted line is from IP_CORE, set the protected bit
    if (m_ip_table[appl_id])
        ins_line->m_protected = true;
    // if from CPU_CORE, clear the protecteed bit
    else
        ins_line->m_protected = false;

//    ins_line->m_pd = m_pd[appl_id];
    ins_line->m_rd = m_pd[appl_id];
	//for check real reuse distance
	ins_line->m_last_set_access_time = m_access_counter[set_id];
	ins_line->m_hit_count = 0;
}

/*
   to get locking effect, find the available room for inserted line
   before victim selection phase.
 */
cache_entry_c* pdp_cache_c::is_there_a_room(int set, int appl_id) {
    int margin[m_num_cores]; // allocations - occupied

    for (int ii = 0; ii < m_num_cores; ii++)
        margin[ii] = m_allocation[ii];


    for (int ii = 0; ii < m_assoc; ii++) {
        cache_entry_c *line = &(m_set[set]->m_entry[ii]);
        if (line->m_valid) {
            margin[line->m_appl_id]--;
        }
	/*	else 
			return line;*/
			
    }

    if (margin[appl_id] <= 0) {
        // if inserted line is from over-occupied master IP,
        // victim only can select from own expired entry
        int victim_index = -1;
        int max_rrd = -1;
		//cout << " good " <<endl;
        for (int ii = 0; ii < m_assoc; ii++) {
            cache_entry_c *line = &(m_set[set]->m_entry[ii]);
            if (line->m_valid && line->m_appl_id == appl_id && line->m_protected == false)
            {
                // if encounter expired line with maximum rrd, can be victim entry
                if (max_rrd < line->m_rd)
                {
                    max_rrd = line->m_rd;
                    victim_index = ii;
                }
            }
        }

        if (victim_index != -1){
         	cache_entry_c *line = &(m_set[set]->m_entry[victim_index]);
		    return &(m_set[set]->m_entry[victim_index]);
		}
    }
    else {
        // if inserted line is from less-occupied master IP,
        // victim can be select from own and other over-occupied master IP
        int victim_index = -1;
        int max_rrd = -1;
        for (int ii = 0; ii < m_assoc; ii++) {
            cache_entry_c *line = &(m_set[set]->m_entry[ii]);
          /*  if (line->m_valid && line->m_protected == false &&
                    (line->m_appl_id == appl_id || margin[line->m_appl_id] < 0))*/// original
  			if (!line->m_valid) return line;         // less occupied can select invaild line
            if (line->m_valid && line->m_protected == false &&
                    ( margin[line->m_appl_id] < 0||line->m_appl_id == appl_id))
			{
                // if encounter expired line with maximum rrd, can be victim entry
                if (max_rrd < line->m_rd)
                {
                    max_rrd = line->m_rd;
                    victim_index = ii;
                }
            }
        }
//		cout  << " no " << endl;
        if (victim_index != -1){

         	cache_entry_c *line = &(m_set[set]->m_entry[victim_index]);
		    return &(m_set[set]->m_entry[victim_index]);

		}
    }
//	cout << "all protected" << endl;
    return NULL;

}

/*
   victim selection phase of PDP
   
 */
cache_entry_c* pdp_cache_c::find_replacement_line(int set, int appl_id)
{
    // 1. only if there is an available room for inserted line,
    //    can enter this pahse(victim selection)
    assert(is_there_a_room(set, appl_id));

    return is_there_a_room(set, appl_id); 
}

void pdp_cache_c::update_cache_policy(Counter m_cycle){
	static uns32 check = 1;
	uns32 check_point = m_cycle /5000000; // 30m cycle	
	if(m_cycle%1000000==0) cout << m_cycle << endl;
	if (check_point >= check){
		//calculate PD - set PD to average of RDD	 
		// TODO: modified after
		cout << m_cycle;
		cout << ": 30m cycle"<<endl;
		calculate_pd();
		check++;
	}
}
void pdp_cache_c::calculate_way(void){
	int sum_way;
	int IP_count=0;
	
	sum_way = m_assoc;

	for(int ii = 0; ii <m_num_cores ; ii++){
		if(m_ip_table[ii]==0){
			for(int jj = 0; jj < m_assoc ;jj++){
				if(pre_hit_rate[ii]*jj>=10000){
					while(sum_way-jj<1)
						jj--;
					if(jj<-1) jj=-1; //for error
					m_allocation[ii] = jj + 1;	
					sum_way -= (jj + 1);
				}
			}
		}
	}

	for(int ii = 0; ii<m_num_cores ;ii++){
		if(m_ip_table[ii] == 1){
			m_pd[ii]=1024;
			for(int jj = 0; jj < m_assoc ;jj++){
				if(pre_hit_rate[ii]*jj>=10000){
					while(sum_way-jj<1)
						jj--;
					if(jj<-1) jj=-1;
					m_allocation[ii] = jj + 1;	
					sum_way -= (jj + 1);
				}
			}
		}
	}
	if(sum_way>0){
		for(int ii = 0;ii<m_num_cores;ii++){
			m_allocation[ii]++;
			sum_way--;
			if(sum_way==0) break;
		}
	}
	for(int ii=0; ii<m_num_cores;ii++){
		cout <<"way["<<ii<<"] =" <<m_allocation[ii]<<endl;
	}
}

void pdp_cache_c::calculate_pd(void){
	int real_rd;
	
	for(int i = 0 ; i < m_num_cores; i++){
		int N_t =0;
		double max_E_dp = 0;
		for( int j = 0 ; j < m_max_rd;j++)
			N_t += m_rd_counter[i][j];
		cout <<"total access:["<< N_t <<"]"<<endl;
		if(N_t != 0){
			for( int j = 0 ; j< 1024 ;j++){
				int N_i = 0;
				int N_i_mul = 0;
				double E_dp = 0;
				for( int k = 0; k < j ;k++){
					N_i+=m_rd_counter[i][k];
					N_i_mul+=m_rd_counter[i][k]*(k+1);
				}
				E_dp =(double) N_i;
				E_dp = E_dp/(double)(N_i_mul+((N_t - N_i) * (j+m_assoc)));

				if(max_E_dp < E_dp){
					max_E_dp = E_dp;
					m_pd[i] = j;
				}
			}
			cout << i << ":" << m_pd[i] <<":"<< max_E_dp << endl;
		}
		else
			cout << i << "(past):"<< m_pd[i] << endl;

		if(m_ip_table[i]==1) m_pd[i]=1024;
//		pre_hit_rate[i]=max_E_dp*10000;
	}
	
//	calculate_way();
}


void pdp_cache_c::update_line_on_hit(cache_entry_c* line, int set, int appl_id){
	int real_rd=0;
	line->m_rd = m_pd[appl_id];
	cache_set_promotion(set);

	real_rd=m_access_counter[set]-line->m_last_set_access_time;
	if(real_rd >= m_max_rd) real_rd = m_max_rd;
	if(real_rd!=0)  line->m_hit_count = real_rd;
//	m_real_counter[appl_id][real_rd]++;
	line->m_last_set_access_time=m_access_counter[set];


	//TODO : ip -> unprotect cpu -> protect
	if (m_ip_table[appl_id] == IP_CORE) {
		line->m_protected = true;
	}
	else {
		line->m_protected = true;
	}
}

void pdp_cache_c::update_cache_on_miss(int set, int appl_id){
	cache_set_promotion(set); 
}


void pdp_cache_c::cache_set_promotion(int set){
    for (int ii = 0; ii < m_assoc; ii++) {
        cache_entry_c *line = &(m_set[set]->m_entry[ii]);
        if(line->m_valid) {
			line->m_rd--;
            if(line->m_rd <= 0){
                line->m_rd = m_max_rd*2;
                line->m_protected = false;
            }
        }
    }
}

bool pdp_cache_c::bypass_monitor(Addr addr, int m_appl_id) 
{
        Addr line_addr;
        Addr tag;
        int set;

        find_tag_and_set(addr, &tag, &set);
        line_addr = base_cache_line(addr);

        if (is_there_a_room(set, m_appl_id) == NULL){
         	//cout << "bypass:"<< m_appl_id<<endl;
		    return true;
		}
        else
            return false;
}

