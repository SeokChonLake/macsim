#ifndef HMC_TYPES_H
#define HMC_TYPES_H

#include <string>

typedef enum HMC_Type_enum
{
    HMC_NONE=0,
    HMC_CAS_equal_16B=1,
    HMC_CAS_zero_16B=2,
    HMC_CAS_greater_16B=3,
    HMC_CAS_less_16B=4,
    HMC_ADD_16B=5,
    HMC_ADD_8B=6,
    HMC_ADD_DUAL=7,
    HMC_SWAP=8,
    HMC_BIT_WR=9,
    HMC_AND=10,
    HMC_NAND=11,
    HMC_OR=12,
    HMC_XOR=13,
    HMC_FP_ADD=14,
    HMC_COMP_greater=15,
    HMC_COMP_less=16,
    HMC_COMP_equal=17,
    HMC_CANDIDATE=18,
    NUM_HMC_TYPES=19
} HMC_Type;


class hmc_type_c
{
    public:    
    static std::string HMC_Type2String(HMC_Type hmctype)
    {
        switch(hmctype)
        {
            case HMC_NONE: return std::string("HMC_NONE");
            case HMC_CAS_equal_16B: return std::string("HMC_CAS_equal_16B");
            case HMC_CAS_zero_16B:  return std::string("HMC_CAS_zero_16B");
            case HMC_CAS_greater_16B:return std::string("HMC_CAS_greater_16B");
            case HMC_CAS_less_16B:  return std::string("HMC_CAS_less_16B");
            case HMC_ADD_16B:       return std::string("HMC_ADD_16B");
            case HMC_ADD_8B:        return std::string("HMC_ADD_8B");
            case HMC_ADD_DUAL:      return std::string("HMC_ADD_DUAL");
            case HMC_SWAP:          return std::string("HMC_SWAP");
            case HMC_BIT_WR:        return std::string("HMC_BIT_WR");
            case HMC_AND:           return std::string("HMC_AND");
            case HMC_NAND:          return std::string("HMC_NAND");
            case HMC_OR:            return std::string("HMC_OR");
            case HMC_XOR:           return std::string("HMC_XOR");
            case HMC_FP_ADD:        return std::string("HMC_FP_ADD");
            case HMC_COMP_greater:  return std::string("HMC_COMP_greater");
            case HMC_COMP_less:     return std::string("HMC_COMP_less");
            case HMC_COMP_equal:    return std::string("HMC_COMP_equal");
            case HMC_CANDIDATE:     return std::string("HMC_CANDIDATE");
            case NUM_HMC_TYPES:
            default: return std::string("");
        }
    }
    static HMC_Type  HMC_String2Type(std::string name)
    {
        if (name=="HMC_CAS_equal_16B") return HMC_CAS_equal_16B;
        else if (name=="HMC_CAS_zero_16B") return HMC_CAS_zero_16B;
        else if (name=="HMC_CAS_greater_16B") return HMC_CAS_greater_16B;
        else if (name=="HMC_CAS_less_16B") return HMC_CAS_less_16B;
        else if (name=="HMC_ADD_16B")      return HMC_ADD_16B;
        else if (name=="HMC_ADD_8B")       return HMC_ADD_8B;
        else if (name=="HMC_ADD_DUAL")     return HMC_ADD_DUAL;
        else if (name=="HMC_SWAP")         return HMC_SWAP;
        else if (name=="HMC_BIT_WR")       return HMC_BIT_WR;
        else if (name=="HMC_AND")          return HMC_AND;
        else if (name=="HMC_NAND")         return HMC_NAND;
        else if (name=="HMC_OR")           return HMC_OR;
        else if (name=="HMC_XOR")          return HMC_XOR;
        else if (name=="HMC_FP_ADD")       return HMC_FP_ADD;
        else if (name=="HMC_COMP_greater") return HMC_COMP_greater;
        else if (name=="HMC_COMP_less")    return HMC_COMP_less;
        else if (name=="HMC_COMP_equal")   return HMC_COMP_equal;
        else if (name=="HMC_CANDIDATE")    return HMC_CANDIDATE;
        else return HMC_NONE;
    }

};

#define HMC_EVENT_COUNT(core_id, hmctype) \
switch(hmctype) {\
case HMC_NONE: break;\
case HMC_CAS_equal_16B: STAT_CORE_EVENT(core_id, HMC_INST_COUNT_CAS_equal_16B);break;\
case HMC_CAS_zero_16B:  STAT_CORE_EVENT(core_id, HMC_INST_COUNT_CAS_zero_16B);break;\
case HMC_CAS_greater_16B:STAT_CORE_EVENT(core_id, HMC_INST_COUNT_CAS_greater_16B);break;\
case HMC_CAS_less_16B:  STAT_CORE_EVENT(core_id, HMC_INST_COUNT_CAS_less_16B);break;\
case HMC_ADD_16B:       STAT_CORE_EVENT(core_id, HMC_INST_COUNT_ADD_16B);break;\
case HMC_ADD_8B:        STAT_CORE_EVENT(core_id, HMC_INST_COUNT_ADD_8B);break;\
case HMC_ADD_DUAL:      STAT_CORE_EVENT(core_id, HMC_INST_COUNT_ADD_DUAL);break;\
case HMC_SWAP:          STAT_CORE_EVENT(core_id, HMC_INST_COUNT_SWAP);break;\
case HMC_BIT_WR:        STAT_CORE_EVENT(core_id, HMC_INST_COUNT_BIT_WR);break;\
case HMC_AND:           STAT_CORE_EVENT(core_id, HMC_INST_COUNT_AND);break;\
case HMC_NAND:          STAT_CORE_EVENT(core_id, HMC_INST_COUNT_NAND);break;\
case HMC_OR:            STAT_CORE_EVENT(core_id, HMC_INST_COUNT_OR);break;\
case HMC_XOR:           STAT_CORE_EVENT(core_id, HMC_INST_COUNT_XOR);break;\
case HMC_FP_ADD:        STAT_CORE_EVENT(core_id, HMC_INST_COUNT_FP_ADD);break;\
case HMC_COMP_greater:  STAT_CORE_EVENT(core_id, HMC_INST_COUNT_COMP_greater);break;\
case HMC_COMP_less:     STAT_CORE_EVENT(core_id, HMC_INST_COUNT_COMP_less);break;\
case HMC_COMP_equal:    STAT_CORE_EVENT(core_id, HMC_INST_COUNT_COMP_equal);break;\
case NUM_HMC_TYPES: break;}


/*
std::string HMC_Type_str[]=
{
    "HMC_NONE", // not really a HMC inst
    "HMC_CAS_equal_16B",
    "HMC_CAS_zero_16B",
    "HMC_CAS_greater_16B",
    "HMC_CAS_less_16B",
    "HMC_ADD_16B",
    "HMC_ADD_8B",
    "HMC_ADD_DUAL",
    "HMC_SWAP",
    "HMC_BIT_WR",
    "HMC_AND",
    "HMC_NAND",
    "HMC_OR",
    "HMC_XOR",
    "HMC_FP_ADD",
    "HMC_COMP_greater",
    "HMC_COMP_less",
    "HMC_COMP_equal",
    ""          //last element must be empty
};
*/
#endif
