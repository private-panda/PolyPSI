#include "m_hash.h"

uint32_t hash_function(uint32_t element, uint8_t i){
    // uint32_t hash_value = 0;
    uint32_t primes[3] = {4126878211, 2658490831, 3168826777};
    std::hash<uint32_t> int_hash;
    return int_hash(element^primes[i]);
}