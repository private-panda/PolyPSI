#pragma once

#include <cstdint> 
#include <vector>
#include <bitset>
#include <iostream>
#include "defs.h"
#include <math.h>

// //to compress the element into a narrow domain
// uint32_t domain_hashing(std::bitset<sigma> element, uint32_t n1, uint32_t n2, uint32_t lambda){
//     // std::hash<std::vector<bool>> bits_hash;
//     uint32_t desired_length = 2*ceil(log2(n1+n2)) + lambda -1;
//     if(desired_length >= sigma){
//         return 0;
//     }

//     // auto image = bits_hash(element);
    

//     return 0;
// }

// std::pair<uint32_t, uint32_t> permutation_hashing(const uint32_t element, uint32_t size){
//     uint32_t left_bits_num = floor(log2(size));
//     uint32_t left_bits = element>>(sigma-left_bits_num);
//     uint32_t zeros_bit = (1ull<<31)-1;
//     zeros_bit >>= left_bits_num;
//     uint32_t right_bits = element & zeros_bit; //keep the right sigma-log(m) bits

//     std::hash<uint32_t> int_hash;
//     uint32_t location = left_bits ^ int_hash(right_bits);
//     return std::move(std::pair<uint32_t, uint32_t>{location, right_bits});
// }
//to get the hash value of an element by using hash function i
uint32_t hash_function(uint32_t element, uint8_t i);

std::vector<std::vector<uint32_t>> simple_hashing(std::vector<uint32_t> &set, uint32_t m, uint8_t hash_num, uint32_t bin_height){
    std::vector<std::vector<uint32_t>> result(m);
    uint32_t left_bits=0, right_bits=0;
    uint32_t hash_bits_num = (uint32_t)ceil(log2(hash_num));
    uint32_t table_bits_num = (uint32_t)floor(log2(m));
    uint32_t zeros_bit = uint32_t(-1)>>table_bits_num;
    uint32_t location = 0;
    
    for(auto &element: set){
        left_bits = element>>(sigma - table_bits_num);
        right_bits = element & zeros_bit;
        for(uint8_t j=0; j<hash_num; j++){
            location = left_bits ^ hash_function(right_bits, j);
            result[location%m].emplace_back((j<<(sigma-table_bits_num)) ^ right_bits);
        }
    }

    //insert dummy
    uint32_t dummy = 1ul << (sigma-table_bits_num+hash_bits_num); // 2^(sigma')
    for(uint i=0; i<result.size(); i++){
        while(result[i].size() < bin_height){
            result[i].push_back(dummy);
        }
    }
    return result;
}

