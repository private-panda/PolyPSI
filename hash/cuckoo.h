#pragma once

#include <vector>
#include <cstdint>
#include <math.h>
#include <iostream>
#include "hash.h"
#include "../tools/item.h"


class Cuckoo{
public:
    uint32_t size;
    uint32_t size_bit_length;
    uint32_t size_mask; //since the polymodulus is with the power of 2, we can use bitwise & instead of %

    //whether a location is occupied
    std::vector<bool> status;
    std::vector<item> bins;

    //number of hash functions
    uint8_t hash_num;
    // uint32_t hash_bits;

    //record the current hash index
    std::vector<uint8_t> hash_indices;

    uint32_t modulus;

    // uint32_t slice_length;
    // uint64_t slice_mask;


public:
    //plain_modulus to support field arithmetic operations
    Cuckoo(uint32_t size, uint8_t hash_num, uint32_t modulus){
        this->size = size;
        this->size_bit_length = (uint32_t)log2(size);
        this->size_mask = ((uint32_t)1<<size_bit_length)-1;
        this -> hash_num = hash_num;
        // hash_bits = (uint32_t)ceil(log2(hash_num));
        // this -> hash_bits = (uint8_t) ceil(log2(hash_num));
        // locations.resize(size, std::vector<uint32_t>(hash_num, 0));
        bins.resize(size);
        status.resize(size, false);
        hash_indices.resize(size, hash_num);
        
        this->modulus = modulus;

        // slice_length = (uint32_t)ceil(element_length*1.0/3);
        // slice_mask = ((uint64_t)1 << slice_length) - 1;
    }
    ~Cuckoo(){
    }


    // void init_locations(const std::vector<uint32_t> &set){
    //     uint32_t mask = (1<<mask_len)-1;
    
    //     for(int i=0; i<set.size(); i++){
    //         for(uint32_t j=0; j<hash_num; j++){
    //             r_bits[i] = set[i]&mask;
    //             l_bits[i] = set[i]>>mask_len;
    //             //we need 2 bits for the hash function index in l_bits[i]
    //             l_bits[i] <<= 2;
    //             l_bits[i] |= j;
    //             locations[i][j] = (uint32_hash(l_bits[i])^r_bits[i]) % size;
    //         }
    //     }
    // }
    uint32_t insert(item element){
        // uint32_t mask = (1<<mask_len)-1;

        int max_iterations = size;
        uint8_t hash_index = 0;
        // uint32_t location = g_hash(element, hash_index, modulus) % size;
        uint32_t location = g_hash(element, hash_index, modulus) & size_mask;

        int count = 0;

        item evicted_element;
        uint8_t evicted_hash_index = 0;

        while(status[location] == true && count<max_iterations){
            //prepare the next hash for the evicted element
            evicted_element = bins[location];
            evicted_hash_index = hash_indices[location];

            bins[location] = element;
            hash_indices[location] = hash_index;

            element = evicted_element;
            hash_index = (evicted_hash_index+1)%hash_num;
            location = g_hash(element, hash_index, modulus) & size_mask;

            count++;
        }
        //cycle occurs, return failure
        if(count == max_iterations){
            return -1;
        }

        bins[location] = element;
        hash_indices[location] = hash_index;
        status[location] = true;

        return location;
    }

    bool insert_set(const std::vector<item> &set){
        for(uint32_t i=0; i<set.size(); i++){
            if(insert(set[i]) == -1){
                return false;
            }
            // uint32_t location = insert(set[i]);
            // if(location != -1){
            //     std::cout<<status[location]<<", ";
            // }
        }
        // for(int j=0; j<size; j++){
        //     std::cout<<status[j]<<",";
        // }
        return true;
    }

    std::pair<uint32_t, uint32_t> query(const item &element){
        // uint32_t mask = (1<<mask_len)-1;
        // l_part <<= hash_bits;
        uint32_t location = 0, location_1d=0, location_2d=0;
        for(int i=0; i<hash_num; i++){
            location = g_hash(element, i, modulus);
            location_1d = location&size_mask; location >>= size_bit_length;
            location_2d = location;
            if(status[location_1d] == true && element == bins[location_1d]){
                return std::make_pair(location_1d, location_2d);
            }
        }
        //cannot find the element
        return std::make_pair(-1, -1);
    }
    
    void clear(){
        for(int i=0; i<size; i++){
            status[i] = false;
        }
    }

    uint32_t num_elements(){
        int count = 0;
        for(int i=0; i<size; i++){
            if(status[i] == true){
                count++;
            }
        }
        return count;
    }

    bool remove(const item &element){
        uint32_t location = query(element).first;
        if(location == -1){
            return false;
        }
        status[location] = false;
        return true;
    }
};