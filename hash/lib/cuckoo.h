#pragma once

#include <vector>
#include <cstdint>
#include <math.h>
#include <iostream>
#include "hash.h"

class Cuckoo{
public:
    uint32_t size;
    //to store the left bits
    std::vector<uint64_t> l_bits;
    //the right bits to find the location to store r_bits
    uint64_t *r_bits;
    //wether a location of the table is occupied
    bool *status;

    //the mask to get the left bits and right bits
    uint32_t mask_len;
    uint32_t mask;

    //number of hash functions
    uint32_t hash_num;
    uint8_t hash_bits;
    // std::hash<uint32_t> uint32_hash;
    //the locations for different elements and hash functions
    // std::vector<std::vector<uint32_t>> locations;
    //record the current hash index
    std::vector<uint64_t> hash_indices;

    uint32_t left_bits_length;

public:
    Cuckoo(uint32_t size, uint32_t hash_num, uint32_t element_length){
        this->size = size;
        mask_len = (uint32_t) log2(size);
        mask = (1 << mask_len) - 1;

        l_bits.resize(size, -1);
        r_bits = new uint64_t[size]{};
        status = new bool[size]{false};

        this -> hash_num = hash_num;
        this -> hash_bits = (uint8_t) ceil(log2(hash_num));
        // locations.resize(size, std::vector<uint32_t>(hash_num, 0));
        hash_indices.resize(size, 0);

        //the bit length of each element
        this -> left_bits_length = element_length-mask_len;
    }
    ~Cuckoo(){
        delete []r_bits;
        delete []status;
    }

    uint32_t get_location(uint64_t element, uint32_t hash_index){
        // uint32_t mask = (1<<mask_len)-1;
        uint32_t r_part = element&mask;
        uint32_t l_part = element>>mask_len;
        //revserve hash_bits for the hash index
        // l_part <<= hash_bits;
        // l_part |= hash_index;
        // return (uint32_hash(l_part|hash_index)^r_part) % size;
        return (hasher((hash_index<<left_bits_length)|l_part)^r_part) % size;
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
    uint32_t insert(uint64_t element){
        // uint32_t mask = (1<<mask_len)-1;

        int max_iterations = size;
        uint32_t location = get_location(element, 0);
        int count = 0;

        // uint32_t evicted_element = 0;
        uint32_t r_part = element&mask, l_part = element>>mask_len;
        // l_part <<= hash_bits;
        uint32_t hash_index = 0;

        uint32_t tmp_l_part = 0, tmp_r_part = 0;
        uint32_t tmp_hash_index = 0;

        while(status[location] == true && count<max_iterations){
            //prepare the next hash for the evicted element
            tmp_l_part = l_bits[location];
            tmp_r_part = r_bits[location];
            tmp_hash_index = (hash_indices[location]+1)%hash_num;

            l_bits[location] = l_part;
            r_bits[location] = r_part;
            hash_indices[location] = hash_index;
            // std::cout<<location<<", "<<l_part<<", "<<r_part<<hash_index<<std::endl;

            l_part = tmp_l_part;
            r_part = tmp_r_part;
            hash_index = tmp_hash_index;
            
            // location = (uint32_hash(l_part|hash_index)^r_part)%size;
            location = (hasher((hash_index<<left_bits_length)|l_part)^r_part)%size;

            count++;
        }
        //cycle occurs, return failure
        if(count == max_iterations){
            return -1;
        }
        l_bits[location] = l_part;
        r_bits[location] = r_part;
        hash_indices[location] = hash_index;
        status[location] = true;

        return location;
    }

    bool insert_set(const std::vector<uint64_t> &set){

        for(uint32_t i=0; i<set.size(); i++){
            if(insert(set[i]) == -1){
                return false;
            }
        }
        for(uint32_t i=0; i<size; i++){
            //add the hash index to the left bits
            // std::cout<<l_bits[i]<<","<<(uint32_t)hash_indices[i]<<": ";
            l_bits[i] =  (hash_indices[i]<<left_bits_length) | l_bits[i];
            // std::cout<<l_bits[i]<<", "<<(int)hash_indices[i]<<", "<<l_bits[i]<<" ;  "<<std::endl;
        }

        return true;
    }

    uint32_t query(const uint64_t &element){
        // uint32_t mask = (1<<mask_len)-1;
        uint32_t r_part = element&mask, l_part = element>>mask_len;
        l_part <<= hash_bits;
        uint32_t location = 0;
        for(uint32_t i=0; i<hash_num; i++){
            // location = (uint32_hash(l_part|i)^r_part)%size;
            // location = (hasher(l_part|i)^r_part)%size;
            location = (hasher((i<<left_bits_length)|l_part)^r_part)%size;

            // std::cout<<":"<<location<<", "<<status[location]<<", "<<(int)hash_indices[location]<<", "<<l_bits[location]<<", "<<(l_part|i)<<std::endl;
            if(status[location] == true && l_bits[location]==(l_part|i)){ //after insert_set, we have already contatenated the hash index i
                return location;
            }
        }
        //cannot find the element
        return -1;
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

    bool remove(const uint64_t &element){
        uint32_t location = query(element);
        if(location == -1){
            return false;
        }
        status[location] = false;
        return true;
    }
};