#pragma once

#include <vector>
#include <cstdint>
#include "m_hash.h"
#include <random>

class CuckooHashing{
private:
    uint32_t size; // the number of elements the table contains
    uint32_t capacity; // the maximum number of elements, i.e., the size of hash table
    uint8_t hash_num; // the number of hash functions
    uint32_t max_iteration; // the maximum iterations for fear of loops
private:
    uint32_t table_bits_num; // for using permutation-based hashing
    uint8_t hash_bits_num; //show which hash function is used, and store in the location
    
public:
    struct Entry{
        uint32_t left_bits;
        uint32_t right_bits;
        uint8_t hash_index;
        bool status;

        Entry(){};
        Entry(const Entry &entry){ //const is very important
            this->left_bits = entry.left_bits;
            this->right_bits = entry.right_bits;
            this->hash_index = entry.hash_index;
            this->status = entry.status;
        };
        public:
        Entry& operator=(const Entry& entry){
            this->hash_index = entry.hash_index;
            this->left_bits = entry.left_bits;
            this->right_bits = entry.right_bits;
            this->status = entry.status;
            return *this;
        };

    };
// private:
public:
    std::vector<Entry> table;
public:
    std::vector<uint32_t> table_real; //in table, the hash function index and rights bits are separate, now concatenate them into an uint32_t
    // std::vector<uint32_t> table; // store the right bits and the hash function index
    // std::vector<uint32_t> table_original; //store the original element (i.e., the element before permutation-based hashing) for eviction
    // std::vector<bool> location_status; //if a location is set, set the its status as true; otherwise false.
public:
    CuckooHashing(uint32_t capacity, uint8_t hash_num, uint32_t max_iteration);
    // CuckooHashing(uint32_t capacity, uint8_t hash_num, uint32_t max_iteration){
    //     this->capacity = capacity;
    //     table.resize(capacity);
    //     table_real.resize(capacity);
    //     // location_status.resize(capacity); // the default value in each location is false
    //     this->size = 0;
    //     this->hash_num = hash_num; 
    //     this->hash_bits_num = ceil(log2(hash_num));
    //     this->max_iteration = max_iteration;

    //     // std::cout<<"initate cuckoo hashing successfully"<<std::endl;
    // }
    // CuckooHashing(const std::vector<uint32_t> &set, uint32_t size, uint8_t hash_num){
    //     this->size = size;
    //     table.resize(size);
    //     this->hash_num = hash_num;
    // }
    //insert an element into the table, permutation hashing is used
    bool insert(const uint32_t element);
    //query and return the location the element is stored
    uint32_t query(const uint32_t element);
    bool remove(const uint32_t element);
    
    //insert the set into the cuckoo hash table, permutation hashing is used
    bool insert_set(const std::vector<uint32_t> &set);

    //get the real table
    std::vector<uint64_t> get_real_table();

    //insert dummy element into the table
    void insert_dummies();

    //get the load factor of cuckoo hash table
    double get_load_factor();
};