#include "cuckoo.h"

 CuckooHashing::CuckooHashing(uint32_t capacity, uint8_t hash_num, uint32_t max_iteration){
    this->capacity = capacity;
    table.resize(capacity);
    table_real.resize(capacity);
    // location_status.resize(capacity); // the default value in each location is false
    this->size = 0;
    this->hash_num = hash_num; 
    this->hash_bits_num = ceil(log2(hash_num));
    this->table_bits_num = floor(log2(capacity));
    this->max_iteration = max_iteration;
    // std::cout<<"initate cuckoo hashing successfully"<<std::endl;
}

bool CuckooHashing::insert(const uint32_t element){

    Entry entry;
    entry.left_bits = (element>>(sigma-table_bits_num)); //get the left bits
    // std::cout<<"find left bits: "<<sigma<<", "<<table_bits_num<<", "<<(sigma-table_bits_num)<<", "<<(element>>(sigma-table_bits_num))<<std::endl;
    uint32_t zeros_bit = -1; //(1ull<<32)-1; zero_bits is initiated as 0xffffffff
    zeros_bit >>= table_bits_num;
    entry.right_bits = (element & zeros_bit); //keep the right sigma-log(m) bits
    entry.status = false;
    entry.hash_index = 0;
    // std::cout<<"insert "<<element<<": ";
    // std::cout<<entry.left_bits<<", "<<entry.right_bits<<std::endl;

    // std::cout<<"before eviction: "<<std::endl;
    uint32_t location = 0; //the location the right bits and hash function index will be stored
    uint32_t start_hash_index = 0, hash_index_bound = 0;
    Entry evicted_entry = entry;
    evicted_entry.status = true;
    for(uint32_t i=0; i<max_iteration; i++){
        hash_index_bound = start_hash_index+hash_num;
        for(start_hash_index; start_hash_index<hash_index_bound; start_hash_index++){
            location = evicted_entry.left_bits ^ hash_function(evicted_entry.right_bits, start_hash_index%hash_num);
            location %= capacity; //the location should be within capacity
            if(table[location].status == false){
                // right_bits <<= hash_bits_num;
                // right_bits ^= j; //concatenate the right bits with the hash function index
                // table[location].left_bits = left_bits;
                // table[location].right_bits = right_bits;
                // table[location].hash_index = j;
                table[location] = evicted_entry;
                table[location].hash_index = start_hash_index%hash_num;
                table[location].status = true;
                // std::cout<<"location: "<<location<<", "<<element<<std::endl;
                // std::cout<<table[location].left_bits<<", "<<table[location].right_bits<<", "<<(uint32_t)table[location].hash_index<<", "<<table[location].status<<std::endl;
                size++; //increase the size of the table by 1
                return true;
            }
        }
        //if no empty place found, start to evict the entry about the first hash location; then iteratively evict elements
        location = evicted_entry.left_bits ^ hash_function(evicted_entry.right_bits, (start_hash_index)%hash_num);
        location %= capacity;
        // evicted_entry = table[location];
        // std::swap(evicted_entry, table[location]);
        Entry tmp = evicted_entry; //swap(evicted_entry, table[location]);
        evicted_entry = table[location];
        table[location] = tmp;

        // table[location].status= true;//set the status of this location as true
        table[location].hash_index = (start_hash_index)%hash_num;

        start_hash_index = (evicted_entry.hash_index+1)%hash_num; // the next hash function to be explored
    }
    
    //cannot find a proper location
    return false;
}

uint32_t CuckooHashing::query(const uint32_t element){
    uint32_t left_bits = element>>(sigma - table_bits_num); //get the left bits
    uint32_t zeros_bit = -1;
    zeros_bit >>= table_bits_num;
    uint32_t right_bits = element & zeros_bit; //keep the right sigma-log(m) bits

    uint32_t location = 0;
    for(uint8_t j=0; j<hash_num; j++){
        location = left_bits ^ hash_function(right_bits, j);
        location %= capacity;
        std::cout<<"location: "<<location<<std::endl;
        std::cout<<(uint32_t)table[location].hash_index<<", "<<table[location].right_bits<<": ";
        std::cout<<(uint32_t)j<<", "<<right_bits<<std::endl;
        if(table[location].hash_index==j && table[location].right_bits==right_bits){
            return location; //return the location of element
        }
    }
    return -1; //return 0xfffffff, the table capacity should not be this value
}

bool CuckooHashing::remove(const uint32_t element){
    uint32_t left_bits = element>>(sigma - table_bits_num); //get the left bits
    uint32_t zeros_bit = -1;
    zeros_bit >>= table_bits_num;
    uint32_t right_bits = element & zeros_bit; //keep the right sigma-log(m) bits

    uint32_t location = 0;
    for(uint32_t j=0; j<hash_num; j++){
        location = left_bits ^ hash_function(right_bits, j);
        location %= capacity;
        if(table[location].hash_index==j && table[location].right_bits==right_bits){
            table[location].status = false; //put this location as unfilled
            size--;
            return true; //delete success
        }
    }

    return false; // fail to find this element in the table, thus cannot remove it
}

bool CuckooHashing::insert_set(const std::vector<uint32_t> &set){
    bool insert_flag = true;
    for(auto &element: set){
        insert_flag = insert(element);
        if(insert_flag == false){
            std::cerr<<"insert "<<element<<" failed"<<std::endl;
            return false;
        }
    }
    return insert_flag;
}

// void CuckooHashing::insert_dummies(){
//     std::random_device rd;
//     for(auto &entry: table){
//         if(entry.status == false){
//             entry.hash_index = rd()%hash_num;
//             entry.right_bits = rd()%(1<<(sigma-table_bits_num));
//         }
//     }   
// }
std::vector<uint64_t> CuckooHashing::get_real_table(){
    std::vector<uint64_t> result(capacity);
    uint32_t dummy1 = (1ul << (sigma-table_bits_num+hash_bits_num))-1;
    for(uint32_t i=0; i<capacity; i++){
        if(table[i].status == true){
            result[i] = (table[i].hash_index<<(sigma-table_bits_num)) ^ table[i].right_bits;
            // table_real[i] = (table[i].hash_index<<(sigma-table_bits_num)) ^ table[i].right_bits; //contanetate the hash index with the rights bits
        }else{
            //if not set, put dummy elements into the table
            // table_real[i] = dummy1; //2^(sigma')-1
            result[i] = dummy1;
        }
    }
    return result;
}

double CuckooHashing::get_load_factor(){
    return size*1.0/capacity;
}