#include "bitVec.h"
#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main(){
    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // uint64_t modulus_pow_3 = (uint64_t)plain_modulus*plain_modulus*plain_modulus;
    std::uniform_int_distribution<uint32_t> distrib;
    uint32_t modulus = 65537;
    
    vector<uint32_t> vec(32*17);
    for(int i=0; i<vec.size(); i++){
        vec[i] = distrib(gen) % modulus;
    }
    
    bitVector bit_vec(vec.size(), modulus);
    for(int i=0; i<vec.size(); i++){
        bit_vec.put(vec[i]);
    }
    // cout<<endl;
    // for(int i=0; i<10; i++){
    //     std::cout<<bit_vec.buffer_remain_bits[i]<<", ";
    // }
    // cout<<endl;
    for(int i=0; i<bit_vec.size; i++){
        cout<<vec[i]<<", ";
        if((i+1)%bit_vec.bits_unit_per_uint32_t == 0){
            cout<<hex<<bit_vec.buffer_remain_bits[i/bit_vec.bits_unit_per_uint32_t]<<std::endl;
        }
    }
    bit_vec.set_count_zero();
    
    vector<vector<uint32_t>> vec_restore(32, vector<uint32_t>(17));

    bit_vec.get(vec_restore);
    for(int i=0; i<32; i++){
        for(int j=0; j<17; j++){
            if(vec_restore[i][j] != vec[i*17+j]){
                cout<<"find unequal: " <<i<<","<<j<<", "<< vec_restore[i][j] <<", "<<vec[i*17+j]<<std::endl;
                break;
            }
        }
    }

    bit_vec.free_buffer();
    return 0;
}