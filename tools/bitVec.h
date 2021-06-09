#pragma once

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

class bitVector{
public:
    uint32_t size;
    
    uint32_t modulus_bits_num;
    uint32_t modulus_bytes_len;
    uint32_t modulus_remain_bits_len;
    uint32_t modulus_remain_bits_mask;

    uint32_t buffer_bytes_size;
    uint16_t* buffer_bytes; //for both 65537 and 1638574, the byte length is 2, so we use uint16_t
    uint32_t buffer_remain_bits_size;
    uint32_t* buffer_remain_bits; //uint32_t to avoid data type conversion
    uint32_t tmp_remain_bits;

    uint32_t byte_count = 0;
    uint32_t bit_count = 0;
    uint32_t inter_byte_count = 0;

    uint32_t bits_unit_per_uint32_t;

public:
    bitVector(uint32_t size, const uint32_t &modulus){
        this->size = size;

        modulus_bits_num = (uint32_t)log2(modulus)+1; //since the plain_modulus must not be divided by 2, we add it by 1 directly.
        modulus_bytes_len = modulus_bits_num / 8;
        modulus_remain_bits_len = modulus_bits_num % 8; //for modulus=65537, remian_bits is 1; for modulus=163841, remain_bits is 2.
        modulus_remain_bits_mask = (1<<modulus_remain_bits_len)-1;

        buffer_bytes_size = size * modulus_bytes_len;
        buffer_bytes = (uint16_t *)malloc(buffer_bytes_size);
        buffer_remain_bits_size = size * modulus_remain_bits_len /8;
        buffer_remain_bits = (uint32_t *)malloc(buffer_remain_bits_size); //size must be divided by 32

        bits_unit_per_uint32_t = 32/modulus_remain_bits_len;
    }

    //to insert a field element
    void put(const uint32_t &c){
        tmp_remain_bits <<= modulus_remain_bits_len;
        tmp_remain_bits ^= (c & modulus_remain_bits_mask);
        inter_byte_count++;
        if(inter_byte_count == bits_unit_per_uint32_t){
            //reverse it to keep the bits in get() in order
            while(inter_byte_count){
                buffer_remain_bits[bit_count] <<= modulus_remain_bits_len;
                buffer_remain_bits[bit_count] ^= (tmp_remain_bits & modulus_remain_bits_mask);
                tmp_remain_bits >>= modulus_remain_bits_len;
                inter_byte_count--;
            }
            bit_count++;
        }
        buffer_bytes[byte_count] = (c >> modulus_remain_bits_len);
        byte_count++;
    }
    void put(const std::vector<std::vector<uint32_t>> &vec_2d){
        uint32_t range_1d = vec_2d.size();
        uint32_t range_2d = vec_2d[0].size();
        for(int i = 0; i < range_1d; i++){
            for(int j = 0; j < range_2d; j++){
                put(vec_2d[i][j]);
            }
        }
    }
    void get(std::vector<std::vector<uint32_t>> &vec_2d){
        uint32_t range_1d = vec_2d.size();
        uint32_t range_2d = vec_2d[0].size();
        for(int i=0; i<range_1d; i++){
            for(int j=0; j<range_2d; j++){
                vec_2d[i][j] = buffer_bytes[byte_count];
                byte_count++;
                vec_2d[i][j] <<= modulus_remain_bits_len;
                vec_2d[i][j] ^= (buffer_remain_bits[bit_count] & modulus_remain_bits_mask);
                buffer_remain_bits[bit_count] >>= modulus_remain_bits_len;
                inter_byte_count++;
                if(inter_byte_count == bits_unit_per_uint32_t){
                    inter_byte_count = 0;
                    bit_count++;
                }
            }
            // std::cout<<byte_count<<", "<<bit_count<<"; ";
        }   
    }
    void set_count_zero(){
        byte_count = 0;
        bit_count = 0;
        inter_byte_count = 0;
    }
    void free_buffer(){
        free((void *)buffer_bytes);
        free((void *)buffer_remain_bits);
    }
};