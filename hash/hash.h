#pragma once

#include "./lib/BLAKE3/c/blake3.h"
// #include "blake3.h"
// #include <stdio.h>
#include <unistd.h>
#include "../tools/item.h"
#include <array>

// uint32_t hasher(uint32_t &element){
//       // Initialize the hasher.
//   blake3_hasher hasher;
//   blake3_hasher_init(&hasher);

//   blake3_hasher_update(&hasher, &element, 4);

//   uint32_t output=0;
//   blake3_hasher_finalize(&hasher, (uint8_t*)&output, sizeof(output));

//   return output;
// }

uint32_t hasher(uint64_t element){
      // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  blake3_hasher_update(&hasher, (uint8_t*)&element, 8);

  uint64_t output=0;
  blake3_hasher_finalize(&hasher, (uint8_t*)&output, sizeof(output));

  return output;
}

void hasher(std::array<uint64_t, 2> &element, std::array<uint64_t, 2> &output){
      // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  blake3_hasher_update(&hasher, (uint8_t*)element.data(), 16);

  blake3_hasher_finalize(&hasher, (uint8_t*)output.data(), 16);

}

uint32_t g_hash(item element, uint8_t hash_index, uint32_t modulus){
    uint64_t tmp_r = element.r; // cast from uint32_t to uint64_t
    uint32_t mul_r_l = (tmp_r * element.l)%modulus;
    if(hash_index == 0){
      //each slice, either l, m, or r, should meet l*r+m < (modulus-1)*(modulus-1)+(modulus-1) < 2**32, then we can set the data type 
      //as uint32_t and we do not need to do modulo for l*r. Otherwise, it needs data type uint64_t as the data type.
        return (mul_r_l + element.l + element.r)%modulus;
    }else if(hash_index == 1){
        return (mul_r_l + element.m + element.r + 1)%modulus;
    }else if(hash_index == 2){
      return (mul_r_l + element.r + 2)%modulus;
    }
    std::cout<<"the hash index is out of range"<<std::endl;
    return modulus;
}