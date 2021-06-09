#pragma once

#include <vector>
#include <cstdint>
#include <iostream>
#include "seal/seal.h"

//here W is the window, d is the bin size, d is the multiplicative depth
void get_window(std::vector<uint32_t> &W, const uint32_t &b, uint32_t &d){
    if(d == 0){
        for(int i=1; i<=b; i++){
            W.push_back(i);
        }
        return;
    }
    uint32_t l_b = 0; // the bit length of b
    uint32_t tmp_b = b;
    while(tmp_b){
        tmp_b >>= 1;
        l_b++;
    }
    if(l_b <= d){
        d = l_b -1;
    }
    // std::cout<<"l_b: "<<l_b<<std::endl;

    uint32_t prefix_b = 0, postfix_b = 0;
    uint32_t postfix_mask = (1<<(l_b-d-1))-1;
    postfix_b = b & postfix_mask;
    prefix_b = b >> (l_b-d-1);

    // std::cout<<"prefix: "<<prefix_b<<"\t"<<"postfix: "<<postfix_b<<std::endl;
    // std::cout<<"(1<<(d+1))-1): "<<(1<<(d+1))-1<<std::endl;

    uint32_t tmp_bound = 1<<(l_b-d-1);
    if(prefix_b == (1<<(d+1))-1){
        for(int i=1; i<=tmp_bound+postfix_b; i++){
            W.push_back(i);
            // std::cout<<i<<std::endl;
        }
    }else{
        for(int i=1; i<=tmp_bound; i++){
            W.push_back(i);
        }
    }

    for(int i=l_b-d; i<l_b-1; i++){
        W.push_back(1<<i);
    }
    
    //check the second bit of b is 1 or not to decide whether we need the the base power 1<<(l_b-1)
    if((b>>(l_b-2)) & 1){
        W.push_back(1<<(l_b-1));
    }
}

void de_windowing(std::vector<seal::Ciphertext> &ctxt_y_powers_full, const std::vector<uint32_t> &W, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys){
    uint32_t window_size = W.size();
    uint32_t b = ctxt_y_powers_full.size();
    
    int i=0;
    // std::cout<<"de-windowing: "<<window_size<<std::endl;
    while((i+1) == W[i]){
        // std::cout<<i<<", ";
        i++;
    }
    uint32_t next_start_pos = i;
    uint32_t base = 1;
    while(i != 1){
        base <<= 1;
        i >>= 1;
    }
    
    // std::cout<<"find the base"<<std::endl;

    uint32_t depth = window_size-next_start_pos;
    uint32_t max_base = (base<<depth);
    // std::cout<<"depth: "<<depth<<",  max_base: "<<max_base<<std::endl;
    if(max_base*2 <= b){
        evaluator.square(ctxt_y_powers_full[max_base-1], ctxt_y_powers_full[max_base*2-1]);
        evaluator.relinearize_inplace(ctxt_y_powers_full[max_base*2-1], relin_keys);
        // std::cout<<max_base*2<<", ";
        depth++;
        max_base <<= 1;
    }

    uint32_t next_base = base<<1;
    for(int i=0; i<depth; i++){
        for(int j=next_start_pos; j<next_base-1; j++){
            evaluator.multiply(ctxt_y_powers_full[base-1], ctxt_y_powers_full[j-base], ctxt_y_powers_full[j]);
            evaluator.relinearize_inplace(ctxt_y_powers_full[j], relin_keys);

            // std::cout<<j<<", ";
        }
        next_start_pos = next_base;
        next_base <<= 1;
    }

    // std::cout<<"next_start_pos: "<<next_start_pos<<std::endl;
    for(int i=next_start_pos; i<b; i++){
        evaluator.multiply(ctxt_y_powers_full[max_base-1], ctxt_y_powers_full[i-max_base], ctxt_y_powers_full[i]);
        evaluator.relinearize_inplace(ctxt_y_powers_full[i], relin_keys);
        // std::cout<<i<<", ";
    }

    // std::cout<<std::endl;
}