#pragma once

#include <vector>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include "seal/seal.h"

//here W is the window, d is the bin size, ell is the circuit depth
void get_window(std::vector<uint32_t> &W, const uint32_t &b, uint32_t &ell){
    if(ell == 0){
        for(int i=1; i<=b; i++){
            W.push_back(i);
        }
        return;
    }
    uint32_t log2_b = 0;
    uint32_t power2_ell = 1<<ell;
    uint32_t tmp_b = b;
    while(tmp_b){
        log2_b++;
        tmp_b >>= 1;
    }
    log2_b--;
    uint32_t upper_bound_j = log2_b/ell;
    uint32_t tmp_ans = 0;
    for(int i=1; i <= power2_ell-1; i++){
        for(int j=0; j<= upper_bound_j; j++){
            tmp_ans = i*(1<<(ell*j));
            if(tmp_ans <= b){
                W.emplace_back(tmp_ans);
            }
        }
    }
    sort(W.begin(), W.end());
}
void de_windowing(std::vector<seal::Ciphertext> &ctxt_y_powers_full, std::vector<uint32_t> &W, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys){
    W.push_back(0);
    uint32_t cur_next = 0;
    for(int i=0; i<ctxt_y_powers_full.size(); i++){
        // std::cout<<i<<", "<< W[cur_next] <<", "<<cur_next<<", ";
        if((i+1) == W[cur_next]){
            cur_next++;
        }else{
            // std::cout<<i<<":("<<W[cur_next-1]-1<<","<<i - W[cur_next-1]<<"), ";
            evaluator.multiply(ctxt_y_powers_full[W[cur_next-1]-1], ctxt_y_powers_full[i - W[cur_next-1]], ctxt_y_powers_full[i]);
            evaluator.relinearize_inplace(ctxt_y_powers_full[i], relin_keys);
            // ctxt_y_powers_full[i] = ctxt_y_powers_full[W[cur_next-1]-1] + ctxt_y_powers_full[i - W[cur_next-1]];
        }
        // std::cout<<i<<":("<<W[cur_next-1]-1<<","<<i - W[cur_next-1]<<"), "<<std::endl;
    }
    // std::cout<<std::endl;
}
