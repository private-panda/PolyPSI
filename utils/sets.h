#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <random>
#include <unordered_map>
#include <unordered_set>

void set_generate(std::vector<uint64_t> &set, uint32_t n){
    // set.resize(n);
    std::random_device rd;
    uint64_t element = 0;
    std::unordered_set<uint64_t> tmp_set;
    for(int i=0; i<n; i++){
        tmp_set.emplace(rd());
    }
    //there maybe duplicated elements, so we need to generate recursively
    int tmp_size = tmp_set.size();
    while(tmp_size != n){
        tmp_set.emplace(rd());
        if(tmp_set.size() != tmp_size){
            tmp_size++;
        }
    }
    for(auto& element: tmp_set){
        set.emplace_back(element);
    }
}
std::vector<uint64_t> get_intersection(std::vector<uint64_t> &set1, std::vector<uint64_t> &set2){
    // std::unordered_set<uint32_t> uset1(set1.begin(), set1.end());
    std::unordered_set<uint64_t> uset2(set2.begin(), set2.end());
    std::vector<uint64_t> intersection;
    for(auto &element: set1){
        if(uset2.find(element) != uset2.end()){
            intersection.push_back(element);
        }
    }
    return intersection;
}
void set_generate(std::vector<std::array<uint64_t, 2>> &set, uint32_t n){
    set.resize(n);
    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<uint64_t> distrib(0, (uint64_t)18446744073709551615); //18446744073709551615=2^64-1
    std::uniform_int_distribution<uint64_t> distrib; //18446744073709551615=2^64-1

    for(int i=0; i<n; i++){
        set[i][0] = distrib(gen);
        set[i][1] = distrib(gen);
    }
}

std::vector<std::array<uint64_t, 2>> get_intersection(std::vector<std::array<uint64_t, 2>> &set1, std::vector<std::array<uint64_t, 2>> &set2){
    // std::unordered_set<uint32_t> uset1(set1.begin(), set1.end());
    std::unordered_map<uint64_t, uint64_t> uset2;
    for(auto &element: set2){
        uset2[element[0]] = element[1];
    }
    if(uset2.size() != set2.size()){
        std::cout<<"there are duplicated keys"<<std::endl;
    }
    std::vector<std::array<uint64_t, 2>> intersection;
    for(auto &element: set1){
        if(uset2.find(element[0]) != uset2.end() && uset2[element[0]] == element[1]){
            intersection.push_back(element);
        }
    }
    return intersection;
}

