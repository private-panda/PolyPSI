
#include <vector>
#include <cstdint>
#include <math.h>
#include <iostream>
// #include "hash.h"
#include "cuckoo.h"
#include "../utils/sets.h"
#include <unistd.h>
#include "../tools/item.h"
#include <sys/time.h>
#include <thread>

// class Cuckoo{
// public:
//     uint32_t size;
//     //to store the left bits
//     std::vector<uint32_t> l_bits;
//     //the right bits to find the location to store r_bits
//     uint32_t *r_bits;
//     //wether a location of the table is occupied
//     bool *status;

//     //the mask to get the left bits and right bits
//     uint32_t mask_len;
//     uint32_t mask;

//     //number of hash functions
//     uint32_t hash_num;
//     // std::hash<uint32_t> uint32_hash;
//     //the locations for different elements and hash functions
//     // std::vector<std::vector<uint32_t>> locations;
//     //record the current hash index
//     std::vector<uint8_t> hash_indices;

// public:
//     Cuckoo(uint32_t size, uint32_t hash_num){
//         this->size = size;
//         mask_len = (uint32_t) log2(size);
//         mask = (1 << mask_len) - 1;

//         l_bits.resize(size);
//         r_bits = new uint32_t[size]{};
//         status = new bool[size]{false};

//         this -> hash_num = hash_num;
//         // locations.resize(size, std::vector<uint32_t>(hash_num, 0));
//         hash_indices.resize(size, 0);
//     }
//     ~Cuckoo(){
//         delete []r_bits;
//         delete []status;
//     }

//     uint32_t get_location(uint32_t element, uint32_t hash_index){
//         // uint32_t mask = (1<<mask_len)-1;
//         uint32_t r_part = element&mask;
//         uint32_t l_part = element>>mask_len;
//         //revserve 2 bits for the hash index
//         l_part <<= 2;
//         // l_part |= hash_index;
//         // return (uint32_hash(l_part|hash_index)^r_part) % size;
//         return (hasher(l_part|hash_index)^r_part) % size;
//     }
//     // void init_locations(const std::vector<uint32_t> &set){
//     //     uint32_t mask = (1<<mask_len)-1;
    
//     //     for(int i=0; i<set.size(); i++){
//     //         for(uint32_t j=0; j<hash_num; j++){
//     //             r_bits[i] = set[i]&mask;
//     //             l_bits[i] = set[i]>>mask_len;
//     //             //we need 2 bits for the hash function index in l_bits[i]
//     //             l_bits[i] <<= 2;
//     //             l_bits[i] |= j;
//     //             locations[i][j] = (uint32_hash(l_bits[i])^r_bits[i]) % size;
//     //         }
//     //     }
//     // }
//     uint32_t insert(uint32_t element){
//         // uint32_t mask = (1<<mask_len)-1;

//         int max_iterations = size;
//         uint32_t location = get_location(element, 0);
//         int count = 0;

//         uint32_t evicted_element = 0;
//         uint32_t r_part = element&mask, l_part = element>>mask_len;
//         l_part <<= 2;
//         uint32_t hash_index = 0;

//         uint32_t tmp_l_part = 0, tmp_r_part = 0;
//         uint32_t tmp_hash_index = 0;

//         while(status[location] == true && count<max_iterations){
//             //prepare the next hash for the evicted element
//             tmp_l_part = l_bits[location];
//             tmp_r_part = r_bits[location];
//             tmp_hash_index = (hash_indices[location]+1)%hash_num;

//             l_bits[location] = l_part;
//             r_bits[location] = r_part;
//             hash_indices[location] = hash_index;
//             // std::cout<<location<<", "<<l_part<<", "<<r_part<<hash_index<<std::endl;

//             l_part = tmp_l_part;
//             r_part = tmp_r_part;
//             hash_index = tmp_hash_index;
            
//             // location = (uint32_hash(l_part|hash_index)^r_part)%size;
//             location = (hasher(l_part|hash_index)^r_part)%size;

//             count++;
//         }
//         //cycle occurs, return failure
//         if(count == max_iterations){
//             return -1;
//         }
//         l_bits[location] = l_part;
//         r_bits[location] = r_part;
//         hash_indices[location] = hash_index;
//         status[location] = true;

//         return location;
//     }

//     bool insert_set(const std::vector<uint32_t> &set){

//         for(uint32_t i=0; i<set.size(); i++){
//             if(insert(set[i]) == -1){
//                 return false;
//             }
//         }
//         for(uint32_t i=0; i<set.size(); i++){
//             //add the hash index to the left bits
//             std::cout<<l_bits[i]<<","<<(uint32_t)hash_indices[i]<<": ";
//             l_bits[i] |= hash_indices[i];
//             std::cout<<l_bits[i]<<" ;  ";
//         }

//         return true;
//     }

//     uint32_t query(const uint32_t &element){
//         // uint32_t mask = (1<<mask_len)-1;
//         uint32_t r_part = element&mask, l_part = element>>mask_len;
//         l_part <<= 2;
//         uint32_t location = 0;
//         for(uint32_t i=0; i<hash_num; i++){
//             // location = (uint32_hash(l_part|i)^r_part)%size;
//             location = (hasher(l_part|i)^r_part)%size;
//             // std::cout<<":"<<location<<", "<<status[location]<<", "<<l_bits[location]<<", "<<(l_part|i)<<std::endl;
//             if(status[location] == true && l_bits[location]==(l_part|i)){
//                 return location;
//             }
//         }
//         //cannot find the element
//         return -1;
//     }

//     bool remove(const uint32_t &element){
//         uint32_t location = query(element);
//         if(location == -1){
//             return false;
//         }
//         status[location] = false;
//         return true;
//     }
// };
// class item{
// public:
//     uint32_t l;
//     uint32_t m;
//     uint32_t r;
// };
static double getMillies(timeval start, timeval end)
{
	long time1 = (start.tv_sec * 1000000) + (start.tv_usec );
	long time2 = (end.tv_sec * 1000000) + (end.tv_usec );

	return (double)(time2-time1)/1000;
}

void print_my_vec_one_more(std::vector<uint64_t> &vec){
    //declare a vector in the main, the n+1 element will be another value
    for(int i=0; i<vec.size()+1; i++){
        std::cout<<vec[i]<<", ";
    }
    std::cout<<std::endl;

    //declare a vector in the function, the n+1 element will be 0
    std::vector<uint64_t> my_vec(vec.size());
    for(int i=0; i<my_vec.size()+1; i++){
        std::cout<<my_vec[i]<<", ";
    }
    std::cout<<std::endl;

}

int main(){
    // int m=1<<13, n=m/2;
    int m=1<<13, n=5535;

    // int m=200, n= 100;
    // Cuckoo cuckoo(12, 2);
    // int sigma = 46;
    uint8_t hash_num = 3;
    // uint64_t modulus = 41299;
    // uint64_t modulus = 16411; //when m=8096, n=m/2, we have many faliures
    uint64_t modulus = 65537;
    
    Cuckoo cuckoo(m, hash_num, modulus);

    std::vector<item> set(n);

    int trails = 1<<19; // we need 6 minutes ~ 7 minutes for 1<<18
    int fails_count = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    // std::uniform_int_distribution<unsigned long> dist(0, modulus);

    timeval start, stop;
    // gettimeofday(&start, nullptr);
    // for(int j=0; j<n; j++){
    //     set[j].l = dist(gen);
    //     set[j].m = dist(gen);
    //     set[j].r = dist(gen);
    // }
    // gettimeofday(&stop, nullptr);
    // std::cout<<"By generate three times: time to generate "<<n<<" items: "<<getMillies(start, stop)<<std::endl;

    //is faster than the above generating random number three times
    uint64_t modulus_pow_3 = (uint64_t)modulus*modulus*modulus;
    std::uniform_int_distribution<unsigned long long> dist(0, modulus_pow_3-1);
    // std::uniform_int_distribution<unsigned long> dist(0, modulus);

    uint64_t rnd = 0;
    // gettimeofday(&start, nullptr);
    // for(int j=0; j<n; j++){
    //     rnd = dist1(gen);
    //     set[j].l = rnd%modulus; rnd /= modulus;
    //     set[j].m = rnd%modulus; rnd /= modulus;
    //     set[j].r = rnd;
    // }
    // gettimeofday(&stop, nullptr);
    // std::cout<<"By generate once: time to generate "<<n<<" items: "<<getMillies(start, stop)<<std::endl;

    auto fun = [&](int m_trails){
        for(int i=0; i<m_trails; i++){
            for(int j=0; j<n; j++){
                rnd = dist(gen);
                set[j].l = rnd%modulus; rnd /= modulus;
                set[j].m = rnd%modulus; rnd /= modulus;
                set[j].r = rnd%modulus;
                            // rnd = dist(gen);
                // set[j].l = dist(gen);
                // set[j].m = dist(gen);
                // set[j].r = dist(gen);
                // set[j] <<= 32;
                // set[j] ^=  (uint64_t)dist(gen);
            }
            // set_generate(set, n);//{1,2,3,4,5, 6, 100, 8};
            // std::cout<<set[0].l<<set[0].m<<set[0].r<<", ";
            bool flag = cuckoo.insert_set(set);
            // std::cout<<std::endl;
            // for(int j=0; j<cuckoo.size; j++){
            //     std::cout<<cuckoo.status[j]<<",";
            // }
            // std::cout<<std::endl;
            // for(int j=0; j<n; j++){
            //     for(int k=0; k<3; k++){
            //         std::cout<<g_hash(set[j], k, modulus)<<", ";
            //     }
            // }
            // std::cout<<std::endl;
            if(flag == false){
                fails_count++;
            }
            cuckoo.clear();
            // usleep(30000);
        }
    };
    // std::vector<std::thread> thread_vec;
    // thread_vec.emplace_back(std::thread(fun, trails));
    // thread_vec.emplace_back(std::thread(fun, trails));
    int thread_num = 1<<5;
    
    gettimeofday(&start, nullptr);
    for(int i=0; i<thread_num; i++){
        std::thread m_thread(fun, trails);
        m_thread.join();
    }
    gettimeofday(&stop, nullptr);
    std::cout<<"m, n: "<<m<<","<<n<<std::endl;
    std::cout<<"after "<<thread_num*trails<<" iterations, it takes "<<getMillies(start, stop)/1000<<" s"<<std::endl;
    std::cout<<"failed count: "<<fails_count<<std::endl;

    std::cout<<"("<<set[0].l<<","<<set[0].m<<","<<set[0].r<<")\t"<<std::endl;
    auto location = cuckoo.query(set[0]);
    auto location_1d = location.first;
    auto location_2d = location.second;
    std::cout<<location_1d<<","<<location_2d<<std::endl;
    uint32_t g_location0 = g_hash(set[0], 0, modulus) % m;
    uint32_t g_location1 = g_hash(set[0], 1, modulus) % m;
    uint32_t g_location2 = g_hash(set[0], 2, modulus) % m;
    std::cout<<"status for g0: "<<cuckoo.status[g_location0]<<std::endl;
    std::cout<<"("<<cuckoo.bins[g_location0].l<<","<<cuckoo.bins[g_location0].m<<","<<cuckoo.bins[g_location0].r<<")\t"<<std::endl;
    
    std::cout<<"status for g1: "<<cuckoo.status[g_location1]<<std::endl;
    std::cout<<"("<<cuckoo.bins[g_location1].l<<","<<cuckoo.bins[g_location1].m<<","<<cuckoo.bins[g_location1].r<<")\t"<<std::endl;

    std::cout<<"status for g2: "<<cuckoo.status[g_location2]<<std::endl;
    std::cout<<"("<<cuckoo.bins[g_location2].l<<","<<cuckoo.bins[g_location2].m<<","<<cuckoo.bins[g_location2].r<<")\t"<<std::endl;
    // cuckoo.bins[g_location0]
    
    // for(int i=0; i<cuckoo.size; i++){
    //     std::cout<<cuckoo.status[i]<<",";
    // }
    // std::cout<<std::endl;

    // int count = 0;
    // for(int i=0; i<cuckoo.size; i++){
    //     if(cuckoo.hash_indices[i]!=hash_num){
    //         count++;
    //     }
    //     std::cout<<(uint32_t)cuckoo.hash_indices[i]<<",";
    // }
    // std::cout<<std::endl;
    // std::cout<<count<<std::endl;
    // std::cout<<true<<", "<<false<<std::endl;
    // std::cout<<cuckoo.num_elements()<<std::endl;
    

    // for(auto e: set){
    //     auto location = cuckoo.query(e);
    //     if(location.first == -1 && location.second == -1){
    //         std::cout<<"("<<e.l<<","<<e.m<<","<<e.r<<")\t";
    //     }
    // }

    // set_generate(set, n);//{1,2,3,4,5, 6, 100, 8};
    // // for(auto e: set){
    // //     std::cout<<e<<", ";
    // // }
    // // std::cout<<std::endl;
    // bool flag = cuckoo.insert_set(set);
    // if(flag == false){
    //     std::cout<<"insert set failed!"<<std::endl;
    //     // return -1;
    // }
    // // for(auto e: cuckoo.bins){
    // //     std::cout<<e<<", ";
    // // }
    // // std::cout<<std::endl;
    // std::cout<<cuckoo.num_elements()<<std::endl;
    // // std::cout<<"cuckoo table: "<<std::endl;
    // // for(int i=0; i<cuckoo.size; i++){
    // //     std::cout<<cuckoo.l_bits[i]<<", "<<cuckoo.r_bits[i]<<", "<<(int)cuckoo.hash_indices[i]<<", "<<cuckoo.status[i]<<std::endl;
    // // }
    // // cuckoo.insert(9);
    // std::cout<<"start query the set: "<<std::endl;
    // int i=0;
    // int count = 0;
    // for(auto element: set){
    //     // std::cout<<i<<": ";
    //     if(cuckoo.query(element) == -1){
    //         std::cout<<"query "<<element<<" failed, "<<std::endl;
    //         count++;
    //         // return false;
    //     }else{
    //         // std::cout<<element<<", ";
    //     }
    //     // i++;
    // }
    // std::cout<<"Failure times: "<<count<<std::endl;
    // std::cout<<std::endl;
    // std::cout<<"start query an item: "<<std::endl;

    // uint32_t location = cuckoo.query(99);

    // if(location == -1){
    //     std::cout<<"query "<<10<<" failed"<<std::endl;
    // }else{
    //     std::cout<<location<<std::endl;
    //     std::cout<<cuckoo.status[10]<<", "<<cuckoo.l_bits[location]<<std::endl;
    //     std::cout<<" query "<<10<<" succeed"<<std::endl;
    // }
    // std::cout<<std::endl;

    // // bool *new_test = new bool[5]{false};
    // // for(int i=0; i<5; i++){
    // //     std::cout<<new_test[i]<<", ";
    // // }
    // std::cout<<std::endl;
    // std::hash<uint32_t> uint32_hash;
    // std::cout<<uint32_hash(100)<<std::endl;

    // std::cout<<sizeof(item)<<std::endl;

    // int ff = 10;
    // uint64_t vec[ff];// = {0};
    // for(int i=0; i<ff; i++){
    //     std::cout<<vec[i]<<", ";
    // }
    // std::cout<<std::endl;

    //conclusion, when the data types are different, (l-m_+modu) is different from modu+l-m_). Caution: divide in uint is dangerous
    // uint32_t l = 21034, m_ = 24611;
    // uint64_t modu = 65537;
    // std::cout<<(l-m_)<<", "<<(l-m_)%modu<<", "<<(l-m_+modu)<<", "<<(modu+l-m_)<<","<<(l-m_+modu)%modu<<std::endl;
    // uint64_t ll = l, mm = m_;
    // std::cout<<(ll-mm)<<", "<<(ll-mm)%modu<<", "<<(ll-mm+modu)<<", "<<(ll-mm+modu)%modu<<std::endl;
    // uint32_t modu1 = modu;
    // std::cout<<(l-m_)<<", "<<(l-m_)%modu1<<", "<<(l-m_+modu1)<<", "<<(l-m_+modu1)%modu1<<std::endl;

    // std::vector<uint64_t> my_vec = {1,2,3,4};
    // std::vector<uint64_t> my_vec1(4);
    
    // for(int i=0; i<5; i++){
    //     std::cout<<my_vec1[i]<<", ";
    // }
    // std::cout<<std::endl;

    // for(int i=0; i<4; i++){
    //     my_vec1[i] = my_vec[i];
    // }
    // std::cout<<my_vec[4]<<std::endl;
    // std::cout<<my_vec1[4]<<std::endl;

    // print_my_vec_one_more(my_vec);

    // uint32_t abc[3];
    // std::cout<<abc[0]<<", ";
    // abc[0] <<= 1;
    // std::cout<<abc[0]<<std::endl;

    // uint32_t num = 65486;
    // uint64_t tmp = 65486;
    // uint32_t m_modulus = 65537;
    // std::cout<<(tmp*num)%m_modulus<<", "<<(num*tmp)%m_modulus<<std::endl;
    // std::cout<<sizeof(tmp*num)<<", "<<sizeof(num*tmp)<<std::endl;
    return 0;
}