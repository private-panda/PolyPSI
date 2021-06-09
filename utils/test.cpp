#include <bitset>
#include <iostream>
#include <vector>
#include <NTL/ZZ.h>
#include <algorithm>
#include <chrono>
#include <random>
#include "sets.h"

using namespace std;
struct Entry{
    uint32_t left_bits;
    uint32_t right_bits;
    uint8_t hash_index;
    bool status;
    
    Entry(){};
    Entry(const Entry &entry){
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
//the first method for combination
//from https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
void comb(int N, int K)
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    int count = 0;
    // print integers and permute bitmask
    do {
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            // if (bitmask[i]) std::cout << " " << i;
            if(bitmask[i]){
                count++;
            }
        }
        // std::cout << std::endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    cout<<"count: "<<count<<endl;
}
//the second method for combination
template<typename T>
void pretty_print(const T& c, int combo)
{
    int n = c.size();
    for (int i = 0; i < n; ++i) {
        if ((combo >> i) & 1)
            cout << c[i] << ' ';
    }
    cout << endl;
}

template<typename T>
void combo(const T& c, int k)
{
    int n = c.size();
    int combo = (1 << k) - 1;       // k bit sets
    while (combo < 1<<n) {

        pretty_print(c, combo);

        int x = combo & -combo;
        int y = combo + x;
        int z = (combo & ~y);
        combo = z / x;
        combo >>= 1;
        combo |= y;
    }
}
int main(){
    std::vector<bool> vec{1, 0, 1, 1, 0}; 
    void *vec_ptr = nullptr;
    std::vector<int> vec1{1, 0, 1, 1, 0}; 
    vec_ptr = vec1.data();

    std::bitset<64> bs;
    std::cout<<sizeof(bs)<<std::endl;

    // NTL::ZZ seed;
    // RandomLen(seed, 32);
    // SetSeed(seed);

    // NTL::ZZ rnd1 = NTL::GenPrime_ZZ(32, 80);
    // NTL::ZZ rnd2 = NTL::GenPrime_ZZ(32, 80);
    // NTL::ZZ rnd3 = NTL::GenPrime_ZZ(32, 80);
    // std::cout<<rnd1<<", "<<rnd2<<", "<<rnd3<<endl;

    vec.resize(10);
    for(auto e: vec){
        cout<<e<<", ";
    }
    cout<<endl;

    uint32_t a = -1;
    cout<<hex<<a<<endl;
    a >>= 5;
    cout<<hex<<a<<endl;

    vector<Entry> entry_vec(10);
    entry_vec.resize(20);
    cout<<dec<<entry_vec.size()<<endl;

    chrono::microseconds d(0);
    chrono::high_resolution_clock::time_point t0, t1;
    t0 = chrono::high_resolution_clock::now();
    comb(4,2);
    t1 = chrono::high_resolution_clock::now();
    d = chrono::duration_cast<chrono::microseconds>(t1-t0);
    cout<<"time for combination: "<<d.count()<<endl;

    // NTL::ZZX zzx(NTL::ZZ(1L));
    // for(auto zz: zzx){
    //     cout<<zz<<", ";
    // }
    // long p=11;
    // long x = 0, y=3, z=-3;
    // x -= (p-1); y-=(p-1); z-=(p-1);
    // x %= p; y %= p; z %= p;
    // cout<<x<<", "<<y<<", "<<z<<endl;
    // cout<<(x+2*p-2)%p <<", "<<(y+2*p-2)%p <<", "<<(z+2*p-2)%p<<endl;
    // cout<< (x-(x+2*p-2)%p)%p <<", "<<(y-(y+2*p-2)%p)%p<<","<<(z-(z+2*p-2)%p)%p<<endl;

    ulong p=11;
    ulong x = 0, y=3, z=10;
    // x -= (p-1); y-=(p-1); z-=(p-1);
    // x %= p; y %= p; z %= p;
    cout<<x<<", "<<y<<", "<<z<<endl;
    cout<<(x+p-1)%p <<", "<<(y+p-1)%p <<", "<<(z+p-1)%p<<endl;
    cout<< (x-(x+p-1)%p)%p <<", "<<(y-(y+p-1)%p)%p<<","<<(z-(z+p-1)%p)%p<<endl;

    ulong t = long(-10);
    cout<<t<<endl;

    random_device rd;
    auto r1 = rd();
    auto r2 = rd();
    cout<<r1<<", "<<r2<<endl;
    cout<<(r1+r2-2*r1*r2)<<", "<<(r1^r2)<<endl;

    cout<<-3%11<<endl;
    cout<<3%11<<endl;

    std::vector<uint64_t> my_vec;
    set_generate(my_vec, 1<<20);
    cout<<my_vec.size()<<endl;

    return 0;
}