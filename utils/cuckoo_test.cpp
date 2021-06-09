#include "cuckoo.h"

using namespace std;

int main(){
    cout<<"start to initiate cuckoo hashing"<<endl;
    CuckooHashing ch(10, 3, 10);
    cout<<"iniate cuckoo hashing"<<endl;
    random_device rd;
    vector<uint32_t> vec(6);
    for(auto &e: vec){
        e = rd();
    }
    vec[0] = 3686967744;
    //before insertion
    cout<<"before insertion"<<endl;
    for(auto &e: vec){
        // const uint32_t t = e;
        cout<<e<<", ";
        bool flag = ch.insert(e);
        if(flag == false){
            cout<<endl<<"insert "<<e<<" failed"<<endl;
        }
    }
    
    cout<<endl;
    cout<<"before query"<<endl;
    for(auto &e: vec){
        uint32_t location = ch.query(e);
        if(location == -1){
            cout<<"query "<<e<<" failed"<<endl;
            continue;
        }
        // cout<<location<<", ";
        // cout<<ch.table[location].left_bits<<", "<<ch.table[location].right_bits<<", "<<(uint32_t)ch.table[location].hash_index<<":"<<endl;
    }
    cout<<endl;
    uint32_t location = ch.query(rd());
    cout<<hex<<location<<endl;

    double load_factor = ch.get_load_factor();
    cout<<load_factor<<endl;

    return 0;
}