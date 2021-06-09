#include "windowing.h"
#include <iostream>

using namespace std;

int main(){
    uint32_t b = 8; 
    uint32_t d = 1;
    vector<uint32_t> W;
    get_window(W, b, d);

    for(int i=0; i<W.size(); i++){
        cout<<W[i]<<", ";
    }
    cout<<endl;

    return 0;
}