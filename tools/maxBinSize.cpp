#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <vector>

using namespace std;
using namespace NTL;

RR comb(int n, int i){
    RR ans(1.0);
    RR::SetPrecision(1000);
    int j = 0;
    while(j<i){
        ans *= RR(n-j);
        j += 1;
    }
    while(j>1){
        ans /= RR(j);
        j -= 1;
    }
    return ans;
}

RR max_bin_prob(int n, double b, int max_b){
    
    RR ans(0.0);
    RR tmp(0.0);
    // cout<<ans<<endl;
    //decimal precision
    RR::SetPrecision(1000);
    // ans += (comb(n, i)*(1/b**i)*((1-1/b)**(n-i)))
    RR for_comb(1.0);
    RR for_b_i(1.0);
    //to get rid of using division that losses precision, we use a vector here
    vector<RR> for_1_b_n_i(max_b, RR(1.0));
    RR tmp_1_b_n_i = RR(1.0)-RR(1.0)/RR(b);
    for(int i=0; i<n-max_b+1; i++){
        for_1_b_n_i[0] *= tmp_1_b_n_i;
    }
    
    for(int i=1; i<max_b; i++){
        for_1_b_n_i[i] = for_1_b_n_i[i-1]*tmp_1_b_n_i;
    }
    // cout<<for_1_b_n_i<<", "<<endl;
    int i = 1;
    while(i <= max_b){
        ans += (for_comb * for_b_i * for_1_b_n_i[max_b-i]);
        // cout<<ans<<", "<<endl;
        // cout<<(for_comb * for_b_i * for_1_b_n_i[max_b-i])<<endl;
        // tmp *= for_b_i;
        // tmp *= for_1_b_n_i;

        for_comb *= RR(n-i+1)/RR(i);
        for_b_i /= RR(b);
        // for_1_b_n_i /= RR(1.0-1.0/b);

        // cout<<endl;
        i += 1;
        // cout<<ans<<", "<<endl;
        if(i == max_b/4){
            cout<<"\tprogress ... 1/4"<<endl;
            cout<<ans<<endl;
        }
        if(i == max_b/2){
            cout<<"\tprogress ... 1/2"<<endl;
            cout<<ans<<endl;
        }
        if(i == max_b*3/4){
            cout<<"\tprogress ... 3/4"<<endl;
            cout<<ans<<endl;

        }
        if(i == max_b*9/10){
            cout<<"\tprogress ... 9/10"<<endl;
            cout<<ans<<endl;

        }
    }
    cout<<ans<<endl;

    ans = b*(1-ans);
    cout<<"The ans: "<<ans<<endl;
    if(ans > 1){
        cout<<"Warning: the probability is greater than 1."<<endl;
    }

    return -log(ans)/log(2);
}

RR max_bin_prob1(int n, int b, int max_b){
    cout<<"starting ..."<<endl;
    RR ans(0.0);
    RR tmp(0.0);
    // cout<<ans<<endl;
    //decimal precision
    RR::SetPrecision(100000);
    // ans += (comb(n, i)*(1/b**i)*((1-1/b)**(n-i)))
    int i=0;
    while(i < max_b){
        tmp = comb(n, i);
        for(int j=0; j<i; j++){
            tmp /= RR(b);
        }
        for(int j=0; j<n-i; j++){
            tmp *= RR(1-1.0/b);
            // if(j<10)
            // cout<<tmp<<", ";
        }
        cout<<tmp<<endl;
        ans += tmp;
        i += 1;
        // cout<<ans<<", "<<endl;
        if(i == max_b/4){
            cout<<"\t"<<ans;
            cout<<"\tprogress ... 1/4"<<endl;
        }
        if(i == max_b/2){
            cout<<"\t"<<ans;
            cout<<"\tprogress ... 1/2"<<endl;
        }
        if(i == max_b*3/4){
            cout<<"\t"<<ans;
            cout<<"\tprogress ... 3/4"<<endl;
        }
        if(i == max_b*9/10){
            cout<<"\t"<<ans;
            cout<<"\tprogress ... 9/10"<<endl;
        }
    }
    // cout<<ans<<endl;
    ans = b*(1-ans);
    cout<<"The ans is: "<<ans<<endl;

    return -log(ans)/log(2);
}

//inputs: the bit length of the tag, the partition number, the number of hash functions
RR tag_false_positive_probability(int tag_len, int alpha, int k){
    RR::SetPrecision(1000);

    RR p = power2_RR(tag_len);
    p = 1/p;
    RR ans = power(RR(1.0)-power(1-p, alpha), k);

    return -log(ans)/log(2);

}

int main(){
    // int n=(1<<25)+(1<<15);
    // int max_bin_size = 716;
    // int n=3*(1<<13);

    // int b=1024;
    // int max_b=17412;
    // double root_len=(long long)1 << 36; //take the hash table length m=2^13, the actual stored root length is 31-13=18
    // double xx = ((long long)1<<31)
    // cout<<"xx: "<<xx<<endl;
    
    RR::SetOutputPrecision(20);
    // cout<<comb(100, 10)<<endl;
    // cout<<n<<", "<<root_len<<", "<<max_bin_size<<endl;
    // // cout<<"ffffff"<<endl;
    // RR max_prob = max_bin_prob(n, root_len, max_bin_size);
    // std::cout<<"the max bin probability is: "<<max_prob<<std::endl;
    // // cout<<"\t"<<max_prob<<endl;

    // int tag_len = 23; // we aim to find a balanced root_len-log2(m) = tal_len
    // int alpha = max_bin_size;
    // int k = 2;
    // std::cout<<"tag_len: "<<tag_len << ", "<< "alpha: "<<alpha<< ", " << "k: "<<k<<endl;
    // RR tag_fpp = tag_false_positive_probability(tag_len, alpha, k);
    // std::cout<<"tag_fpp: "<<tag_fpp<<std::endl;

    int n=14;
    // int n = 13247;
    int max_bin_size = 5;
    // int max_bin_size = ;

    long omega = 65537; //the number of bins
    RR max_prob = max_bin_prob(n, omega, max_bin_size);
    std::cout<<"max bin size: "<<max_bin_size<<std::endl;
    std::cout<<"the max_prob: "<<max_prob<<std::endl;
    return 0;
}