#include "seal/seal.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include "hash/cuckoo.h"
// #include "tools/poly_coeffs.h"
#include "tools/x_coeffs.h"
#include "./utils/sets.h"
#include "./tools/bitVec.h"
#include "./tools/windowing.h"
#include <set>

#include "utils/client.h"
#include "utils/server.h"

#include <NTL/ZZ.h>

#include <fstream>


//in this version, we expand each item from 32-bit to 128-bit by using virtual bloom filter and labels
using namespace std;
using namespace seal;

static double getMillies(timeval start, timeval end)
{
	long time1 = (start.tv_sec * 1000000) + (start.tv_usec );
	long time2 = (end.tv_sec * 1000000) + (end.tv_usec );

	return (double)(time2-time1)/1000;
}

class poly_parms{
public:
    int nx;
    int ny;
    uint32_t modulus;
    uint32_t omega;
    // uint32_t nu;
    uint32_t B;
    uint32_t alpha;
    uint32_t b;
    uint32_t d; // the multiplicative depth
    std::vector<uint32_t> W; // the window

    int hash_num;
};

void print_poly_parms(const poly_parms &p_parms){
    std::cout<<"nx, ny: "<<p_parms.nx<<", "<<p_parms.ny<<std::endl;
    std::cout<<"modulus: "<<p_parms.modulus<<"; "<<"omega: "<<p_parms.omega<<std::endl;
    std::cout<<"B: "<<p_parms.B<<"; alpha: "<<p_parms.alpha<<std::endl;
    std::cout<<"d: "<<p_parms.d<<std::endl;
    std::cout<<"The window: "<<std::endl;
    for(int i=0; i<p_parms.W.size(); i++){
        std::cout<<p_parms.W[i]<<", ";
    }
    std::cout<<std::endl;
}

void read_random_set(std::vector<std::array<uint64_t, 2>> &set, std::string file_name, const uint32_t set_size){
    set.resize(set_size);
    auto myfile = std::fstream(file_name, std::ios::in | std::ios::binary);
    if(!myfile.is_open()){
        std::cout<<"cannot open the file"<<std::endl;
        exit(0);
    }
    uint32_t bytes_uint64_t = sizeof(uint64_t);

    for(int i=0; i<set_size; i++){
        myfile.read((char*)&set[i][0], bytes_uint64_t);
        myfile.read((char*)&set[i][1], bytes_uint64_t);
    }
    myfile.close();
}

void poly_low_sender(std::vector<std::array<uint64_t, 2>> &X, const poly_parms &p_parms){

    //1. In the preprocessing phase, we do DH-OPRF to process the items. Reference: Fast Secure Computation of Set Intersection, Stanislaw Jarecki and Xiaomin Liu
    //First, sample a random key in G(p).
    uint32_t ny = p_parms.ny, nx = p_parms.nx;
    // NTL::ZZ p; NTL::GenPrime(p, 2048);

    
    uint32_t plain_modulus = p_parms.modulus;
    uint32_t omega = p_parms.omega; // the hash table length
    // uint32_t ell= sigma_vbf/3; //the slice bit length, 2^ell should be no less than omega
    //to do batching, omega should be no less than 4096 in the default polynomial modulus security setting

    uint32_t nu = (uint32_t)ceil(plain_modulus*1.0/omega); // the number of sub-bins in each bin

    uint32_t B = p_parms.B, alpha = p_parms.alpha;
    uint32_t b = p_parms.B/p_parms.alpha; //alpha divides B
    
    std::vector<uint32_t> W = p_parms.W; //the multiplicative depth

    uint32_t hash_num = p_parms.hash_num;

    timeval start, stop;

/*for the sender setup*******************************************************************/
    //for different set size 3*nx*2={2^16, 2^20, 2^24}, its corresponding different hash index length 2^ell={8192, 16384, 65536}, 
    // we get the max bin size B (with the statistical security parameter lambda=40). 
    //for different max bin size B, after slicing, the maximum number of items with the same root x^(r) alpha is:
    //alpha means we need to partition each sub-bin into at least alpha sub-sub-bins to keep x^(r) distinct.

    // uint64_t P_l_coeffs[omega][nu][B];
    // uint64_t P_m_coeffs[omega][nu][B];
    //declare a 3d vector
    gettimeofday(&start, nullptr);
    // std::cout<<"before P_l_coffs"<<std::endl;

    std::vector<std::vector<std::vector<uint32_t>>> P_l_coeffs(omega, std::vector<std::vector<uint32_t>>(nu, std::vector<uint32_t>(B, 0)));
    
    {
        // std::vector<item> X_vbf(nx*2);
        // item x_simple_hashing[omega][nu][B];//to hash
        std::vector<std::vector<std::vector<item>>> x_simple_hashing(omega, std::vector<std::vector<item>>(nu, std::vector<item>(B)));

        {
            // uint32_t hash_sub_bin_rear[omega][nu] = {0}; // initiate as 0, to record the rear of each sub-bin
            std::vector<std::vector<uint32_t>> hash_sub_bin_rear(omega, std::vector<uint32_t>(nu, 0));
            
            uint32_t omega_mask = omega - 1;
            uint32_t omega_bit_length = (uint32_t)log2(omega);

            uint32_t plain_modulus_bit_length = (uint32_t)log2(plain_modulus)+1; //since the modulus is a prime, we directly add 1
            uint32_t plain_modulus_bit_mask = (1<<plain_modulus_bit_length) - 1;

            std::array<uint64_t, 2> hash_output;
            // uint32_t vbf_index = 0;
            uint32_t rear = 0;
            item x_vbf_item;

            // uint32_t zero_count = 0;
            // uint32_t one_count = 0;
            // std::cout<<"before simple hashing: "<<nx<<std::endl;
            for(int i=0; i<nx; i++){
                // if(i%(1<<10)==0)
                // std::cout<<i<<",";
                hasher(X[i], hash_output); //do VBF
                uint32_t g_hash_index = 0, g_hash_1d_index = 0, g_hash_2d_index = 0;
                for(int j=0; j<2; j++){
                    //do slicing
                    // hash_output[j] &= 8*plain_modulus*plain_modulus*plain_modulus;
                    // x_vbf_item.r = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus; hash_output[j] /= plain_modulus; 
                    // x_vbf_item.m = hash_output[j]%plain_modulus; hash_output[j] /= plain_modulus; 
                    // x_vbf_item.l = hash_output[j]%plain_modulus;
                    x_vbf_item.r = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus; hash_output[j] >>= plain_modulus_bit_length; 
                    x_vbf_item.m = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus; hash_output[j] >>= plain_modulus_bit_length; 
                    x_vbf_item.l = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus;
                    // std::cout<<X_vbf[vbf_index].r<<","<<X_vbf[vbf_index].m<<","<<X_vbf[vbf_index].l<<"::";

                    for(int k=0; k<hash_num; k++){
                        g_hash_index = g_hash(x_vbf_item, k, plain_modulus);
                        g_hash_1d_index = g_hash_index & omega_mask; g_hash_index >>= omega_bit_length;
                        g_hash_2d_index = g_hash_index;
                        

                        rear = hash_sub_bin_rear[g_hash_1d_index][g_hash_2d_index];
                        if(rear == B){
                            std::cout<<"bin overflow"<<std::endl;
                            exit(0);
                        }
                        x_simple_hashing[g_hash_1d_index][g_hash_2d_index][rear] = x_vbf_item;
                        hash_sub_bin_rear[g_hash_1d_index][g_hash_2d_index] ++;
                    }
                }
            }

            // std::cout<<"zero count: "<<zero_count<<std::endl;
            // std::cout<<"one count: "<<one_count<<std::endl;

        // std::cout<<"finishing simple hashing"<<std::endl;
        //insert random dummies for each sub-bin
            // std::random_device rd;
            std::mt19937 gen(2948917759); //Standard mersenne_twister_engine seeded with rd()
            uint64_t modulus_pow_3 = (uint64_t)plain_modulus*plain_modulus*plain_modulus;
            std::uniform_int_distribution<uint64_t> distrib(0, modulus_pow_3 - 1);
            std::vector<item> dummies(b);
            uint64_t rnd = 0;
            for(int i=0; i<b; i++){
                rnd = distrib(gen);
                dummies[i].r = rnd%plain_modulus; rnd /= plain_modulus;
                dummies[i].m = rnd%plain_modulus; rnd /= plain_modulus;
                if(rnd == 0){
                    dummies[i].l = 1; // the dummy of the cuckoo hashing is 0, then the tag for the interpolation polynomial cannot be 0
                }else{
                    dummies[i].l = rnd;
                }
            }

            for(int i=0; i<omega; i++){
                for(int j=0; j<nu; j++){
                    rear = hash_sub_bin_rear[i][j];
                    for(int t=rear; t<B; t++){
                        x_simple_hashing[i][j][t] = dummies[t%b];
                    }
                }
            }

            
        }

        // std::cout<<"finishing padding"<<std::endl;

        //partition each sub-bin into alpha distinct sub-sub-bins
        {
            std::set<uint32_t> tmp_set;
            uint32_t tmp_size = 0;
            uint32_t tmp_count = 0;
            for(int i=0; i<omega; i++){
                for(int j=0; j<nu; j++){
                    for(int k=0; k<alpha-1; k++){
                        for(int t=0; t<b; t++){
                            tmp_set.emplace(x_simple_hashing[i][j][k*b+t].r);
                        }
                        tmp_size = tmp_set.size();
                        tmp_set.clear();
                        if(tmp_size < b){ //there are duplicated roots in this sub-sub-bin
                            uint64_t next_start = (k+1)*b; //start to find from next sub-sub-bin
                            for(int t=0; t<b && tmp_size<b; t++){
                                if(tmp_set.find(x_simple_hashing[i][j][k*b+t].r) != tmp_set.end()){ //it is the duplicated root
                                    //find a distinct root in the later sub-sub-sins and swap it
                                    while(next_start<B){
                                        if(tmp_set.find(x_simple_hashing[i][j][next_start].r) == tmp_set.end()){
                                            std::swap(x_simple_hashing[i][j][k*b+t], x_simple_hashing[i][j][next_start]);
                                            break;
                                        }
                                        next_start++;
                                    }
                                    tmp_size++;
                                }
                                tmp_set.emplace(x_simple_hashing[i][j][k*b+t].r);
                            }
                        }
                    }
                }
            }
        }


        polynomial_from_points(P_l_coeffs, x_simple_hashing, b, plain_modulus);

        gettimeofday(&stop, nullptr);
        std::cout<<"finishing setup, takes time (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;

    }

    gettimeofday(&start, nullptr);

    Server server;
    
    seal::EncryptionParameters  parms;
    server.channel.cRecv(parms);
    // std::cout<<"receive parms"<<std::endl;
    // sender_receive_data_size += parms.load(parms_stream);
    seal::SEALContext context(parms);

    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    RelinKeys relin_keys;
    
    server.channel.cRecv(context, relin_keys);
    std::cout<<"relinear key size: "<<relin_keys.size()<<std::endl;
    // sender_receive_data_size += relin_keys.load(context, rlk_stream);

    seal::Ciphertext ctxt_y_l;
    server.channel.cRecv(context, ctxt_y_l);
    // server.channel.cRecv(context, ctxt_y_m);

    uint32_t window_size = W.size();
    std::vector<seal::Ciphertext> ctxt_y_rs(b-1);
    for(int i=0; i<window_size; i++){
        server.channel.cRecv(context, ctxt_y_rs[W[i]-1]);
        // std::cout<<"receive ctxt y^"<<W[i]<<std::endl;
    }
    // std::cout<<"before dewindowing"<<std::endl;
    de_windowing(ctxt_y_rs, W, evaluator, relin_keys);
    std::cout<<"finish receiving y_powers and de-windowing"<<std::endl;
    // sender_receive_data_size += ctxt_y.load(context, data_stream);


    
    seal::Ciphertext ctxt_P_l;
    seal::Ciphertext ctxt_tmp;

    std::vector<uint64_t> rnd_mask(omega);
    std::random_device rd;
    seal::Plaintext ptxt_rnd;

    bool send_flag = false;
    seal::Ciphertext ctxt_send_tmp;

    seal::Plaintext ptxt_P_l;
    std::vector<uint64_t> tmp_P_l_coeff(omega);
    for(int j=0; j<nu; j++){
        for(int k=0; k<alpha; k++){
            //generate a random mask
            // for(int i=0; i<omega; i++){
            //     rnd_mask[i] = rd();
            // }
            // batch_encoder.encode(rnd_mask, ptxt_rnd);

            // Ciphertext ctxt_P_l, ctxt_P_m;
            //handle P_l
            ctxt_P_l = ctxt_y_l;

            for(int i=0; i<omega; i++){
                rnd_mask[i] = rd();

                tmp_P_l_coeff[i] = P_l_coeffs[i][j][k*b];
            }
            batch_encoder.encode(rnd_mask, ptxt_rnd);
            batch_encoder.encode(tmp_P_l_coeff, ptxt_P_l);
            evaluator.add_plain_inplace(ctxt_P_l, ptxt_P_l); //y*c1+c0

            for(int t=1; t<b; t++){
                for(int i=0; i<omega; i++){
                    tmp_P_l_coeff[i] = P_l_coeffs[i][j][k*b+t];
                }
                batch_encoder.encode(tmp_P_l_coeff, ptxt_P_l);
                evaluator.multiply_plain(ctxt_y_rs[t-1], ptxt_P_l, ctxt_tmp); //y*c1

                evaluator.relinearize_inplace(ctxt_tmp, relin_keys);
                evaluator.add_inplace(ctxt_P_l, ctxt_tmp);
            }

            
            if(send_flag){
                evaluator.multiply_inplace(ctxt_P_l, ctxt_send_tmp);
                evaluator.relinearize_inplace(ctxt_P_l, relin_keys);

                auto data_context = context.first_context_data();
                while(data_context->next_context_data()){
                    evaluator.mod_switch_to_next_inplace(ctxt_P_l);
                    data_context = data_context->next_context_data();
                }

                server.channel.cSend(ctxt_P_l);
                send_flag = false;
            }else{
                ctxt_send_tmp = ctxt_P_l;
                send_flag = true;
            }
            // server.channel.cSend(ctxt_P_m);
        }
    }

    gettimeofday(&stop, nullptr);
    std::cout<<"sender computes the encrypted solution (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;

    std::cout<<"Communication (MB): S->R: "<<server.channel.getDataSizeSent()<<"\t";
    std::cout<<"R->S: "<<server.channel.getDataSizeReceived()<<std::endl;
}

std::vector<std::array<uint64_t, 2>> poly_low_receiver(std::vector<std::array<uint64_t, 2>>& Y, const poly_parms &p_parms){
    uint32_t ny = p_parms.ny, nx = p_parms.nx;
    uint32_t plain_modulus = p_parms.modulus;
    uint32_t omega = p_parms.omega; // the hash table length

    // std::cout<<"plain_modulus&omega: "<<plain_modulus<<", "<<omega<<std::endl;

    uint32_t nu = (uint32_t)ceil(plain_modulus*1.0/omega); // the number of sub-bins in each bin

    uint32_t B = p_parms.B, alpha = p_parms.alpha;
    // std::cout<<alpha<<std::endl;
    uint32_t b = p_parms.b; //alpha divides B

    std::vector<uint32_t> W = p_parms.W;

    uint32_t hash_num = p_parms.hash_num;
    
    timeval start, stop;
    gettimeofday(&start, nullptr);

    std::vector<item> Y_vbf(ny*2);
    Cuckoo y_cuckoo(omega, hash_num, plain_modulus);

    {
        uint32_t omega_mask = omega - 1;
        uint32_t omega_bit_length = (uint32_t)log2(omega);

        uint32_t plain_modulus_bit_length = (uint32_t)log2(plain_modulus)+1; //since the modulus is a prime, we directly add 1
        uint32_t plain_modulus_bit_mask = (1<<plain_modulus_bit_length) - 1;

        std::array<uint64_t, 2> hash_output;
        int vbf_index = 0;
        for(int i=0; i<ny; i++){
            hasher(Y[i], hash_output); //VBF
            for(int j=0; j<2; j++){
                //do slicing
                Y_vbf[vbf_index].r = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus; hash_output[j] >>= plain_modulus_bit_length; 
                Y_vbf[vbf_index].m = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus; hash_output[j] >>= plain_modulus_bit_length; 
                Y_vbf[vbf_index].l = (hash_output[j]&plain_modulus_bit_mask) % plain_modulus;
                
                vbf_index++;
            }
        }
        
        bool flag = y_cuckoo.insert_set(Y_vbf);
        if(flag == false){
            std::cout<<"cuckoo hashing failure, exit(-1)"<<std::endl;
            exit(-1);
        }
    }
    std::cout<<"finish cuckoo hashing"<<std::endl;
    //fhe operations
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = omega;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    // if(nx == (1<<20) && ny == (1<<11)){
    //     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 32, 32})); //for ny=2^11, nx=2^20
    // }else if(nx == (1<<20) && ny == (1<<10)){
    //     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 25, 25, 25})); //for ny=2^10, nx=2^20
    // }else if(nx == (1<<20) && ny == (1<<12)){
        // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 32, 32, 32, 32})); //for ny=2^12, nx=2^20
    // }else if(nx == (1<<24) && ny == (1<<12)){

    // }
    if(nx == (1<<16)){
        if(ny == (1<<10)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 24, 24, 24}));
            // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree, seal::sec_level_type::tc128));


        }else if(ny == (1<<11)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 25, 25, 25, 25})); 

        }else if(ny == (1<<12)){
            // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 32, 32, 32, 32})); 
            // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 26, 27, 27})); 
            // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 32, 32})); 
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 26, 26})); 



        }

        if(ny == 5535){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 26, 26})); 
        }else if(ny == 11041){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 26, 26})); 
        }
    }else if(nx == (1<<20)){
        if(ny == (1<<10)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 24, 24, 25}));
        }else if(ny == (1<<11)){
            // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 32, 32})); 
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 25, 25, 26, 26})); 

        }else if(ny == (1<<12)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 27, 27, 27})); 

        }

        if(ny == 5535){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 27, 27, 27})); 
        }else if(ny == 11041){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 27, 27, 27})); 
        }

    }else if(nx == (1<<24)){
        if(ny == (1<<10)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 24, 24, 25}));
        }else if(ny == (1<<11)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {25, 25, 26, 26, 26}));
        }else if(ny == (1<<12)){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {26, 26, 26, 27, 27, 27}));
        }

        if(ny == 5535){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {32, 32, 32, 33, 33, 33})); 

        }else if(ny == 11041){
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {33, 33, 33, 33, 34, 34})); 
        }
    }

    std::cout<<"The coeff modulus: (";
    for(int i=0; i<parms.coeff_modulus().size(); i++){
        std::cout<<parms.coeff_modulus()[i].value()<<", "<<parms.coeff_modulus()[i].uint64_count()<<".  ";
    }
    std::cout<<")"<<std::endl;
    
    parms.set_plain_modulus(plain_modulus);

    // receiver_send_data_size += parms.save(parms_stream); //save the encryption parameters to help the sender get the context
    SEALContext context(parms);
    std::cout << "Parameter validation (success): " << context.parameter_error_message() << std::endl;

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    // std::cout<<"generate the secret key"<<std::endl;
    // PublicKey public_key;
    // keygen.create_public_key(public_key);
    // RelinKeys relin_keys;
    Serializable<RelinKeys> relin_keys = keygen.create_relin_keys();
        // std::cout<<"create relinear keys"<<std::endl;
        // receiver_send_data_size += relin_keys.save(rlk_stream); //save the relinear keys for the sender to relinear the ciphertext after multiplication
        // Serializable<RelinKeys> relin_keys = keygen.create_relin_keys();
        // std::cout<<"after relinear key"<<std::endl;

    Encryptor encryptor(context, secret_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    BatchEncoder batch_encoder(context);

    
    // auto my_location = y_cuckoo.query(Y_vbf[0]);
    // std::cout<<"the first stored vbf item: ("<<y_cuckoo.bins[my_location.first].l<<", "<<y_cuckoo.bins[my_location.first].m<<", "
    //             <<y_cuckoo.bins[my_location.first].r<<")"<<std::endl;
    // std::cout<<"after FHE preparation"<<std::endl;
    uint32_t window_size = W.size();

    std::vector<uint64_t> y_powers_r(omega), y_l(omega), y_m(omega);
    for(int i=0; i<omega; i++){
        y_powers_r[i] = y_cuckoo.bins[i].r;
        y_l[i] = plain_modulus - y_cuckoo.bins[i].l;
    }
    Plaintext ptxt_y_l, ptxt_y_m;
    batch_encoder.encode(y_l, ptxt_y_l);
    seal::Serializable<seal::Ciphertext> ctxt_y_seria_l = encryptor.encrypt_symmetric(ptxt_y_l);

    // std::cout<<y_powers_r[my_location.first]<<", ";
    Plaintext ptxt_y_r;
    batch_encoder.encode(y_powers_r, ptxt_y_r);
    
    seal::Serializable<seal::Ciphertext> ctxt_y_seria = encryptor.encrypt_symmetric(ptxt_y_r);

    gettimeofday(&stop, nullptr);
    std::cout<<"receiver setup time (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;


    gettimeofday(&start, nullptr);

    Client client;
    client.connect();

    client.channel.cSend(parms);
    // std::cout<<"send parms"<<std::endl;
    client.channel.cSend(relin_keys);
    client.channel.cSend(ctxt_y_seria_l);
    // client.channel.cSend(ctxt_y_seria_m);

    client.channel.cSend(ctxt_y_seria);
    // std::cout<<"send ctxt y^1"<<std::endl;

    for(int i=1; i<window_size; i++){
        if(W[i] == W[i-1]+1){
            for(int j=0; j<omega; j++){
                y_powers_r[j] = (y_powers_r[j] * y_cuckoo.bins[j].r) % plain_modulus;
            }
        }else if(W[i] == W[i-1]*2){
            for(int j=0; j<omega; j++){
                y_powers_r[j] = (y_powers_r[j] * y_powers_r[j]) % plain_modulus;
            }
        }else{//the case all the bits of the prefix are '1', the current power is a power of 2
            for(int j=0; j<omega; j++){
                for(int k=W[i-1]; k<W[i]; k++){
                    y_powers_r[j] = (y_powers_r[j] * y_cuckoo.bins[j].r) % plain_modulus;
                }
            }
        }
        // std::cout<<y_powers_r[my_location.first]<<", ";
        batch_encoder.encode(y_powers_r, ptxt_y_r);
        ctxt_y_seria = encryptor.encrypt_symmetric(ptxt_y_r);
        client.channel.cSend(ctxt_y_seria);
        // std::cout<<"send ctxt y^"<<i+1<<std::endl;
    }
    // std::cout<<std::endl;

    bool compress_ctxt_flag = true;
    if(compress_ctxt_flag){
        alpha /= 2;
    }
    std::vector<seal::Ciphertext> result_ciphertexts(nu*alpha);
    uint32_t bin_index = 0;

    bool recv_flag = false;
    
    for(int j=0; j<nu; j++){
        for(int k=0; k<alpha; k++){
            //receive the corresponding answer
            client.channel.cRecv(context, result_ciphertexts[bin_index]);

            // client.channel.cRecv(context, ctxt_P_m);

            std::cout<<decryptor.invariant_noise_budget(result_ciphertexts[bin_index])<<", ";

            bin_index++;
        }
    }
    std::cout<<std::endl;
    gettimeofday(&stop, nullptr);
    std::cout<<"receiver online time (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;

    gettimeofday(&start, nullptr);
    std::vector<std::array<uint64_t, 2>> intersection;
    // alpha /= 2; //have compressed the answer ciphertext by 2
    std::vector<std::vector<uint64_t>> ans_table_P_l(nu*alpha, std::vector<uint64_t>(omega));
    // std::vector<std::vector<uint64_t>> ans_table_P_m(nu*alpha, std::vector<uint64_t>(omega));

    // uint32_t bin_index = 0;
    Plaintext ptxt_P_l;
    // Ciphertext ctxt_P_l;
    bin_index = 0;
    for(int j=0; j<nu; j++){
        for(int k=0; k<alpha; k++){
            //receive the corresponding answer
            // client.channel.cRecv(context, ctxt_P_l);
            // // client.channel.cRecv(context, ctxt_P_m);

            // std::cout<<decryptor.invariant_noise_budget(ctxt_P_l)<<", ";
            // receiver_receive_data_size += ctxt_P_l.load(context, data_stream);
            // receiver_receive_data_size += ctxt_P_m.load(context, data_stream);
            decryptor.decrypt(result_ciphertexts[bin_index], ptxt_P_l);
            // decryptor.decrypt(ctxt_P_m, ptxt_P_m);

            batch_encoder.decode(ptxt_P_l, ans_table_P_l[bin_index]);
            // batch_encoder.decode(ptxt_P_m, ans_table_P_m[bin_index]);

            bin_index++;
        }
    }
    
    // gettimeofday(&stop, nullptr);
    // std::cout<<"receiver online time (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;
    
    // gettimeofday(&start, nullptr);

    for(int i=0; i<Y_vbf.size(); i += 2){
        std::pair<uint32_t, uint32_t> location1 = y_cuckoo.query(Y_vbf[i]);
        bool flag1 = false;
        // for(int k=0; k<alpha; k++){
        //     if(ans_table_P_l[alpha*location1.second+k][location1.first] == 0 && ans_table_P_m[alpha*location1.second+k][location1.first] == 0){
        //         flag1 = true;
        //     }
        // }
        for(int k=0; k<alpha; k++){
            if(ans_table_P_l[alpha*location1.second+k][location1.first] == 0){
                flag1 = true;
            }
        }

        std::pair<uint32_t, uint32_t> location2 = y_cuckoo.query(Y_vbf[i+1]);
        bool flag2 = false;
        // for(int k=0; k<alpha; k++){
        //     if(ans_table_P_l[alpha*location2.second+k][location2.first] == 0 && ans_table_P_m[alpha*location2.second+k][location2.first] == 0){
        //         flag2 = true;
        //     }
        // }
        for(int k=0; k<alpha; k++){
            if(ans_table_P_l[alpha*location2.second+k][location2.first] == 0){
                flag2 = true;
            }
        }

        if(flag1 && flag2){
            intersection.emplace_back(Y[i/2]);
        }
    }

    gettimeofday(&stop, nullptr);
    std::cout<<"receiver query time (s): "<<getMillies(start, stop)*1.0/1000<<std::endl;

    std::cout<<"communication (MB), S->R: "<<client.channel.getDataSizeReceived()<<"\t";
    std::cout<<"R->S: "<<client.channel.getDataSizeSent()<<std::endl;

    return intersection;
}


int main(int argc, char **argv)
{
    uint32_t role = atoi(argv[1]);
    int log_nx = atoi(argv[2]);
    int log_ny = atoi(argv[3]);

    poly_parms p_parms;
    
    p_parms.nx = 1<<log_nx;
    p_parms.ny = 1<<log_ny;
    if(log_ny == 5535){
        p_parms.ny = 5535;
    }else if(log_ny == 11041){
        p_parms.ny = 11041;
    }
    
    if(log_nx == 28){
        p_parms.modulus = 163841;
        p_parms.B = 10692; //10692 = 2^2 * 3^5 * 11
        //alpha_min = 10;
        p_parms.alpha = 44;
    }else if(log_nx == 24){
        p_parms.modulus = 65537;
        p_parms.B = 1880; //1880 = 2^3 * 5 * 47
        p_parms.alpha = 40;
    }else if(log_nx == 20){
        p_parms.modulus = 65537;
        p_parms.B = 162; 
        p_parms.alpha = 6;

        // alpha_min = 6;
    }else if(log_nx == 16){
        p_parms.modulus = 65537;
        p_parms.B = 36; 

        p_parms.alpha = 6;
        // alpha_min = 5;
    }
    p_parms.b = p_parms.B/p_parms.alpha;

    std::cout<<"logx ,y :"<<log_nx<<", "<<log_ny<<std::endl;

    //to do batching, omega should be no less than 4096 in the default polynomial modulus security setting
    if(log_ny == 12){
        p_parms.omega = (1<<14);
        p_parms.d = 2;
    }else if(log_ny == 11){
        p_parms.omega = (1<<13);
        p_parms.d = 2;
    }else if(log_ny == 10){
        p_parms.omega = (1<<12);
        p_parms.d = 0;
    }else if(log_ny == 8){
        p_parms.omega = (1<<12); //enable key switching and reserve enough noise budgets
        p_parms.d = 0;
    }

    if(log_ny == 5535){
        p_parms.omega = (1<<14);
        p_parms.d = 2;
    }else if(log_ny == 11041){
        p_parms.omega = (1<<15);
        p_parms.d = 2;
    }
    

    get_window(p_parms.W, p_parms.b-1, p_parms.d); // we aim to get y^1...(b-1)

    p_parms.hash_num = 3;

    print_poly_parms(p_parms);

    bool file_based = false;
    if(role == 0){
        std::vector<std::array<uint64_t, 2>> X;
        if(file_based){
            std::string file_name("../dataset/sender.binary");
            read_random_set(X, file_name, p_parms.nx);
        }else{
            set_generate(X, p_parms.nx);
        }
        X[0][0] = 12345; X[0][1] = 54321;
        X[1][0] = 11111; X[1][1] = 22222; 

        std::cout<<"The size of X: "<<X.size()<<std::endl;
        poly_low_sender(X, p_parms);
    }else if(role == 1){
        std::vector<std::array<uint64_t, 2>> Y;
        if(file_based){
            std::string file_name("../dataset/receiver.binary");
            read_random_set(Y, file_name, p_parms.ny);
        }else{
            set_generate(Y, p_parms.ny);
        }
        Y[0][0] = 12345; Y[0][1] = 54321;
        Y[1][0] = 11111; Y[1][1] = 22222; 

        auto get_intersection = poly_low_receiver(Y, p_parms);
        std::cout<<"Intersection get from HE: "<<get_intersection.size()<<std::endl;
        for(auto& element: get_intersection){
            std::cout<<element[0]<<element[1]<<", ";
        }
        std::cout<<std::endl;
    }else if(role == 2){
        //do offline psi
        std::vector<std::array<uint64_t, 2>> X, Y;
        set_generate(Y, p_parms.ny);
        set_generate(X, p_parms.nx);
        Y[0][0] = 12345; Y[0][1] = 54321;
        
        std::cout<<"The size of Y: "<<Y.size()<<endl;
        std::cout<<"The size of X: "<<X.size()<<endl;

        //put some common elements
        for(int i=0; i<10; i++){
            X[i] = Y[i];
        }
        std::vector<std::array<uint64_t, 2>> expected_intersection = get_intersection(X, Y);
        std::cout<<"The expected intersection: "<<expected_intersection.size()<<std::endl;
        for(auto element: expected_intersection){
            std::cout<<element[0]<<element[1]<<", ";
        }
        std::cout<<std::endl;

        std::chrono::high_resolution_clock::time_point time_start, time_stop;
        time_start = std::chrono::high_resolution_clock::now();
        // std::vector<std::array<uint64_t, 2>> get_intersection = offline_fast_poly_psi(Y, X);
        time_stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_dur = time_stop-time_start;
        cout<<"offline time(s): "<<time_dur.count()<<endl;

        // std::cout<<"Intersection get from HE: "<<get_intersection.size()<<std::endl;
        // for(auto element: get_intersection){
        //     std::cout<<element[0]<<element[1]<<", ";
        // }
        // std::cout<<std::endl;
    }

    return 0;
}