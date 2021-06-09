// #include "poly_coeffs.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZX.h>
#include <iostream>
#include <sys/time.h>
#include <vector>
#include "x_coeffs.h"
#include "polyFFT.h"
#include "item.h"
#include <random>

using namespace std;
// Timing routines
static double getMillies(timeval start, timeval end)
{
	long time1 = (start.tv_sec * 1000000) + (start.tv_usec );
	long time2 = (end.tv_sec * 1000000) + (end.tv_usec );

	return (double)(time2-time1)/1000;
}

void testLabel(){
    timeval begin, end;
    double t = 0.0;

    uint64_t m = 200;
    // uint64_t plain_modulus = 8519681;
    uint64_t plain_modulus = 65537;

    // uint64_t plain_modulus = 274877906951;

    NTL::Vec<NTL::ZZ> roots, labels, result;

    roots.SetLength(m);
    labels.SetLength(m);

    for(int i=0; i<m; i++){
        roots[i] = NTL::RandomBits_ZZ(16) % plain_modulus;
        labels[i] = NTL::RandomBits_ZZ(16) % plain_modulus;
        // cout<<alpha[i]<<", ";
        // mul_res *= alpha[i];
        // add_res += alpha[i];
    }
    // for(int i=0; i<m; i++){
    //     std::cout<<roots[i]<<", ";
    // }
    // std::cout<<std::endl;
    // for(int i=0; i<m; i++){
    //     std::cout<<labels[i]<<", ";
    // }
    // std::cout<<std::endl;

//test which function with the labels are faster
    int count = 1<<4;

    NTL::ZZ_pX poly;
    NTL::ZZ_p::init(NTL::ZZ(plain_modulus));
    NTL::Vec<NTL::ZZ_p> my_roots;
    NTL::Vec<NTL::ZZ_p> my_labels;
    my_roots.SetLength(m);
    my_labels.SetLength(m);
    for(int i=0; i<m; i++){
        my_roots[i] = to_ZZ_p(roots[i]);
        my_labels[i] = to_ZZ_p(labels[i]);
    }
    uint64_t temp_coeffs = 0;
    gettimeofday(&begin, NULL);
    for(int i=0; i<count; i++){
        interpolate(poly, my_roots, my_labels);
        for(int i=0; i<m; i++){
            BytesFromZZ((unsigned char*)&temp_coeffs, rep(coeff(poly, i)), sizeof(uint32_t));
            // temp_coeffs = rep(coeff(poly, i));
            // std::cout<<temp_coeffs<<", "<<coeff(poly, i)<<":  ";
        }
    }

    // std::cout<<std::endl;
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs by NTL::interplot: "<<t<<std::endl;
    std::cout<<"the degree of this poly: "<<deg(poly)<<std::endl;

    // long degree = deg(poly);
    // for(int i=0; i<=degree; i++){
    //     std::cout<<coeff(poly, i)<<", ";
    // }
    // std::cout<<endl;

    NTL::ZZ plain_modulus_zz;
    plain_modulus_zz = plain_modulus;

    std::vector<uint64_t> m_coeffs(m);
    
    gettimeofday(&begin, NULL);
    for(int i=0; i<count; i++){
        interpolate_zp(poly, my_roots.data(), my_labels.data(), (long)m-1, (int)1, plain_modulus_zz);

        // for(int i=0; i<m; i++){
        //     BytesFromZZ((unsigned char*)&temp_coeffs, rep(coeff(poly, i)), sizeof(uint32_t));
        //     // temp_coeffs = rep(coeff(poly, i));
        //     // std::cout<<temp_coeffs<<", "<<coeff(poly, i)<<":  ";
        // }
    }

    // std::cout<<std::endl;
    
    // NTL::ZZ m_zz = coeff(poly);
    
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs by spot interpolate: "<<t<<std::endl;
    std::cout<<"the degree of this poly: "<<deg(poly)<<std::endl;

    long degree = deg(poly);
    // for(int i=0; i<=degree; i++){
    //     std::cout<<coeff(poly, i)<<", ";
    // }
    // std::cout<<endl;

    NTL::Vec<NTL::ZZ_p> inter_coeffs;
    inter_coeffs.SetLength(m);
    for(int i=0; i<m; i++){
        inter_coeffs[i] = coeff(poly, i);
    }

    std::vector<uint64_t> xs(m), ys(m), coeffs(m);
    for(int i=0; i<m; i++){
        xs[i] = roots[i]%plain_modulus;
        ys[i] = labels[i]%plain_modulus;
    }
    gettimeofday(&begin, NULL);
    for(int i=0; i<count; i++){
        polynomial_from_points(xs, ys, coeffs, plain_modulus);
    }
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs by LU: "<<t<<std::endl;

    // // std::vector<uint64_t> xs(m), ys(m), coeffs(m);
    // gettimeofday(&begin, NULL);
    // dvand(coeffs, xs, ys, plain_modulus);
    // gettimeofday(&end, NULL);
    // t = getMillies(begin, end);
    // std::cout<<"time for computing the coeffs by new dvand: "<<t<<std::endl;
            // std::vector<uint64_t> xs(m), ys(m), coeffs(m);
    std::vector<uint64_t> my_xs(m), my_ys(m), my_coeffs(m);
    for(int i=0; i<m; i++){
        my_xs[i] = roots[i]%plain_modulus;
        my_ys[i] = labels[i]%plain_modulus;
    }
    gettimeofday(&begin, NULL);
    dvand(my_coeffs, my_xs, my_ys, plain_modulus);
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs by new dvand: "<<t<<std::endl;
    // for(int i=0; i<m; i++){
    //     std::cout<<my_coeffs[i]<<", ";
    // }
    // std::cout<<std::endl;


    std::cout<<"Now check the results: "<<std::endl;
    for(int i=0; i<coeffs.size(); i++){
        if(my_coeffs[i] != coeffs[i]){
            std::cout<<"\t something is wrong!"<<std::endl;
            // std::cout<<"\t"<<i<<": "<<my_result[i]<<", "<<coeffs[i]<<endl;
            break;
        }
    }
    // for(int i=0; i<m; i++){
    //     std::cout<<my_result[i]<<", ";
    // }
    // std::cout<<std::endl;
    // for(int i=0; i<m; i++){
    //     std::cout<<coeffs[i]<<", ";
    // }
    std::cout<<std::endl;
    // for(int i=0; i<m; i++){
    //     std::cout<<my_coeffs[i]<<", ";
    // }
    // std::cout<<std::endl;



    std::cout<<"The expected y[0]: "<<labels[0]<<", "<<ys[0]<<std::endl;
    NTL::ZZ ans;
    ans = my_coeffs[0];
    for(int i=1; i<m; i++){
        ans += my_coeffs[i]*NTL::power(roots[0], i);
    }
    // ans += NTL::power(roots[0], m);
    std::cout<<"my_result for dvand: "<<ans%plain_modulus<<std::endl;
    
    ans = coeffs[0];
    for(int i=1; i<m; i++){
        ans += coeffs[i]*NTL::power(roots[0], i);
    }
    // ans += NTL::power(roots[0], m);
    std::cout<<"coeffs for poly: "<<ans%plain_modulus<<std::endl;

    ans = my_coeffs[0];
    for(int i=1; i<m; i++){
        ans += my_coeffs[i]*NTL::power(roots[0], i);
    }
    // ans += NTL::power(roots[0], m);
    std::cout<<"my_coeffs for poly: "<<ans%plain_modulus<<std::endl;

    NTL::ZZ_p int_ans=inter_coeffs[0];
    // int_ans = to_ZZ_p(ZZ(0));
    for(int i=1; i<m; i++){
        int_ans += inter_coeffs[i]*power(my_roots[0], i);
    }
    std::cout<<"the ans from NTL::interplot: "<<int_ans<<std::endl;

    gettimeofday(&begin, NULL);
    polynomial_from_roots(xs, coeffs, plain_modulus);
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs with no labels: "<<t<<std::endl;

    gettimeofday(&begin, NULL);
    BuildFromRoots(poly, my_roots);
    gettimeofday(&end, NULL);
    t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs with no labels: "<<t<<std::endl;
}

void testNoLabel(){
        // double alpha[5]={0, 1, -1, 2, -2};
    // double f[5] = {-5, -3, -15, 39, -9};
    //c[5] should be {-5, 4, -7, 2, 3}
    timeval begin, end;

    // NTL::Vec<NTL::ZZ> alpha;//={1091432, 1470396, 1103006, 347220, 952904};
    // std::vector<uint64_t> alpha_int={1091432, 1470396};//, 1103006};//, 347220, 952904};
    // alpha.SetLength(alpha_int.size());
    // for(int i=0; i<alpha_int.size(); i++){
    //     alpha[i] = alpha_int[i];
    // }

    // int n = alpha_int.size();
    // alpha.SetLength(n);

    int n = 2000;
    // int n = 2000;
    NTL::Vec<NTL::ZZ> alpha;
    alpha.SetLength(n);
    std::vector<uint64_t> alpha_int(n);
    // alpha_int.resize(n);
    // NTL::ZZ mul_res(1);
    // NTL::ZZ add_res(0);
    for(int i=0; i<n; i++){
        long rnd = NTL::RandomBnd(1<<30);
        alpha[i] = rnd;
        alpha_int[i] = rnd; 
        // cout<<alpha[i]<<", ";
        // mul_res *= alpha[i];
        // add_res += alpha[i];
    }
    // cout<<endl;
    // cout<<"the mul: "<<mul_res<<endl;
    // cout<<"add res: "<<add_res<<endl;
    
    NTL::Vec<NTL::ZZ> f(alpha);
    // f.SetLength(n);
    // for(int i=0; i<3; i++){
    //     for(int j=0; j<f.length(); j++){
    //         f[j] = 
    //     }
    // }
    for(int j=0; j<f.length(); j++){
        power(f[j], f[j], n);
    }
    //if n is even, f should be negated; otherwise do not need
    if(n%2 == 0){
        NTL::negate(f, f);
    }

    NTL::Vec<NTL::ZZ> result;
    result.SetLength(n);

    uint64_t plain_modulus = 8519681;


    std::vector<uint64_t> result_int;
    gettimeofday(&begin, NULL);
    polynomial_from_roots(alpha_int, result_int, plain_modulus);
    gettimeofday(&end, NULL);

    auto t = getMillies(begin, end);
    std::cout<<"time for computing the coeffs: "<<t<<std::endl;
    // for(int i=0; i<n; i++){
    //     std::cout<<result_int[i]<<", ";
    // }
    // std::cout<<endl;

    // for(int i=0; i<n; i++){
    //     if((result[i]%plain_modulus) != result_int[i]){
    //         cout<<"find unequal!"<<endl;
    //     }
    // }
    // std::cout<<std::endl;
}

void testModInv(){
    uint64_t plain_modulus = 8519681;
    for(int i=0; i<1000; i++){
        long rd = NTL::RandomBnd(1<<30);
        if(NTL::InvMod(rd, plain_modulus) != modinv(rd, plain_modulus)){
            std::cout<<"something is wrong"<<std::endl;
            break;
        }
    }
    std::cout<<std::endl;

    std::vector<uint32_t> rds(100000);
    for(int i=0; i<rds.size(); i++){
        rds[i] = NTL::RandomBnd(1<<30);
    }
    timeval begin, end;
    gettimeofday(&begin, nullptr);
    for(int i=0; i<rds.size(); i++){
        NTL::InvMod(rds[i], plain_modulus);
    }
    gettimeofday(&end, nullptr);
    std::cout<<"time for ntl xgcd: "<<getMillies(begin, end)<<std::endl;

    gettimeofday(&begin, nullptr);
    for(int i=0; i<rds.size(); i++){
        modinv(rds[i], plain_modulus);
    }
    gettimeofday(&end, nullptr);
    std::cout<<"time for ntl exponent: "<<getMillies(begin, end)<<std::endl;
}

void test_poly_from_roots(){
    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<uint64_t> distrib(0, (uint64_t)18446744073709551615); //18446744073709551615=2^64-1
    std::uniform_int_distribution<uint64_t> distrib; //18446744073709551615=2^64-1
    uint64_t rnd = 0;

    uint32_t omega = 2;
    uint32_t nu = 2;
    uint32_t B = 8;
    uint32_t b = B;
    uint64_t modulus = 65537;
    std::vector<std::vector<std::vector<uint32_t>>> P_l_coeffs(omega, std::vector<std::vector<uint32_t>>(nu, std::vector<uint32_t>(B, 0)));
    std::vector<std::vector<std::vector<uint32_t>>> P_m_coeffs(omega, std::vector<std::vector<uint32_t>>(nu, std::vector<uint32_t>(B, 0)));

    std::vector<std::vector<std::vector<item>>> simple_hashing(omega, std::vector<std::vector<item>>(nu, std::vector<item>(B)));
    for(int i=0; i<omega; i++){
        for(int j=0; j<nu; j++){
            for(int k=0; k<b; k++){
                rnd = distrib(gen);
                simple_hashing[i][j][k].r = rnd % modulus; rnd /= modulus;
                simple_hashing[i][j][k].m = rnd % modulus; rnd /= modulus;
                simple_hashing[i][j][k].l = rnd % modulus;
            }
        }
    }
    simple_hashing[0][0][0].r = 65486; simple_hashing[0][0][0].l = 38777;
    simple_hashing[0][0][1].r = 54916; simple_hashing[0][0][1].l = 45578;
    simple_hashing[0][0][2].r = 6111; simple_hashing[0][0][2].l = 13409;
    simple_hashing[0][0][3].r = 61803; simple_hashing[0][0][3].l = 48223;
    simple_hashing[0][0][4].r = 37741; simple_hashing[0][0][4].l = 45541;
    simple_hashing[0][0][5].r = 31480; simple_hashing[0][0][5].l = 54503;
    simple_hashing[0][0][6].r = 14225; simple_hashing[0][0][6].l = 7384;
    simple_hashing[0][0][7].r = 49488; simple_hashing[0][0][7].l = 30636;
    std::cout<<"the polynomial coeffs should be: "<<std::endl;
    std::vector<uint64_t> my_x(b), my_y_l(b), my_y_m(b);
    // for(int k=0; k<b+1; k++){
    //     std::cout<<my_x[k]<<", ";
    // }
    // std::cout<<std::endl;
    // for(int k=0; k<b+1; k++){
    //     std::cout<<my_y_l[k]<<", ";
    // }
    // std::cout<<std::endl;

    std::vector<uint64_t> coeff_l(b), coeff_m(b);
    for(int i=0; i<omega; i++){
        for(int j=0; j<nu; j++){
            for(int k=0; k<b; k++){
                my_x[k] = simple_hashing[i][j][k].r;
                my_y_m[k] = simple_hashing[i][j][k].m;
                my_y_l[k] = simple_hashing[i][j][k].l;
            }
            // for(int k=0; k<b+1; k++){
            //     std::cout<<my_x[k]<<", ";
            // }
            // std::cout<<std::endl;
            // for(int k=0; k<b+1; k++){
            //     std::cout<<my_y_l[k]<<", ";
            // }
            // std::cout<<std::endl;

            polynomial_from_points(my_x, my_y_l, coeff_l, modulus);
            // if(i==0 && j==0){
            //     uint64_t ans = 0;
            //     uint64_t x = my_x[0];
            //     uint64_t x_power = 1;
            //     for(int k=0; k<b; k++){
            //         ans = (ans + (coeff_l[k]*x_power) % modulus) % modulus;
            //         x_power = (x_power*x) % modulus;
            //     }
            //     std::cout << std::endl << ans <<"&" <<my_y_l[0]<<std::endl;
            // }
            for(int k=0; k<b; k++){
                std::cout<<coeff_l[k]<<", ";
            }
            std::cout<<endl;
        }
    }
    std::cout<<"poly_m: "<<std::endl;
    for(int i=0; i<omega; i++){
        for(int j=0; j<nu; j++){
            for(int k=0; k<b; k++){
                my_x[k] = simple_hashing[i][j][k].r;
                my_y_m[k] = simple_hashing[i][j][k].m;
                my_y_l[k] = simple_hashing[i][j][k].l;
            }
            polynomial_from_points(my_x, my_y_m, coeff_m, modulus);

            // for(int k=0; k<b; k++){
            //     uint64_t ans = 0;
            //     uint64_t x = my_x[k];
            //     uint64_t x_power = 1;
            //     for(int t=0; t<b; t++){
            //         ans = (ans + (coeff_m[t]*x_power) % modulus) % modulus;
            //         x_power = (x_power*x) % modulus;
            //     }
            //     std::cout << std::endl << ans <<"&" <<my_y_m[k]<<"\t";
            // }
            // std::cout<<std::endl;

            for(int k=0; k<b; k++){
                std::cout<<coeff_m[k]<<", ";
            }            
            std::cout<<endl;        
        }
    }


    std::cout<<"we get the polynomial: "<<std::endl;
    polynomial_from_points(P_l_coeffs, P_m_coeffs, simple_hashing, b, modulus);

    for(int i=0; i<omega; i++){
        for(int j=0; j<nu; j++){
            for(int k=0; k<b; k++){
                std::cout<<P_l_coeffs[i][j][k]<<", ";
            }
            std::cout<<endl;
        }
    }
    std::cout<<std::endl;
    for(int i=0; i<omega; i++){
        for(int j=0; j<nu; j++){
            for(int k=0; k<b; k++){
                std::cout<<P_m_coeffs[i][j][k]<<", ";
            }
            std::cout<<endl;
        }
    }

    // for(int i=0; i<omega; i++){
    //     for(int j=0; j<nu; j++){
    //         for(int k=0; k<b; k++){
    //             uint64_t ans_l = 0, ans_m = 0;
    //             uint64_t x = simple_hashing[i][j][k].r;
    //             uint64_t x_power = 1;
    //             for(int t=0; t<b; t++){
    //                 ans_l = (ans_l + (P_l_coeffs[i][j][t]*x_power) % modulus) % modulus;
    //                 ans_m = (ans_m + (P_m_coeffs[i][j][t]*x_power) % modulus) % modulus;

    //                 x_power = (x_power*x) % modulus;
    //             }
    //             std::cout << ans_l <<"&" <<simple_hashing[i][j][k].l<<"\t";
    //             std::cout << ans_m <<"&" <<simple_hashing[i][j][k].m<<"\t\t";
    //         }
    //         std::cout<<std::endl;
    //         // if(true){
    //         // // if(true){
    //         //     uint64_t ans_l = 0, ans_m;
    //         //     uint64_t x = simple_hashing[1][0][0].r;
    //         //     uint64_t x_power = 1;
    //         //     for(int k=0; k<b; k++){
    //         //         ans_l = (ans_l + (P_l_coeffs[i][j][k]*x_power) % modulus) % modulus;
    //         //         ans_m = (ans_m + (P_m_coeffs[i][j][k]*x_power) % modulus) % modulus;

    //         //         x_power = (x_power*x) % modulus;
    //         //     }
    //         //     std::cout << std::endl << ans_l <<"&" <<simple_hashing[1][0][0].l<<std::endl;
    //         //     std::cout << std::endl << ans_m <<"&" <<simple_hashing[1][0][0].m<<std::endl;
    //         // }
    //     }
    // }

}

int main(){
    // testLabel();

    // int m=2;
    // std::vector<double> x(m);
    // std::vector<uint64_t> alpha = {8143845, 3461675};
    // std::vector<uint64_t> b = {5903452, 1787124};
    // dvand(x, alpha, b);

    // for(int i=0; i<x.size(); i++){
    //     std::cout<<x[i]<<", ";
    // }
    // std::cout<<std::endl;

    test_poly_from_roots();

    // std::vector<uint64_t> roots = {65486, 733, 54364, 36329, 15850, 15272, 55884, 31550};
    // std::vector<uint64_t> tags = {8142, 47191, 26723, 36524, 15230, 52097, 39589, 27735};
    // std::vector<uint64_t> coeffs(8);
    // polynomial_from_points(roots, tags, coeffs, 65537);
    // for(int i=0; i<8; i++){
    //     std::cout<<coeffs[i]<<", ";
    // }
    // std::cout<<std::endl;



}