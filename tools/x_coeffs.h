#pragma once
#include <cstdint>
#include <vector>
// #include <NTL/ZZ.h>
#include "item.h"

//we can directly use PowerMod in NTL
/* Iterative Function to calculate (x^y)%p in O(log y) */
uint64_t powerMod(uint64_t x, uint64_t y, uint64_t p)  
{  
    uint64_t res = 1;     // Initialize result  
  
    x = x % p; // Update x if it is more than or  
                // equal to p 
   
    if (x == 0) return 0; // In case x is divisible by p; 
  
    while (y > 0)  
    {  
        // If y is odd, multiply x with result  
        if (y & 1)  
            res = (res*x) % p;  
  
        // y must be even now  
        y = y>>1; // y = y/2  
        x = (x*x) % p;  
    }  
    return res;
} 

void polynomial_from_roots(std::vector<uint64_t> &roots, std::vector<uint64_t> &coeffs, uint64_t modulus) {
    coeffs.clear();
    coeffs.resize(roots.size() + 1);
    coeffs[0] = 1;

    for (size_t i = 0; i < roots.size(); i++) {
        // multiply coeffs by (x - root)
        uint64_t neg_root = modulus - (roots[i] % modulus);

        for (size_t j = i + 1; j > 0; j--) {
            coeffs[j] = (coeffs[j - 1] + (neg_root * coeffs[j]) % modulus) % modulus;
        }
        coeffs[0] = (coeffs[0] * neg_root) % modulus;
    }

    return;
}







uint64_t modexp(uint64_t base, uint64_t exponent, uint64_t modulus) {
    uint64_t result = 1;
    while (exponent > 0) {
        if (exponent & 1) {
            result = (result*base)%modulus;
        }
        base = (base*base)%modulus;
        exponent = (exponent >> 1);
    }
    return result;
}

uint64_t modinv(uint64_t x, uint64_t modulus) {
    return modexp(x, modulus - 2, modulus);
}


void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      // if (a < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      // if (b < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
long invMod(long a, long n)
{
  if(a == 0){
    return 0;
  }
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) {
      // InvModError("InvMod: inverse undefined");
      std::cout<<"the gcd is not 1 or 0"<<std::endl;
      std::cout<<a<<", "<<n<<std::endl;
      return -1;
   }
   if (s < 0)
      return s + n;
   else
      return s;
}

//use modinv is too expensive, we suggest not
void dvand (std::vector<uint64_t> &x, const std::vector<uint64_t> &alpha, const std::vector<uint64_t> &b, uint64_t modulus)

//****************************************************************************80
//
//  Purpose:
//
//    DVAND solves a Vandermonde system A' * x = b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ake Bjorck, Victor Pereyra,
//    Solution of Vandermonde Systems of Equations,
//    Mathematics of Computation,
//    Volume 24, Number 112, October 1970, pages 893-903.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double DVAND[N], the solution of the linear system.
//
//reference website: https://cenit.github.io/jburkardt/vandermonde/vandermonde.html
{
  int j;
  int k;
  // int n=0;
  int n = alpha.size();

  x = b;

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = (((modulus + x[j] - x[j-1] )%modulus) * invMod( modulus + alpha[j] - alpha[j-k-1], modulus )) % modulus;
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = (modulus + x[j] - (alpha[k] * x[j+1])%modulus) % modulus;
    }
  }

}



void polynomial_from_points(std::vector<uint64_t> &xs,
                            std::vector<uint64_t> &ys,
                            std::vector<uint64_t> &coeffs,
                            const uint64_t modulus)
{
    // assert(xs.size() == ys.size());
    coeffs.clear();
    coeffs.resize(xs.size());

    if (xs.size() == 0) {
        return;
    }

    // at iteration i of the loop, basis contains the coefficients of the basis
    // polynomial (x - xs[0]) * (x - xs[1]) * ... * (x - xs[i - 1])
    std::vector<uint64_t> basis(xs.size());
    basis[0] = 1;

    // at iteration i of the loop, ddif[j] contains the divided difference
    // [ys[j], ys[j + 1], ..., ys[j + i]]. thus initially, when i = 0,
    // ddif[j] = [ys[j]] = ys[j]

    std::vector<uint64_t> ddif = ys;
    // for(int i=0; i<xs.size()+1; i++){
    //   std::cout<<ddif[i]<<", ";
    // }
    // std::cout<<std::endl<<xs.size()<<std::endl;
    // for(int i=0; i<xs.size()+1; i++){
    //   std::cout<<xs[i]<<", ";
    // }
    // std::cout<<std::endl;

    for (size_t i = 0; i < xs.size(); i++) {
        for (size_t j = 0; j < i + 1; j++) {
            coeffs[j] = (coeffs[j] + (ddif[0]*basis[j])%modulus) % modulus;
        }

        if (i < xs.size() - 1) {
            // update basis: multiply it by (x - xs[i])
            // uint64_t neg_x = modulus - (xs[i] % modulus);
            uint64_t neg_x = modulus - xs[i];

            for (size_t j = i + 1; j > 0; j--) {
                basis[j] = (basis[j - 1] + (neg_x*basis[j])%modulus) % modulus;
            }
            basis[0] = (basis[0]*neg_x)%modulus;

            // update ddif: compute length-(i + 1) divided differences
            for (size_t j = 0; j + i + 1 < xs.size() + 1; j++) {
                // dd_{j,j+i+1} = (dd_{j+1, j+i+1} - dd_{j, j+i}) / (x_{j+i+1} - x_j)
                uint64_t num = (ddif[j + 1] - ddif[j] + modulus) % modulus;
                // std::cout<<(j+1)<<".."<<j<<".."<<ddif[j + 1]<<".."<<ddif[j]<<".."<<num<<"> ";

                uint64_t den = (xs[j + i + 1] - xs[j] + modulus) % modulus;
                // std::cout<<(j+i+1)<<".."<<j<<".."<<xs[j + i + 1]<<".."<<xs[j]<<".."<<den<<"> ";
                // ddif[j] = (num*modinv(den, modulus))%modulus;
                ddif[j] = (num*invMod(den, modulus))%modulus;
                // std::cout<<num<<"."<<den<<"."<<invMod(den, modulus)<<"."<<ddif[j]<<", ";
                
            }
            // std::cout<<std::endl;
        }
        // std::cout<<"basis: "<<std::endl;
        // for(int j=0; j<xs.size(); j++){
        //   std::cout<<basis[j]<<",";
        // }
        // std::cout<<std::endl;
        // std::cout<<"ddif: "<<std::endl;
        // for(int j=0; j<xs.size(); j++){
        //   std::cout<<ddif[j]<<",";
        // }
        // std::cout<<std::endl;
    }
}

void polynomial_from_points(std::vector<std::vector<std::vector<uint32_t>>>  &P_l_coeffs,
                            std::vector<std::vector<std::vector<uint32_t>>>  &P_m_coeffs,
                            std::vector<std::vector<std::vector<item>>> &x_simple_hashing,
                            const uint32_t &b,
                            const uint32_t &modulus)
{
    int omega = x_simple_hashing.size();
    int nu = x_simple_hashing[0].size();
    int B = x_simple_hashing[0][0].size();
    int alpha = B/b;

    int start_pos = 0;
    int cur_pos = 0;
    uint64_t inv = 0;
    uint64_t tmp_r = 0;
    uint32_t tmp_rl = 0, tmp_rm = 0;
    for(int ii=0; ii<omega; ii++){
      // if((ii+1)%2048 == 0){
      //   std::cout<<(ii+1)<<std::endl;
      // }
      for(int jj = 0; jj<nu; jj++){
        for(int kk=0; kk<B; kk++){
          P_l_coeffs[ii][jj][kk] = x_simple_hashing[ii][jj][kk].l;
          P_m_coeffs[ii][jj][kk] = x_simple_hashing[ii][jj][kk].m;
        }
        start_pos = 0;
        for(int kk = 0; kk<alpha; kk++){
          int j, k;
          for ( k = 0; k < b - 1; k++ )
          {
            for ( j = b - 1; k < j; j-- )
            {
              cur_pos = start_pos+j;
              inv = invMod( (modulus + x_simple_hashing[ii][jj][cur_pos].r - x_simple_hashing[ii][jj][cur_pos-k-1].r) % modulus, modulus );
              P_l_coeffs[ii][jj][cur_pos] = (inv * (modulus + P_l_coeffs[ii][jj][cur_pos] - P_l_coeffs[ii][jj][cur_pos-1])) % modulus;
              P_m_coeffs[ii][jj][cur_pos] = (inv * (modulus + P_m_coeffs[ii][jj][cur_pos] - P_m_coeffs[ii][jj][cur_pos-1])) % modulus;
            }
          }

          for ( k = b - 2; 0 <= k; k-- )
          {
            for ( j = k; j < b - 1; j++ )
            {
              tmp_r = x_simple_hashing[ii][jj][start_pos+k].r;
              cur_pos = start_pos+j;

              // tmp_rl = (tmp_r*P_l_coeffs[ii][jj][cur_pos+1]) % modulus;
              // if(P_l_coeffs[ii][jj][cur_pos] < tmp_rl){
              //   P_l_coeffs[ii][jj][cur_pos] = modulus + P_l_coeffs[ii][jj][cur_pos] - tmp_rl;
              // }else{
              //   P_l_coeffs[ii][jj][cur_pos] = P_l_coeffs[ii][jj][cur_pos] - tmp_rl;
              // }
              // tmp_rm = (tmp_r*P_m_coeffs[ii][jj][cur_pos+1]) % modulus;
              // if(P_m_coeffs[ii][jj][cur_pos] < tmp_rm){
              //   P_m_coeffs[ii][jj][cur_pos] = modulus + P_m_coeffs[ii][jj][cur_pos] - tmp_rm;
              // }else{
              //   P_m_coeffs[ii][jj][cur_pos] = P_m_coeffs[ii][jj][cur_pos] - tmp_rm;
              // }
              P_l_coeffs[ii][jj][cur_pos] = (modulus + P_l_coeffs[ii][jj][cur_pos] - (tmp_r*P_l_coeffs[ii][jj][cur_pos+1]) % modulus) % modulus;
              P_m_coeffs[ii][jj][cur_pos] = (modulus + P_m_coeffs[ii][jj][cur_pos] - (tmp_r*P_m_coeffs[ii][jj][cur_pos+1]) % modulus) % modulus;
            }
          }
          start_pos += b;
        }
      }
    }
}

void polynomial_from_points(std::vector<std::vector<std::vector<uint32_t>>>  &P_l_coeffs,
                            std::vector<std::vector<std::vector<item>>> &x_simple_hashing,
                            const uint32_t &b,
                            const uint32_t &modulus)
{
    int omega = x_simple_hashing.size();
    int nu = x_simple_hashing[0].size();
    int B = x_simple_hashing[0][0].size();
    int alpha = B/b;

    
    for(int ii=0; ii<omega; ii++){
      for(int jj = 0; jj<nu; jj++){
        for(int kk=0; kk<B; kk++){
          P_l_coeffs[ii][jj][kk] = x_simple_hashing[ii][jj][kk].l;
        }
        for(int kk = 0; kk<alpha; kk++){
          int j, k;
          for ( k = 0; k < b - 1; k++ )
          {
            for ( j = b - 1; k < j; j-- )
            {
              uint64_t inv = invMod( (modulus + x_simple_hashing[ii][jj][kk*b+j].r - x_simple_hashing[ii][jj][kk*b+j-k-1].r) % modulus, modulus );
              P_l_coeffs[ii][jj][kk*b+j] = (inv * (modulus + P_l_coeffs[ii][jj][kk*b+j] - P_l_coeffs[ii][jj][kk*b+j-1])) % modulus;
            }
          }

          for ( k = b - 2; 0 <= k; k-- )
          {
            for ( j = k; j < b - 1; j++ )
            {
              uint64_t tmp_r = x_simple_hashing[ii][jj][kk*b+k].r;
              // P_l_coeffs[ii][jj][kk*b+j] = (modulus + P_l_coeffs[ii][jj][kk*b+j] - (tmp_r*P_l_coeffs[ii][jj][kk*b+j+1]) % modulus) % modulus;
              // P_m_coeffs[ii][jj][kk*b+j] = (modulus + P_m_coeffs[ii][jj][kk*b+j] - (tmp_r*P_m_coeffs[ii][jj][kk*b+j+1]) % modulus) % modulus;
              // if()
              P_l_coeffs[ii][jj][kk*b+j] = (modulus + P_l_coeffs[ii][jj][kk*b+j] - (tmp_r*P_l_coeffs[ii][jj][kk*b+j+1]) % modulus) % modulus;
            }
          }
        }
      }
    }
}


// void polynomial_from_points(std::vector<std::vector<uint64_t>> &coeffs, std::vector<std::vector<uint64_t>> &xs, std::vector<std::vector<uint64_t>> &ys,
//                             const uint64_t sub_sub_bin_size,
//                             const uint64_t modulus)
// {
//     assert((xs.size() == ys.size()) && (xs[0].size() == ys[0].size()));

//     if (xs.size() == 0) {
//         return;
//     }

//     int w = xs.size();
//     int B = xs[0].size();
//     int b = sub_sub_bin_size;
//     int alpha = B/b;

//     coeffs.resize(w);
//     for(int i=0; i<w; i++){
//       coeffs[i].resize(B);
//     }

//     // at iteration i of the loop, basis contains the coefficients of the basis
//     // polynomial (x - xs[0]) * (x - xs[1]) * ... * (x - xs[i - 1])

//     for(int ii=0; ii< w; ii++){
//       for(int jj=0; jj<alpha; jj++){
//             std::vector<uint64_t> basis(b);
//             basis[0] = 1;

//             // at iteration i of the loop, ddif[j] contains the divided difference
//             // [ys[j], ys[j + 1], ..., ys[j + i]]. thus initially, when i = 0,
//             // ddif[j] = [ys[j]] = ys[j]
//             std::vector<uint64_t> ddif(b);
//             for(int i=0; i<b; i++){
//               ddif[i] = ys[ii][jj*b+i];
//             }

//             for (int i = 0; i < b; i++) {
//                 for (int j = 0; j < i + 1; j++) {
//                     coeffs[ii][jj*b+j] = (coeffs[ii][jj*b+j] + (ddif[0]*basis[j])%modulus) % modulus;
//                 }

//                 if (i < b - 1) {
//                     // update basis: multiply it by (x - xs[i])
//                     uint64_t neg_x = modulus - (xs[ii][jj*b+i] % modulus);

//                     for (int j = i + 1; j > 0; j--) {
//                         basis[j] = (basis[j - 1] + (neg_x*basis[j])%modulus) % modulus;
//                     }
//                     basis[0] = (basis[0]*neg_x)%modulus;

//                     // update ddif: compute length-(i + 1) divided differences
//                     for (int j = 0; j + i + 1 < b + 1; j++) {
//                         // dd_{j,j+i+1} = (dd_{j+1, j+i+1} - dd_{j, j+i}) / (x_{j+i+1} - x_j)
//                         uint64_t num = (ddif[j + 1] - ddif[j] + modulus) % modulus;
//                         uint64_t den = (xs[ii][jj*b+j+i+1] - xs[ii][jj*b+j] + modulus) % modulus;
//                         // ddif[j] = (num*modinv(den, modulus))%modulus;
//                         ddif[j] = (num*invMod(den, modulus))%modulus;
//                     }
//                 }
//             }
//       }
//     }

// }


//use lagrane inpterpolation to get the polunomial coeffs. However, supposing n is the vector size, it requests one more element vec[n]=0. Therefore,
//the following function is not right.
// void polynomial_from_points(std::vector<std::vector<std::vector<uint64_t>>>  &P_l_coeffs,
//                             std::vector<std::vector<std::vector<uint64_t>>>  &P_m_coeffs,
//                             std::vector<std::vector<std::vector<item>>> &x_simple_hashing,
//                             const uint32_t &b,
//                             const uint64_t &modulus)
// {


//     int omega = x_simple_hashing.size();
//     int nu = x_simple_hashing[0].size();
//     int B = x_simple_hashing[0][0].size();
//     int alpha = B/b;

//     std::vector<uint64_t> basis(b);//basis is only related with the root, which is the same for l and m
//     std::vector<uint64_t> ddif_l(b), ddif_m(b);
//     for(int ii=0; ii<omega; ii++){
//       for(int jj = 0; jj<nu; jj++){
//         for(int kk = 0; kk<alpha; kk++){
//               // at iteration i of the loop, basis contains the coefficients of the basis
//             // polynomial (x - xs[0]) * (x - xs[1]) * ... * (x - xs[i - 1])
//             // std::vector<uint64_t> basis(b);
//             basis[0] = 1;

//             // at iteration i of the loop, ddif[j] contains the divided difference
//             // [ys[j], ys[j + 1], ..., ys[j + i]]. thus initially, when i = 0,
//             // ddif[j] = [ys[j]] = ys[j]
//             // std::vector<uint64_t> ddif(b);
//             for(int t=0; t<b; t++){
//               ddif_l[t] = x_simple_hashing[ii][jj][kk*b+t].l;
//               ddif_m[t] = x_simple_hashing[ii][jj][kk*b+t].m;
//             }
//             for(int t=0; t<b; t++){
//               std::cout<<ddif_l[t]<<", ";
//             }
//             std::cout<<std::endl<<b<<std::endl;
//             for(int t=0; t<b; t++){
//               std::cout<<x_simple_hashing[ii][jj][kk*b+t].r<<", ";
//             }
//             std::cout<<std::endl;

//             for (int i = 0; i < b; i++) {
//                 for (int j = 0; j < i + 1; j++) {
//                     P_l_coeffs[ii][jj][kk*b+j] = (P_l_coeffs[ii][jj][kk*b+j] + (ddif_l[0]*basis[j])%modulus) % modulus;
//                     P_m_coeffs[ii][jj][kk*b+j] = (P_m_coeffs[ii][jj][kk*b+j] + (ddif_m[0]*basis[j])%modulus) % modulus;
//                 }

//                 if (i < b - 1) {
//                     // update basis: multiply it by (x - xs[i])
//                     uint64_t neg_x = modulus - x_simple_hashing[ii][jj][kk*b+i].r;

//                     for (int j = i + 1; j > 0; j--) {
//                         basis[j] = (basis[j - 1] + (neg_x*basis[j])%modulus) % modulus;
//                     }
//                     basis[0] = (basis[0]*neg_x)%modulus;

//                     // update ddif: compute length-(i + 1) divided differences
//                     for (int j = 0; j + i + 1 < b + 1; j++) {
//                         // dd_{j,j+i+1} = (dd_{j+1, j+i+1} - dd_{j, j+i}) / (x_{j+i+1} - x_j)
//                         uint64_t num_l = (modulus + ddif_l[j + 1] - ddif_l[j]) % modulus;
//                         uint64_t den = (modulus + x_simple_hashing[ii][jj][kk*b+j+i+1].r - x_simple_hashing[ii][jj][kk*b+j].r) % modulus;
//                         std::cout<<(j+i+1)<<".."<<j<<".."<<x_simple_hashing[ii][jj][kk*b+j+i+1].r<<".."<<x_simple_hashing[ii][jj][kk*b+j].r<<".."<<den<<"> ";
//                         // ddif[j] = (num*modinv(den, modulus))%modulus;
//                         auto inv = invMod(den, modulus);
//                         ddif_l[j] = (num_l*inv)%modulus;
//                         std::cout<<num_l<<"."<<den<<"."<<inv<<"."<<ddif_l[j]<<", ";

//                         uint64_t num_m = (modulus + ddif_m[j + 1] - ddif_m[j]) % modulus;
//                         ddif_m[j] = (num_m*inv)%modulus;
//                     }
//                     std::cout<<std::endl;
//                 }
//                 std::cout<<"basis: "<<std::endl;
//                 for(int j=0; j<b; j++){
//                   std::cout<<basis[j]<<",";
//                 }
//                 std::cout<<std::endl;
//                 std::cout<<"ddif_l: "<<std::endl;
//                 for(int j=0; j<b; j++){
//                   std::cout<<ddif_l[j]<<",";
//                 }
//                 std::cout<<std::endl;
    
//             }
//         }
//       }
//     }
// }
