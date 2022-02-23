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
  if(a == 0 || a == n){
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

// void dvandprg ( int n, double alpha[], double b[], double x[], double c[], 
//   double m[] )

// //****************************************************************************80
// //https://cenit.github.io/jburkardt/vandermonde/vandermonde.cpp
// //  Purpose:
// //
// //    DVANDPRG solves a Vandermonde system A' * x = f progressively.
// //
// //  Discussion:
// //
// //    This function receives the solution to the system of equations A' * x = f
// //    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
// //    and new values alpha(n) and f(n).  It updates the solution.
// //
// //    To solve a system of Nbig equations, this function may be called 
// //    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
// //    current subsystem is returned.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    18 April 2014
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Reference:
// //
// //    Ake Bjorck, Victor Pereyra,
// //    Solution of Vandermonde Systems of Equations,
// //    Mathematics of Computation,
// //    Volume 24, Number 112, October 1970, pages 893-903.
// //
// //  Parameters:
// //
// //    Input, int N, the new order of the matrix, which is 1 
// //    larger than on the previous call.  For the first call, N must be 1.
// //
// //    Input, double ALPHA[N], the parameters that define the matrix.
// //    The values should be distinct.  The value ALPHA(N) has just been
// //    added to the system.
// //
// //    Input, double B[N], the right hand side of the linear system.
// //
// //    Input/output, double X[N].  On input, the first N-1 entries 
// //    contain the solution of the N-1xN-1 linear system.  On output, the 
// //    solution to the NxN linear system.
// //
// //    Input/output, double C[N], M[N].  On input, the first N-1 
// //    entries contain factorization data for the N-1xN-1 linear system.  On 
// //    output, factorization data for the NxN linear system.
// //
// {
//   double cn;
//   int j;
 
//   c[n-1] = b[n-1];
//   for ( j = n - 1; 1 <= j; j-- )
//   {
//     c[j-1] = ( c[j] - c[j-1] ) / ( alpha[n-1] - alpha[j-1] );
//   }

//   if ( n == 1 )
//   {
//     m[n-1] = 1.0;
//   }
//   else
//   {
//     m[n-1] = 0.0;
//   }

//   cn = c[0];
//   x[n-1] = c[0];

//   for ( j = n - 1; 1 <= j; j-- )
//   {
//     m[j] = m[j] - alpha[n-2] * m[j-1];
//     x[n-j-1] = x[n-j-1] + m[j] * cn;
//   }

//   return;
// }


void dvandprg ( int n, const std::vector<uint64_t> &alpha, const std::vector<uint64_t> &b, std::vector<uint64_t> &x, std::vector<uint64_t> &c, 
  std::vector<uint64_t> &m, const uint64_t &modulus)

//****************************************************************************80
//https://cenit.github.io/jburkardt/vandermonde/vandermonde.cpp
//  Purpose:
//
//    DVANDPRG solves a Vandermonde system A' * x = f progressively.
//
//  Discussion:
//
//    This function receives the solution to the system of equations A' * x = f
//    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
//    and new values alpha(n) and f(n).  It updates the solution.
//
//    To solve a system of Nbig equations, this function may be called 
//    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
//    current subsystem is returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2014
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
//    Input, int N, the new order of the matrix, which is 1 
//    larger than on the previous call.  For the first call, N must be 1.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.  The value ALPHA(N) has just been
//    added to the system.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N].  On input, the first N-1 entries 
//    contain the solution of the N-1xN-1 linear system.  On output, the 
//    solution to the NxN linear system.
//
//    Input/output, double C[N], M[N].  On input, the first N-1 
//    entries contain factorization data for the N-1xN-1 linear system.  On 
//    output, factorization data for the NxN linear system.
//
{
  uint64_t cn;
  int j;
 
  c[n-1] = b[n-1];
  for ( j = n - 1; 1 <= j; j-- )
  {
    // c[j-1] = ( c[j] - c[j-1] ) / ( alpha[n-1] - alpha[j-1] );
    c[j-1] = (((modulus + c[j] - c[j-1] )%modulus) * invMod( modulus + alpha[n-1] - alpha[j-1], modulus )) % modulus;
  }

  if ( n == 1 )
  {
    m[n-1] = 1;
  }
  else
  {
    m[n-1] = 0;
  }

  cn = c[0];
  x[n-1] = c[0];

  for ( j = n - 1; 1 <= j; j-- )
  {
    m[j] = (modulus + m[j] - (alpha[n-2] * m[j-1]) % modulus) % modulus;
    x[n-j-1] = (x[n-j-1] + m[j] * cn ) % modulus;
  }

  return;
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

