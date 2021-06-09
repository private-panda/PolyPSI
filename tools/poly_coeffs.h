#pragma once

#include <string.h>
#include <stdint.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
// #include <NTL/vec_ZZ.h>

double *dvand ( int n, double alpha[], double b[] )

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
  double *x= new double[n];

//   x = r8vec_copy_new ( n, b );
  memcpy(x, b, n*sizeof(double));

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = ( x[j] - x[j-1] ) / ( alpha[j] - alpha[j-k-1] );
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = x[j] - alpha[k] * x[j+1];
    }
  }

  return x;
}

//this is what we want
void dvand (NTL::Vec<NTL::ZZ> &x, const NTL::Vec<NTL::ZZ> &alpha, const NTL::Vec<NTL::ZZ> &b)

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
  int n = x.length();
  // for(int i=0; i<alpha.length(); i++){
  //   if(alpha[i] == 0){ //0 is the dummy element of X_hash
  //     n = i;
  //     break;
  //   }
  // }
  // for(int i=0; i<x.length(); i++){
  //   x[i] = 0; // n can be smaller than alpha.length(), then set the last (alpha.length()-n) coefficients as 0
  // }
//   double *x= new double[n];
//   NTL::ZZ *x = new NTL::ZZ[n];
  
//   x = r8vec_copy_new ( n, b );
//   memcpy(x, b, n*sizeof(NTL::ZZ));
//   NTL::VectorCopy(x, b, n);
  x = b;

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = ( x[j] - x[j-1] ) / ( alpha[j] - alpha[j-k-1] );
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = x[j] - alpha[k] * x[j+1];
    }
  }

  return;
}