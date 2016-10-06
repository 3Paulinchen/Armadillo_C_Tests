// Solving -u''=f(x) with finite differences on (0,pi) in 1D
// with the armadillo library
// complie in console with g++ or gcc
// g++ FD_1D.cpp -o FD_1D -O2 -larmadillo && ./FD_1D

// using armadillo library
#include <armadillo>
// for defining pi constant
#define _USE_MATH_DEFINES 
#include "math.h"

// namespace armadillo
using namespace arma ;

int main(int argc, char** argv)
{
  
  // declaring vectors,floats and int
  vec b,f,u;
  float h;
  int n,i;
  // grid points n
  n = 100;
  // sparse matrix a
  sp_mat A = speye<sp_mat>(n,n);
  // x for plot
  vec x = linspace<vec>(0, M_PI, n);
   
  // setting up laplace stencil in 1D
  // A= 1/h^2 [-1 2 -1]
  
  for (i=1;i<n-1;i++){
	  A(i,i+1)=-1; 
	  }
  for (i=1;i<n;i++){
	  A(i,i)=2; 
	  }
  for (i=1;i<n-1;i++){
	  A(i,i-1)=-1; 
	  } 

  h = x(2)-x(1) ;
  A = A/pow(h,2);
  // function f
  f = sin(x) ;
  // sparse solve u=A\f
  spsolve(u, A, f, "lapack");
  // saving u and x to txt
   u.save("u1.txt", raw_ascii);
   x.save("x1.txt",raw_ascii);
  return 0 ;
}
