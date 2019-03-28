# ifndef RBF_H
# define RBF_H

# include <iostream>
# include <cmath>
# include <typd.h>


/*
 * functors for radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */


// pure virtual functor defining interface for radial basis functions
   struct rbf_f
  {
		rbf_f(){}
		virtual ~rbf_f(){}

      virtual REAL_   rbf( REAL_  r ) = 0;   // value of rbf at r
      virtual REAL_   d(   REAL_  r ) = 0;   // gradient of rbf at r
      virtual REAL_   d2(  REAL_  r ) = 0;   // second gradient of rbf at r

      REAL_  operator()( REAL_  r ){ return rbf( r ); }
  };


// multiquadratric radial basis function: y = ( r^2 + r0^2 )^0.5
   struct rbf_multiquadratic : rbf_f
  {
      REAL_  r0;
      REAL_  r02;

      rbf_multiquadratic( REAL_  scale=1. ) : r0(scale), r02( scale*scale ) {}

      REAL_   rbf( REAL_  r ){ return sqrt( r*r + r02 ); }
      REAL_   d(   REAL_  r ){ return r/rbf(r); }
      REAL_   d2(  REAL_  r ){ return ( 1. - r*d(r)/rbf(r) ) / rbf(r); }
  };


// inverse multiquadratric radial basis function: y(r) = ( r^2 + r0^2 )^-0.5
   struct rbf_invmultiquadratic : rbf_f
  {
      REAL_  r0;
      REAL_  r02;

      rbf_invmultiquadratic( REAL_  scale=1. ) : r0(scale), r02( scale*scale ) {}

      REAL_   rbf( REAL_  r ){ return 1./sqrt( r*r + r02 ); }
      REAL_   d(   REAL_  r ){ return -r*pow(rbf(r),3); }
      REAL_   d2(  REAL_  r ){ return ( 3.*r*d(r)/rbf(r) - 1. ) / ( rbf(r)*rbf(r)*rbf(r) ); }
  };


// thin-plate spline radial basis function: y(r) = r^2 * log( r/r0 )
   struct rbf_thinplate : rbf_f
  {
      REAL_  r0;
      REAL_  r02;

      rbf_thinplate( REAL_  scale=1. ) : r0(scale), r02( scale*scale ) {}

      REAL_  rbf( REAL_  r ){ return r <= 0. ? 0. : r*r*log( r/r0 ); }
      REAL_  d(   REAL_  r ){ return r/r02 + 2.*rbf(r)/r; }
      REAL_  d2(  REAL_  r ){ return ( 1./r02 + 2.*d(r)/r - 2.*rbf(r)/(r*r) ); }
  };


// gaussian radial basis function: y(r) = exp( -0.5*r*r/(r0*r0) )
   struct rbf_gaussian : rbf_f
  {
      REAL_  r0;
      REAL_  r02;

      rbf_gaussian( REAL_  scale=1. ) : r0(scale), r02( scale*scale ) {}

      REAL_   rbf( REAL_  r ){ return exp( -0.5*r*r/r02 ); }
      REAL_   d(   REAL_  r ){ return -rbf(r)*r/r0; }
      REAL_   d2(  REAL_  r ){ return ( rbf(r)*( r*r/r0 -1. )/r0 ); }
  };


// biharmonic radial basis function: y(r) = r
   struct rbf_biharmonic : rbf_f
  {
      REAL_  rbf( REAL_  r ){ return r; }
      REAL_  d(    REAL_  r ){ return 1; }
      REAL_  d2(   REAL_  r ){ return 0; }
  };


// triharmonic radial basis function: y(r) = r^3
   struct rbf_triharmonic : rbf_f
  {
      REAL_  rbf( REAL_  r ){ return r*r*r; }
      REAL_  d(   REAL_  r ){ return 3*r*r; }
      REAL_  d2(  REAL_  r ){ return 6*r;   }
  };


# endif
