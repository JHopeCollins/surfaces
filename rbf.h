# ifndef RBF_H
# define RBF_H

# include <iostream>
# include <math.h>


/*
 * functors for radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */


// pure virtual functor defining interface for radial basis functions
   struct rbf_f
  {
      virtual double  rbf( double r ) = 0;   // value of rbf at r
      virtual double  d(   double r ) = 0;   // gradient of rbf at r
      virtual double  d2(  double r ) = 0;   // second gradient of rbf at r

      double operator()( double r ){ return rbf( r ); }
  };


// multiquadratric radial basis function: y = ( r^2 + r0^2 )^0.5
   struct rbf_multiquadratic : rbf_f
  {
      double r0;
      double r02;

      rbf_multiquadratic( double scale=1. ) : r0(scale), r02( scale*scale ) {}

      double  rbf( double r ){ return sqrt( r*r + r02 ); }
      double  d(   double r ){ return r/rbf(r); }
      double  d2(  double r ){ return ( 1. - r*d(r)/rbf(r) ) / rbf(r); }
  };


// inverse multiquadratric radial basis function: y(r) = ( r^2 + r0^2 )^-0.5
   struct rbf_invmultiquadratic : rbf_f
  {
      double r0;
      double r02;

      rbf_invmultiquadratic( double scale=1. ) : r0(scale), r02( scale*scale ) {}

      double  rbf( double r ){ return 1./sqrt( r*r + r02 ); }
      double  d(   double r ){ return -r*pow(rbf(r),3); }
      double  d2(  double r ){ return ( 3.*r*d(r)/rbf(r) - 1. ) / ( rbf(r)*rbf(r)*rbf(r) ); }
  };


// thin-plate spline radial basis function: y(r) = r^2 * log( r/r0 )
   struct rbf_thinplate : rbf_f
  {
      double r0;
      double r02;

      rbf_thinplate( double scale=1. ) : r0(scale), r02( scale*scale ) {}

      double rbf( double r ){ return r <= 0. ? 0. : r*r*log( r/r0 ); }
      double d(   double r ){ return r/r02 + 2.*rbf(r)/r; }
      double d2(  double r ){ return ( 1./r02 + 2.*d(r)/r - 2.*rbf(r)/(r*r) ); }
  };


// gaussian radial basis function: y(r) = exp( -0.5*r*r/(r0*r0) )
   struct rbf_gaussian : rbf_f
  {
      double r0;
      double r02;

      rbf_gaussian( double scale=1. ) : r0(scale), r02( scale*scale ) {}

      double  rbf( double r ){ return exp( -0.5*r*r/r02 ); }
      double  d(   double r ){ return -rbf(r)*r/r0; }
      double  d2(  double r ){ return ( rbf(r)*( r*r/r0 -1. )/r0 ); }
  };


// biharmonic radial basis function: y(r) = r
   struct rbf_biharmonic : rbf_f
  {
      double rbf( double r ){ return r; }
      double d(    double r ){ return 1; }
      double d2(   double r ){ return 0; }
  };


// triharmonic radial basis function: y(r) = r^3
   struct rbf_triharmonic : rbf_f
  {
      double rbf( double r ){ return r*r*r; }
      double d(   double r ){ return 3*r*r; }
      double d2(  double r ){ return 6*r;   }
  };


# endif
