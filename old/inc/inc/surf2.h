# ifndef SURF2_H
# define SURF2_H

#include "vec.h"

// basic 1D surfaces in 2D space

// virtual base surface class defining interface
   struct surf2_t : public virtual vec_t
  {
      surf2_t(){ dims=2; };
     ~surf2_t(){};
   // convert from 1D parameterisation to real space
      virtual void x( double s, double x[2] ) = 0;

   // s is arc length parameter of point of surface closest to x
      virtual double s( double x[2] ) = 0;

   // calculate normal and gradient of surface
      virtual void tangent( double s, double d[2] ) = 0;
      virtual void normal(  double s, double n[2] );

   // q is projection of point p onto surface.
      virtual void project( double p[2], double q[2] ) = 0;

   // returns distance from point p to closest point on surface
      double distance( double p[2] );
  };


// a straight line
   struct line_t : surf2_t
  {
      double origin[2];
      double direction[2];

      line_t(){};
     ~line_t(){};

      line_t( const double o[2], const double d[2] );

   // convert from 1D parameterisation to real space
      void x( double s, double x[2] );

   // s is arc length parameter of point of surface closest to x
      double s( double x[2] );

   // calculate normal and gradient of surface
      void tangent( double s, double d[2] );
      void normal(  double s, double n[2] );

   // q is projection of point p onto surface.
      void project( double p[2], double q[2] );
  };

# endif
