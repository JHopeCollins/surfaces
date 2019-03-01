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
      virtual void normal(  double s, double n[2] )
     {
         double t[2];
         tangent( s, t );
         n[0] = -t[1];
         n[1] =  t[0];
         unit( n );
         return;
     }

   // q is projection of point p onto surface.
      virtual void project( double p[2], double q[2] ) = 0;

   // returns distance from point p to closest point on surface
      inline double distance( double p[2] )
     {
         double   q[2];

         project( p, q );

         q[0]=p[0]-q[0];
         q[1]=p[1]-q[1];

         return length( q );

     }
  };


// a straight line
   struct line_t : surf2_t
  {
      double origin[2];
      double direction[2];

      line_t(){};
     ~line_t(){};

      line_t( const double o[2], const double d[2] )
     {
         origin[0] = o[0];
         origin[1] = o[1];
         direction[0] = d[0];
         direction[1] = d[1];
     }

      void x( double s, double x[2] )
     {
         x[0] = origin[0] + s*direction[0];
         x[1] = origin[1] + s*direction[1];
         return;
     }

      double s( double x[2] )
     {
         double p[2];
         project( x, p );

         double t = ( dot( p, direction ) > 0 ? 1 : -1 );

         p[0]-=origin[0];
         p[1]-=origin[1];

         return t*length(p);
     }

      void tangent( double s, double d[2] )
     {
         double sign = ( direction[0] > 0 ? 1 : -1 );

         d[0] = sign*direction[0];
         d[1] = sign*direction[1];
         return;
     }
      
      void normal( double s, double n[2] )
     {
         double t[2];
         tangent( s, t );
         n[0] = -t[1];
         n[1] =  t[0];
         unit( n );
         return;
     }

      void project( double p[2], double q[2] )
     {
         double t=0;

         normal(t,q);
         t= dot(p,q);

         q[0]*=t;
         q[1]*=t;

         q[0]=p[0]-q[0];
         q[1]=p[1]-q[1];
         return;
     }
  };

# endif
