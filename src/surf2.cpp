#include "surf2.h"

// basic 1D surfaces in 2D space

// returns distance from point p to closest point on surface
   double surf2_t::distance( double p[2] )
  {
      double   q[2];

      project( p, q );

      q[0]=p[0]-q[0];
      q[1]=p[1]-q[1];

      return length( q );

  }

   void surf2_t::normal( double s, double n[2] )
     {
         double t[2];
         tangent( s, t );
         n[0] = -t[1];
         n[1] =  t[0];
         unit( n );
         return;
     }

// a straight line
   line_t::line_t( const double o[2], const double d[2] )
  {
      origin[0] = o[0];
      origin[1] = o[1];
      direction[0] = d[0];
      direction[1] = d[1];
  }

   void line_t::x( double s, double x[2] )
  {
      x[0] = origin[0] + s*direction[0];
      x[1] = origin[1] + s*direction[1];
      return;
  }

   double line_t::s( double x[2] )
  {
      double p[2];
      project( x, p );

      double t = ( dot( p, direction ) > 0 ? 1 : -1 );

      p[0]-=origin[0];
      p[1]-=origin[1];

      return t*length(p);
  }

   void line_t::tangent( double s, double d[2] )
  {
      double sign = ( direction[0] > 0 ? 1 : -1 );

      d[0] = sign*direction[0];
      d[1] = sign*direction[1];
      return;
  }

   void line_t::normal( double s, double n[2] )
  {
      double t[2];
      tangent( s, t );
      n[0] = -t[1];
      n[1] =  t[0];
      unit( n );
      return;
  }

   void line_t::project( double p[2], double q[2] )
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

