#include <pseudosurf.h>


   void pseudo3surf_t::tangent( double *p, double *t0, double *t1 )
  {
      t0[0]=1.;
      t0[1]=0.;
      t0[2]=dFx(p);

      t1[0]=0.;
      t1[1]=1.;
      t1[2]=dFy(p);

      unit( t0 );
      unit( t1 );
  }

   void pseudo3surf_t::normal( double *p, double *n )
  {
      double   *t0=NULL, *t1=NULL;

      t0=new double[dims];
      t1=new double[dims];

      tangent( p, t0, t1 );
      cross(   n, t0, t1 );
      unit(n);

      delete[] t0; t0=NULL;
      delete[] t1; t1=NULL;
  }


// z = a + bx + cy + dxy + ex^2 + fy^2
   double quadratic_t::F(     double x, double y ){ return a + b*x + c*y + d*x*y +   e*x*x +   f*y*y; }
   double quadratic_t::dFx(   double x, double y ){ return     b   +       d * y +2.*e*x;             }
   double quadratic_t::dFy(   double x, double y ){ return           c   + d*x             +2.*f*y;   }
   double quadratic_t::d2Fxx( double x, double y ){ return                        2.*e;               }
   double quadratic_t::d2Fyy( double x, double y ){ return                                  2.*f;     }
   double quadratic_t::d2Fxy( double x, double y ){ return d; }


// z = a*cos(2*pi*x*lx) * b*sin(2*pi*y*ly)
   double sinusoid_t::F(     double x, double y ){ return        a *cos( pi2*x*lx ) *      b *sin( pi2*y*ly ); }
   double sinusoid_t::dFx(   double x, double y ){ return  -pi2lxa *sin( pi2*x*lx ) *      b *sin( pi2*y*ly ); }
   double sinusoid_t::dFy(   double x, double y ){ return        a *cos( pi2*x*lx ) * pi2lyb *cos( pi2*y*ly ); }
   double sinusoid_t::d2Fxx( double x, double y ){ return  -pi2lxa2*cos( pi2*x*lx ) *      b *sin( pi2*y*ly ); }
   double sinusoid_t::d2Fyy( double x, double y ){ return       -a *cos( pi2*x*lx ) * pi2lyb2*sin( pi2*y*ly ); }
   double sinusoid_t::d2Fxy( double x, double y ){ return  -pi2lxa *sin( pi2*x*lx ) * pi2lyb *cos( pi2*y*ly ); }
