#include <iostream>
#include <math.h>

#include "vtx.h"
#include "rbf_interp.h"
#include "pseudosurf.h"

#define  DIMS  2

// z = a + bx + cy + dxy + ex^2 + fy^2
/*
   struct quadratic_t
  {
      double a, b, c, d, e, f;

      quadratic_t( double aa, double bb, double cc, double dd, double ee, double ff )
            : a(aa), b(bb), c(cc), d(dd), e(ee), f(ff) {};

      double F(     double x, double y ){ return a + b*x + c*y + d*x*y +   e*x*x + f*y*y; }
      double dFx(   double x, double y ){ return     b   +       d * y +2.*e*x;           }
      double dFy(   double x, double y ){ return           c   + d*x   +        2.*f*y;   }
      double d2Fxx( double x, double y ){ return                        2.*e;             }
      double d2Fyy( double x, double y ){ return                                2.*f;     }
      double d2Fxy( double x, double y ){ return d; }

      double F(     vtx_t p ){ return F(     p[0], p[1] ); }
      double dFx(   vtx_t p ){ return dFx(   p[0], p[1] ); }
      double dFy(   vtx_t p ){ return dFy(   p[0], p[1] ); }
      double d2Fxx( vtx_t p ){ return d2Fxx( p[0], p[1] ); }
      double d2Fyy( vtx_t p ){ return d2Fyy( p[0], p[1] ); }
      double d2Fxy( vtx_t p ){ return d2Fxy( p[0], p[1] ); }
  };
*/

   int main()
  {
      const int   npts=10;
      const int   ntest=4;

      const double   xmin=-1.0;
      const double   xmax= 1.0;
      const double   ymin=-1.0;
      const double   ymax= 1.0;

      int         i,j,k;
      int         npts2=npts*npts;
      int         ntest2=ntest*ntest;
      vtx_t       *pt=NULL;
      vtx_t       *xt=NULL;
      double     *val=NULL;
      vtx_t              p;

      double      x,y;
      double      z, dx, dy, dxy, dxx, dyy;
      double      zi,dxi,dyi,dxyi,dxxi,dyyi;
      double      g[DIMS], h[DIMS*DIMS];

      quadratic_t    q(1.,2.,3.,4.,5.,6.);
//    std::cout << q.a << ", "
//              << q.b << ", "
//              << q.c << ", "
//              << q.d << ", "
//              << q.e << ", "
//              << q.f << ", " << std::endl;

//    std::cout << q.F( 19., 20056. ) << std::endl;


   // allocate arrays
      val = new double [npts2];
      pt  = new vtx_t  [npts2];
      xt  = new vtx_t  [ntest2];

      dx=(xmax-xmin)/(npts-1);
      dy=(ymax-ymin)/(npts-1);
   // initialise interpolation points
      p[2]=0.;
      for( i=0; i<npts; i++ )
     {
         x=xmin+i*dx;
         p[0]=x;
         for( j=0; j<npts; j++ )
        {
            y=ymin+j*dy;
            p[1]=y;

            k=indx(i,j,npts);
            pt[k]=p;

            val[k] = q.F( p[0], p[1] );
        }
     }

      dx= ( xmax - xmin )/(ntest-1.);
      dx= ( xmax - xmin - dx ) / (ntest-1.);
      dy= ( ymax - ymin )/(ntest-1.);
      dy= ( ymax - ymin - dy ) / (ntest-1.);
   // initialise test points
      p[2]=0.;
      for( i=0; i<ntest; i++ )
     {
         x=xmin+(0.5+i)*dx;
         p[0]=x;
         for( j=0; j<ntest; j++ )
        {
            y=ymin+(0.5+j)*dy;
            p[1]=y;

            k=indx(i,j,ntest);
            xt[k]=p;
//          std::cout << "p= (" << p[0] << ", " << p[1] << ", " << p[2] << ")  |  f(p)= " << val[k] << std::endl;
        }
     }

      rbf_interp<2, vtx_t, double>  interpolation( npts*npts, pt, val );
      i = interpolation.build_weights();
      std::cout << "info: " << i << std::endl;

   // compare interpolation
      std::cout << "z,     zi,     dx,   dxi,   dy,  dyi, dxx,  dxxi, dyy,  dyyi, dxy, dxyi" << std::endl;
      std::cout << std::scientific;
      for( i=0; i<ntest*ntest; i++ )
     {
         p=xt[i];

         interpolation.derivative( p, g );
         interpolation.hessian(    p, h );
         z   = q.F(     p[0], p[1] );
         dx  = q.dFx(   p[0], p[1] );
         dy  = q.dFy(   p[0], p[1] );
         dxx = q.d2Fxx( p[0], p[1] );
         dyy = q.d2Fyy( p[0], p[1] );
         dxy = q.d2Fxy( p[0], p[1] );

         zi   = interpolation( p );
         dxi  = g[0];
         dyi  = g[1];
         dxxi = h[0];
         dyyi = h[3];
         dxyi = h[1];

         printf( "%3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f \n",
                  z,   zi,    dx,   dxi,  dy,   dyi,  dxx,  dxxi, dyy,  dyyi, dxy,  dxyi );
     }

   // deallocate arrays
      delete[]  pt;  pt=NULL;
      delete[]  xt;  xt=NULL;
      delete[] val; val=NULL;

      return 0;
  }
