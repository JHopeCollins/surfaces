#include <iostream>
#include <math.h>

#include "rbf.h"
#include "rbf_interp.h"

#define  DIMS  2

// z = a + bx + cy + dxy + ex^2 + fy^2
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

      double F(     double *p ){ return F(     p[0], p[1] ); }
      double dFx(   double *p ){ return dFx(   p[0], p[1] ); }
      double dFy(   double *p ){ return dFy(   p[0], p[1] ); }
      double d2Fxx( double *p ){ return d2Fxx( p[0], p[1] ); }
      double d2Fyy( double *p ){ return d2Fyy( p[0], p[1] ); }
      double d2Fxy( double *p ){ return d2Fxy( p[0], p[1] ); }
  };

   int main()
  {
      const int   npts=50;
      const int   ntest=4;

      const double   xmin=-1.0;
      const double   xmax= 1.0;
      const double   ymin=-1.0;
      const double   ymax= 1.0;

      const double   scale=0.1;

      int         i,j,k;
      int         npts2=npts*npts;
      int         ntest2=ntest*ntest;
      double     **pt=NULL;
      double     **xt=NULL;
      double     *val=NULL;
      double      p[2]={0.,0.};

      double      x,y;
      double      z, dx, dy, dxy, dxx, dyy;
      double      zi,dxi,dyi,dxyi,dxxi,dyyi;
      double      g[DIMS], h[DIMS*DIMS];

      rbf_f      *rbf=NULL;
      rbf_f    *testf=NULL;
      vec_t      vec(DIMS);

      quadratic_t    square(1.,2.,3.,4.,5.,6.);

//    rbf= new rbf_biharmonic;
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);

      testf= new rbf_gaussian(scale);

   // allocate arrays
      val = new double [npts2];
      pt  = new double*[npts2];
      xt  = new double*[ntest2];

      for( i=0; i<npts2;  i++ ){ pt[i]=NULL; pt[i]= new double[DIMS]; }
      for( i=0; i<ntest2; i++ ){ xt[i]=NULL; xt[i]= new double[DIMS]; }

      dx=(xmax-xmin)/(npts-1);
      dy=(ymax-ymin)/(npts-1);
   // initialise interpolation points
      for( i=0; i<npts; i++ )
     {
         x=xmin+i*dx;
         p[0]=x;
         for( j=0; j<npts; j++ )
        {
            y=ymin+j*dy;
            p[1]=y;

            k=indx(i,j,npts);
            vec.eq( pt[k], p );

            val[k] = square.F( p );
        }
     }

      dx= ( xmax - xmin )/(ntest-1.);
      dx= ( xmax - xmin - dx ) / (ntest-1.);
      dy= ( ymax - ymin )/(ntest-1.);
      dy= ( ymax - ymin - dy ) / (ntest-1.);
   // initialise test points
      for( i=0; i<ntest; i++ )
     {
         x=xmin+(0.5+i)*dx;
         p[0]=x;
         for( j=0; j<ntest; j++ )
        {
            y=ymin+(0.5+j)*dy;
            p[1]=y;

            k=indx(i,j,ntest);
            vec.eq( xt[k], p );
        }
     }

      rbf_interp  interpolation( DIMS, npts*npts, pt, val, rbf );
      interpolation.build_weights();

   // compare interpolation
      std::cout << "z,     zi,     dx,     dxi,    dy,    dyi,  dxx,   dxxi,  dyy,  dyyi,  dxy,  dxyi" << std::endl;
      std::cout << std::scientific;
      for( i=0; i<ntest*ntest; i++ )
     {
         vec.eq( p, xt[i] );

         interpolation.derivative( p, g );
         interpolation.hessian(    p, h );
         z   = square.F(     p );
         dx  = square.dFx(   p );
         dy  = square.dFy(   p );
         dxx = square.d2Fxx( p );
         dyy = square.d2Fyy( p );
         dxy = square.d2Fxy( p );

         zi   = interpolation( p );
         dxi  = g[0];
         dyi  = g[1];
         dxxi = h[0];
         dyyi = h[3];
         dxyi = h[1];

         printf( "%3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f  %3.2f \n",
                  z,   zi,    dx,   dxi,  dy,   dyi,  dxx,  dxxi, dyy,  dyyi, dxy,  dxyi );

//       std::cout << z    << ", " << zi   << ", ";
//       std::cout << dx   << ", " << dxi  << ", ";
//       std::cout << dy   << ", " << dyi  << ", ";
//       std::cout << dxx  << ", " << dxxi << ", ";
//       std::cout << dyy  << ", " << dyyi << ", ";
//       std::cout << dxy  << ", " << dxyi << ", " << std::endl;
     }

   // deallocate arrays
      for( i=0; i<npts*npts;   i++ ){ delete[] pt[i]; pt[i]=NULL; }
      for( i=0; i<ntest*ntest; i++ ){ delete[] xt[i]; xt[i]=NULL; }

      delete[]  pt;  pt=NULL;
      delete[]  xt;  xt=NULL;
      delete[] val; val=NULL;

      delete   rbf;   rbf=NULL;
      delete testf; testf=NULL;

      return 0;
  }
