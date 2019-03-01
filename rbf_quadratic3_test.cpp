#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

#include "rbf.h"
#include "rbf_interp.h"

#define  DIMS  3

// z = a + bx + cy + dxy + ex^2 + fy^2
   struct quadratic_t : public vec_t
  {
      double a, b, c, d, e, f;

      quadratic_t( double aa, double bb, double cc, double dd, double ee, double ff )
            : a(aa), b(bb), c(cc), d(dd), e(ee), f(ff){ dims=DIMS; }

      inline double F(     double x, double y ){ return a + b*x + c*y + d*x*y +   e*x*x + f*y*y; }
      inline double dFx(   double x, double y ){ return     b   +       d * y +2.*e*x;           }
      inline double dFy(   double x, double y ){ return           c   + d*x   +        2.*f*y;   }
      inline double d2Fxx( double x, double y ){ return                        2.*e;             }
      inline double d2Fyy( double x, double y ){ return                                2.*f;     }
      inline double d2Fxy( double x, double y ){ return d; }

      inline double F(     double *p ){ return F(     p[0], p[1] ); }
      inline double dFx(   double *p ){ return dFx(   p[0], p[1] ); }
      inline double dFy(   double *p ){ return dFy(   p[0], p[1] ); }
      inline double d2Fxx( double *p ){ return d2Fxx( p[0], p[1] ); }
      inline double d2Fyy( double *p ){ return d2Fyy( p[0], p[1] ); }
      inline double d2Fxy( double *p ){ return d2Fxy( p[0], p[1] ); }

      void tangent( double *p, double *t0, double *t1 )
     {
         t0[0]=1.;
         t0[1]=0.;
         t0[2]=dFx(p);

         t1[0]=0.;
         t1[1]=1.;
         t1[2]=dFy(p);
     }

      void normal( double *p, double *n )
     {
         double   t0[DIMS], t1[DIMS];

         tangent( p, t0, t1 );
         cross( n, t0, t1 );
     }

  };


   int main()
  {
      const int      npts  = 50;
      const int      npts2 = npts*npts;
      const int      ntest = 5;
      const int      ntest2= ntest*ntest;
      const int      nplot = 51;


      const double   scale=0.01;

      const double   xmin=-1., xmax=1.;
      const double   ymin=-1., ymax=1.;

      int             i,j,k;
      double       **x=NULL;
      double       **n=NULL;
      double      **nt=NULL;
      double     **xt0=NULL;
      double     **xt1=NULL;
      double     **xt2=NULL;

      double    s, ds, off;
      double   dx, dy, xp,yp;
      double   accuracy=10e-8;
      double   *v=NULL, *w=NULL;

      vec_t         vec(DIMS);
      rbf_f        *rbf=NULL;
      quadratic_t   surf(1.,2.,3.,4.,5.,6.);
      surf.construct_levi_civita();

      rbf= new rbf_biharmonic;
//    rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);


   // initialise random number generator
      srand( 42 );

   // allocate arrays
      v   = new double [DIMS];
      w   = new double [DIMS];
      x   = new double*[npts2 ];
      n   = new double*[npts2 ];
      nt  = new double*[ntest2];
      xt0 = new double*[ntest2];
      xt1 = new double*[ntest2];
      xt2 = new double*[ntest2];

      for( i=0; i<npts2;  i++ ){ x[  i]=NULL; x[  i] = new double[DIMS]; }
      for( i=0; i<npts2;  i++ ){ n[  i]=NULL; n[  i] = new double[DIMS]; }
      for( i=0; i<ntest2; i++ ){ nt[ i]=NULL; nt[ i] = new double[DIMS]; }
      for( i=0; i<ntest2; i++ ){ xt0[i]=NULL; xt0[i] = new double[DIMS]; }
      for( i=0; i<ntest2; i++ ){ xt1[i]=NULL; xt1[i] = new double[DIMS]; }
      for( i=0; i<ntest2; i++ ){ xt2[i]=NULL; xt2[i] = new double[DIMS]; }

   // build point/normal cloud
      dx=(xmax-xmin)/(npts-1);
      dy=(ymax-ymin)/(npts-1);
      for( i=0; i<npts; i++ )
     {

         for( j=0; j<npts; j++ )
        {
            k=index(i,j,npts);

            x[k][0]=xmin+i*dx;
            x[k][1]=ymin+j*dy;

            x[k][2]=surf.F(x[k]);

            surf.normal( x[k], n[k] );
        }
     }

/*
   // points to test surface reconstruction at
      ds= ( xmax - xmin )/(ntest-1.);
      ds= ( xmax - xmin - ds ) / (ntest-1.);
      s= xmin+0.5*ds;
      for( i=0; i<ntest; i++ )
     {
         xt0[i][0]=s;
         xt0[i][1]=quad.y(s);
         quad.normal( xt0[i][0], nt[i] );
         s+=ds;
     }

      j=1;
      off=0.5;
   // build points xt1 to test from
   // points are displaced normal off the surface by a random amount
   // x1[i] should project back to x0[i]
      for( i=0; i<ntest; i++ )
     {
      // random number in j*[0-1] to 2 decimal places
         s=off;
//       s = j*( rand() % 100 ) / 100; s*=off;

         vec.mul( v, s, nt[i] );

         vec.add( xt1[i], xt0[i], v );

         j*=-1;
     }


   // build implicit surface interpolation
      rbf_surf surface( DIMS, npts, x, n, scale, accuracy, rbf );
      i = surface.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // write distance function for plotting
     {
      std::ofstream ofile;
      ofile.open( "distance.dat" );

      dx = ( xmax - xmin )/( nplot-1 );
      dy = ( ymax - ymin )/( nplot-1 );

      xp=xmin;
      for( i=0; i<nplot; i++ )
     {
         yp=ymin;
         for( j=0; j<nplot; j++ )
        {
            v[0]=xp;
            v[1]=yp;

            ds=surface.F( v );

            ofile << v[0] << " " << v[1] << " " << ds << std::endl;

            yp+=dy;
        }
         xp+=dx;
         ofile << std::endl;
     }
      ofile.close();
     }


   // build projections
      for( i=0; i<ntest; i++ ){ surface.project( xt1[i], xt2[i] ); }

      std::cout << "projection length errors and distance function" << std::endl;
      std::cout << std::scientific;
      for( i=0; i<ntest; i++ )
     {
         vec.sub( v, xt0[i], xt2[i] );
         ds=vec.length( v )/off;
         s=surface.F( xt2[i] );

         std::cout << ds << ", " << s << ", (" << v[0] << ", " << v[1] << ")" << std::endl;
     }


*/

   // deallocate pointers
      for( i=0; i<npts2;  i++ ){ delete[] x[  i]; x[  i]=NULL; }
      for( i=0; i<npts2;  i++ ){ delete[] n[  i]; n[  i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] nt[ i]; nt[ i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt0[i]; xt0[i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt1[i]; xt1[i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt2[i]; xt2[i]=NULL; }

      delete[]  v;    v=NULL;
      delete[]  w;    w=NULL;
      delete[]  x;    x=NULL;
      delete[]  n;    n=NULL;
      delete[] nt;   nt=NULL;
      delete[] xt0; xt0=NULL;
      delete[] xt1; xt1=NULL;
      delete[] xt2; xt2=NULL;

      delete  rbf; rbf=NULL;

      return 0;
  }

