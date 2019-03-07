#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

#include "rbf.h"
#include "rbf_surf.h"

#define  DIMS  2

   struct quadratic_t
  {
      double   origin[2];
      double   coeff[2];

      vec_t    vec;

      quadratic_t( const double o[2], const double c[2] )
     {
         origin[0]=o[0];
         origin[1]=o[1];

         coeff[0]=c[0];
         coeff[1]=c[1];

         vec=vec_t(2);
     }

      double y( double x )
     {
         double x0=x-origin[0];
         return  coeff[0]*x0*x0 + coeff[1]*x0 + origin[1];
     }

      void tangent( double x, double dy[2] )
     {
         double x0=x-origin[0];
         dy[0]=1.;
         dy[1]= 2*coeff[0]*x0 + coeff[1];
         vec.unit( dy );
     }

      void normal( double x, double n[2] )
     {
         tangent( x, n );
         double t=n[0];
         n[0]=-n[1];
         n[1]= t;
     }
  };


   int main()
  {
      const int      npts  = 1500;
      const int      ntest = 15;
      const int      nplot = 51;

      const double   origin[2]={0.,0.};
      const double   coeff[2] ={1.,0.};

      const double   scale=0.01;

      const double   xmin=-3., xmax=3.;
      const double   ymin=-1., ymax=10.;

      int         i,j;
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
      quadratic_t   quad( origin, coeff );

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
      x   = new double*[npts ];
      n   = new double*[npts ];
      nt  = new double*[ntest];
      xt0 = new double*[ntest];
      xt1 = new double*[ntest];
      xt2 = new double*[ntest];

      for( i=0; i<npts;  i++ ){ x[  i]=NULL; x[  i] = new double[DIMS]; }
      for( i=0; i<npts;  i++ ){ n[  i]=NULL; n[  i] = new double[DIMS]; }
      for( i=0; i<ntest; i++ ){ nt[ i]=NULL; nt[ i] = new double[DIMS]; }
      for( i=0; i<ntest; i++ ){ xt0[i]=NULL; xt0[i] = new double[DIMS]; }
      for( i=0; i<ntest; i++ ){ xt1[i]=NULL; xt1[i] = new double[DIMS]; }
      for( i=0; i<ntest; i++ ){ xt2[i]=NULL; xt2[i] = new double[DIMS]; }

   // build point/normal cloud
      dx = ( xmax - xmin )/( npts -1 );
      s  = xmin;
      for( i=0; i<npts; i++ )
     {
         x[i][0]=s;
         x[i][1]=quad.y(s);
         quad.normal( x[i][0], n[i] );
         s+=dx;
     }

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


   // deallocate pointers
      for( i=0; i<npts;  i++ ){ delete[] x[  i]; x[  i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] n[  i]; n[  i]=NULL; }
      for( i=0; i<ntest; i++ ){ delete[] nt[ i]; nt[ i]=NULL; }
      for( i=0; i<ntest; i++ ){ delete[] xt0[i]; xt0[i]=NULL; }
      for( i=0; i<ntest; i++ ){ delete[] xt1[i]; xt1[i]=NULL; }
      for( i=0; i<ntest; i++ ){ delete[] xt2[i]; xt2[i]=NULL; }

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

