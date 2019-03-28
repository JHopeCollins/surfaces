#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

#include "rbf.h"
#include "vtx.h"
#include "rbf_surf.h"

#define  DIMS  2

   template< typename POS_2, typename POS_1 >
   struct quadratic_surf_t
  {
      POS_2      origin;
      double   coeff[2];

      quadratic_surf_t( POS_2 o, double c[2] )
     {
         origin[0]=o[0];
         origin[1]=o[1];

         coeff[0]=c[0];
         coeff[1]=c[1];
     }

      POS_1 y( POS_1 x )
     {
         POS_1 x0=x-origin[0];
         return  coeff[0]*x0*x0 + coeff[1]*x0 + origin[1];
     }

      void tangent( POS_1 x, POS_2 &dy )
     {
         POS_1 x0=x-origin[0];
         dy[0]=1.;
         dy[1]= 2*coeff[0]*x0 + coeff[1];
         dy/=length(dy);
     }

      void normal( POS_1 x, POS_2 &n )
     {
         tangent( x, n );
         POS_1 t=n[0];
         n[0]=-n[1];
         n[1]= t;
     }
  };


   int main()
  {
      const int      npts  = 150;
      const int      ntest = 15;
      const int      nplot = 51;

      vtx_t     origin(0.,0.,0.);
      double   coeff[2] ={1.,0.};

      const double   scale=0.01;

      const double   xmin=-3., xmax=3.;
      const double   ymin=-1., ymax=10.;

      int         i,j;
      vtx_t       *x=NULL;
      vtx_t       *n=NULL;
      vtx_t      *nt=NULL;
      vtx_t     *xt0=NULL;
      vtx_t     *xt1=NULL;
      vtx_t     *xt2=NULL;
      double    *scales=NULL;

      double    s, ds, off;
      double   dx, dy, xp,yp;
      double   accuracy=10e-8;
      vtx_t    v, w;

      quadratic_surf_t<vtx_t,double>  quad( origin, coeff );

   // initialise random number generator
      srand( 42 );

   // allocate arrays
     {
      x   = new vtx_t [npts ];
      n   = new vtx_t [npts ];
      nt  = new vtx_t [ntest];
      xt0 = new vtx_t [ntest];
      xt1 = new vtx_t [ntest];
      xt2 = new vtx_t [ntest];
      scales=new double[npts];
     }

   // build point/normal cloud
      dx = ( xmax - xmin )/( npts -1 );
      s  = xmin;
      for( i=0; i<npts; i++ )
     {
         x[i][0]=s;
         x[i][1]=quad.y(s);
         quad.normal( x[i][0], n[i] );
         s+=dx;

         scales[i]=scale;
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
   // x1[i] should project back to x0[i]
      for( i=0; i<ntest; i++ )
     {
      // random number in j*[0-1] to 2 decimal places
         s=off;
//       s = j*( rand() % 100 ) / 100; s*=off;

         xt1[i]=xt0[i]+s*nt[i];

         j*=-1;
     }

   // build implicit surface interpolation
      rbf_surf<2, vtx_t, double> surface( npts, x, n, scales, accuracy );
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
         v  = xt0[i] - xt2[i];
         ds = length( v )/off;

         s=surface.F( xt2[i] );

         std::cout << ds << ", " << s << ", (" << v[0] << ", " << v[1] << ")" << std::endl;
     }


   // deallocate pointers
      delete[]  x;    x=NULL;
      delete[]  n;    n=NULL;
      delete[] nt;   nt=NULL;
      delete[] xt0; xt0=NULL;
      delete[] xt1; xt1=NULL;
      delete[] xt2; xt2=NULL;
      delete[] scales; scales=NULL;

      return 0;
  }

