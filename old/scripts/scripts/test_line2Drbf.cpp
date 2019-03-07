#include <iostream>
#include <math.h>
#include <fstream>

#include "rbf.h"
#include "rbf_surf.h"
#include "surf2.h"

   int main()
  {
      const int   dims =2;
      const int   npts =100;
      const int   ntest=5;
      const int   nplot=51;

      const double      origin[2]={0.,0.};
      const double   direction[2]={1.,1.};

      const double   scale=0.1;
      const double   smin =0.0;
      const double   smax =5.0;

      const double   xmin=-1., xmax=6.;
      const double   ymin=-1., ymax=6.;

      int         i,j;
      double      **x=NULL;
      double      **n=NULL;
      double     **xt=NULL;

      double    s, ds;
      double   dx, dy, xp,yp;
      double   accuracy=0.01;
      double   *v=NULL, *w=NULL;

      rbf_f    *rbf=NULL;
      line_t   line( origin, direction );

//    rbf= new rbf_biharmonic;
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);


   // allocate arrays
      v  = new double [dims];
      w  = new double [dims];
      x  = new double*[npts ];
      n  = new double*[npts ];
      xt = new double*[ntest];

      for( i=0; i<npts;  i++ ){ x[ i]=NULL; x[ i] = new double[dims]; }
      for( i=0; i<npts;  i++ ){ n[ i]=NULL; n[ i] = new double[dims]; }
      for( i=0; i<ntest; i++ ){ xt[i]=NULL; xt[i] = new double[dims]; }

   // build point/normal cloud
      ds = ( smax - smin )/( npts -1 );
      s  = smin;
      for( i=0; i<npts; i++ )
     {
         line.x(      s, x[i] );
         line.normal( s, n[i] );
         s+=ds;
     }

   // points to test surface reconstruction at
      ds= ( smax - smin )/(ntest-1.);
      ds= ( smax - smin - ds ) / (ntest-1.);
      s= smin+0.5*ds;
      for( i=0; i<ntest; i++ ){ line.x( s+i*ds, xt[i] ); }


   // build implicit surface interpolation
      rbf_surf surface( dims, npts, x, n, scale, accuracy, rbf );
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


   // test distance function of points on line
      std::cout << "distance function at test points: " << std::endl;
      for( i=0; i<ntest; i++ )
     {
         ds=surface.distance( xt[i] );
         std::cout << ds << ", " << std::endl;
     }
      std::cout << std::endl;

      std::cout << "v" << std::endl;
      v[0]=0.; v[1]=2.;
      ds=surface.distance( v );
      std::cout << ds << std::endl;


   // test tangent vectors
      std::cout << "tangent vector at test points: " << std::endl;
      for( i=0; i<ntest; i++ )
     {
         surface.tangent( xt[i], v );
         std::cout << "(" << v[0] << "," << v[1] << ")" << std::endl;
     }
      std::cout << std::endl;


   // test normal vectors
      std::cout << "normal vector at test points: " << std::endl;
      for( i=0; i<ntest; i++ )
     {
         surface.normal( xt[i], v );
         std::cout << "(" << v[0] << "," << v[1] << ")" << std::endl;
     }
      std::cout << std::endl;


   // test projection onto surface

      std::cout << "project(0,2) ~> (1,1):" << std::endl;
      v[0]=0.; v[1]=2.;
      w[0]=0.; w[1]=0.;

//    std::cout << "distance at start:" << std::endl;
//    ds = surface.distance( v );
//    std::cout << ds << std::endl;

      std::cout << "projected point:" << std::endl;
      surface.project( v, w );
      std::cout << "(" << w[0] << "," << w[1] << ")" << std::endl;
      std::cout << std::endl;

//    std::cout << "distance at end:" << std::endl;
//    ds = surface.distance( w );
//    std::cout << ds << std::endl;

   // deallocate pointers
      for( i=0; i<npts;  i++ ){ delete[] x[ i]; x[ i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] n[ i]; n[ i]=NULL; }
      for( i=0; i<ntest; i++ ){ delete[] xt[i]; xt[i]=NULL; }

      delete[]  v;   v=NULL;
      delete[]  w;   w=NULL;
      delete[]  x;   x=NULL;
      delete[]  n;   n=NULL;
      delete[] xt;  xt=NULL;

      delete  rbf; rbf=NULL;

      return 0;
  }
