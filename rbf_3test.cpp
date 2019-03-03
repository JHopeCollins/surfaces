#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

#include "rbf.h"
#include "rbf_interp.h"
#include "pseudosurf.h"

#define  DIMS  3

   int main()
  {
      const int      npts  = 30;
      const int      npts2 = npts*npts;
      const int      ntest = 4;
      const int      ntest2= ntest*ntest;
      const int      nplot = 51;


      const double   scale=0.05;

      const double   xmin=-3., xmax=3.;
      const double   ymin=-3., ymax=3.;
      const double   zmin=-0., zmax=30.;

      int             i,j,k,np;
      double       **x=NULL;
      double       **n=NULL;
      double      **nt=NULL;
      double     **xt0=NULL;
      double     **xt1=NULL;
      double     **xt2=NULL;

      double           s, ds, off;
      double   dx,dy,dz, xp,yp,zp;
      double       accuracy=10e-10;
      double     *u=NULL, *v=NULL, *w=NULL;

      vec_t         vec(DIMS);
      rbf_f          *rbf=NULL;
      pseudo3surf_t *surf=NULL;

//    rbf= new rbf_biharmonic();
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);

//    surf= new quadratic_t(1.,1.,1.,1.,1.,1.);
      surf= new sinusoid_t( 1.,0.20,1.,0.20);

   // initialise random number generator
      srand( 42 );
      std::cout << std::scientific;

   // allocate arrays
     {
      u   = new double [DIMS];
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
     }

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

            x[k][2]=(*surf).F( x[k] );

           (*surf).normal( x[k], n[k] );
        }
     }

   // points to test surface reconstruction at
      dx= ( xmax - xmin )/(ntest-1.);
      dx= ( xmax - xmin - dx ) / (ntest-1.);
      dy= ( ymax - ymin )/(ntest-1.);
      dy= ( ymax - ymin - dy ) / (ntest-1.);
      for( i=0; i<ntest; i++ )
     {
         for( j=0; j<ntest; j++ )
        {
            k=index(i,j,ntest);

            xt0[k][0]=xmin+i*dx;
            xt0[k][1]=ymin+j*dy;

            xt0[k][2]=(*surf).F(xt0[k]);

           (*surf).normal( xt0[k], nt[k] );
        }
     }

      j=1;
      off=0.1;
   // build points xt1 to test from
   // points are displaced normal off the surface by a random amount
   // x1[i] should project back to x0[i]
      for( i=0; i<ntest2; i++ )
     {
      // random number in j*[0-1] to 2 decimal places
         s=off;
//       s = j*( rand() % 100 ) / 100; s*=off;

         vec.mul( v, s, nt[i] );

         vec.add( xt1[i], xt0[i], v );

         j*=-1;
     }


   // build implicit surface interpolation
      rbf_surf isurf( DIMS, npts2, x, n, scale, accuracy, rbf );
      i = isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // write surface coordinates for plotting
     {
      std::ofstream ofile;
      ofile.open( "surface.dat" );

      np = npts;
      dx = ( xmax - xmin )/np;
      dy = ( ymax - ymin )/np;

      xp=xmin;
      for( i=0; i<np; i++ )
     {
         yp=ymin;
         for( j=0; j<np; j++ )
        {
            v[0]=xp;
            v[1]=yp;
            v[2]=(*surf).F( v );

            ofile << v[0] << " " << v[1] << " " << v[2] << std::endl;

            yp+=dy;
        }
         xp+=dx;
         ofile << std::endl;
     }
      ofile.close();
     }


/*
   // write distance function for plotting
     {
      std::ofstream ofile;
      ofile.open( "distance.dat" );

      dx = ( xmax - xmin )/( nplot-1 );
      dy = ( ymax - ymin )/( nplot-1 );
      dz = ( zmax - zmin )/( nplot-1 );

      xp=xmin;
      for( i=0; i<nplot; i++ )
     {
         yp=ymin;
         for( j=0; j<nplot; j++ )
        {
            zp=zmin;
            for( k=0; k<nplot; k++ )
           {
               v[0]=xp;
               v[1]=yp;
               v[2]=zp;

               ds=isurf.F( v );

               ofile << v[0] << " " << v[1] << " " << v[2] << " " << ds << std::endl;

               zp+=dz;
           }
            yp+=dy;
        }
         xp+=dx;
         ofile << std::endl;
     }
      ofile.close();
     }
*/

   // points to test surface normal at
      std::cout << "x, n, ni" << std::endl;
      dx= ( xmax - xmin )/(ntest-1.);
      dx= ( xmax - xmin - dx ) / (ntest-1.);
      dy=0;
      for( i=0; i<ntest2; i++ )
     {
         std::cout << "(" << xt0[i][0] << ", " << xt0[i][1] << ", " << xt0[i][2] << ") ";

        (*surf).normal( xt0[i], u );
         std::cout << "(" << v[0] << ", " << v[1] << ", " <<v[2] << ") ";
         isurf.dF( xt0[i], v );
         std::cout << "(" << v[0] << ", " << v[1] << ", " <<v[2] << ")  ";

         std::cout << vec.radius( u, v );

         std::cout << std::endl;
     }

   // build projections
      for( i=0; i<ntest2; i++ ){ isurf.project( xt1[i], xt2[i] ); }

      std::cout << "projection length errors and distance function" << std::endl;
      for( i=0; i<ntest2; i++ )
     {
         vec.sub( v, xt0[i], xt2[i] );
         ds=vec.length( v )/off;
         s=isurf.F( xt2[i] );

         std::cout << ds << ", " << s << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }


   // deallocate pointers
     {
      for( i=0; i<npts2;  i++ ){ delete[] x[  i]; x[  i]=NULL; }
      for( i=0; i<npts2;  i++ ){ delete[] n[  i]; n[  i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] nt[ i]; nt[ i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt0[i]; xt0[i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt1[i]; xt1[i]=NULL; }
      for( i=0; i<ntest2; i++ ){ delete[] xt2[i]; xt2[i]=NULL; }

      delete[]  u;    u=NULL;
      delete[]  v;    v=NULL;
      delete[]  w;    w=NULL;
      delete[]  x;    x=NULL;
      delete[]  n;    n=NULL;
      delete[] nt;   nt=NULL;
      delete[] xt0; xt0=NULL;
      delete[] xt1; xt1=NULL;
      delete[] xt2; xt2=NULL;

      delete  rbf;  rbf=NULL;
      delete surf; surf=NULL;
     }

      return 0;
  }

