#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

#include "rbf.h"
#include "rbf_surf.h"

#define  DIMS  3
#define  pi    4.*atan(1.)
#define  pi2   8.*atan(1.)

   double cosineSpacing( int i, int np ){ return 0.5*( 1. - cos( i*pi/np ) ); }


   void PointsOnSharpBranch( int nCircle, int nPlane, int nBranch, double rmin, double rmax, double zmin, double zmax, double scale, double fillet, double phase, double **x, double **n, double *scales )
  {
      
      vec_t  vec(DIMS);
      double dtheta= pi2/nCircle;
      double ctheta, stheta, theta;
      double r,xp,yp,zp, dr, dz;
      int i, j;

      i=0;
      for( int k=0; k<nCircle; k++ )
     {
         theta=k*dtheta+phase;
         ctheta=cos(theta);
         stheta=sin(theta);
   
      // construct points on plane
         zp=0.;
         for( j=0; j<nPlane; j++ )
        {
//          dr = j*(rmax-rmin)/(nPlane-1);

            dr=rmax-rmin;
            dr=dr*cosineSpacing( j, nPlane-1 );

            r  = rmin + fillet + dr;
            xp = r*ctheta;
            yp = r*stheta;
   
            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;
   
            n[i][0]=0.;
            n[i][1]=0.;
            n[i][2]=1.;

            vec.unit( n[i] );

            scales[i]=std::min( scale, r-rmin );
   
            i++;
        }
   
      // construct points on cylinder
         r = rmin;
         xp=r*ctheta;
         yp=r*stheta;
         for( j=0; j<nBranch; j++ )
        {
//          dz = j*(zmax-zmin)/(nBranch-1);

            dz=zmax-zmin;
            dz=dz*cosineSpacing( j, nBranch-1 );

            zp = zmax - fillet - j*dz;
   
            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;
   
            n[i][0]=-xp;
            n[i][1]=-yp;
            n[i][2]=0.;
   
            vec.unit( n[i] );

            scales[i]=std::min( scale, abs( zmax-zp ) );
            scales[i]=std::min( scales[i], rmin );
   
            i++;
        }
     }

  }

   int main()
  {
      const int      nCircle = 12;
      const int      nBranch = 8;
      const int      nPlane  = 8;
      const int      npts = nCircle*( nBranch + nPlane );

      const int      nCircleTest = 5;
      const int      nBranchTest = 3;
      const int      nPlaneTest  = 3;
      const int      nTest = nCircleTest*( nBranchTest + nPlaneTest );

      const double   scale=0.025;

      const double   rmin=1., rmax= 2.;
      const double   zmin=-1., zmax=0.;
      const double   fillet=0.05;

      int					  i;
      double       **x=NULL;
      double       **n=NULL;
      double      **nt=NULL;
      double     **xt0=NULL;
      double     **xt1=NULL;
      double     **xt2=NULL;
      double   *scales=NULL;
      double   *scalesTest=NULL;

      double   r,dr;
      double   accuracy=10e-10;
      double   *u=NULL, *v=NULL, *w=NULL;

      vec_t         vec(DIMS);
      rbf_f         *rbf=NULL;

//    rbf= new rbf_biharmonic();
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);


   // initialise random number generator
      srand( 42 );
      std::cout << std::scientific;


   // allocate arrays
     {
      u   = new double [DIMS];
      v   = new double [DIMS];
      w   = new double [DIMS];
      x   = new double*[npts ];
      n   = new double*[npts ];
      nt  = new double*[nTest];
      xt0 = new double*[nTest];
      xt1 = new double*[nTest];
      xt2 = new double*[nTest];
      scales=new double[npts];
      scalesTest=new double[npts];

      for( i=0; i<npts;  i++ ){ x[  i]=NULL; x[  i] = new double[DIMS]; }
      for( i=0; i<npts;  i++ ){ n[  i]=NULL; n[  i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ nt[ i]=NULL; nt[ i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt0[i]=NULL; xt0[i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt1[i]=NULL; xt1[i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt2[i]=NULL; xt2[i] = new double[DIMS]; }
     }


   // construct points and normals
      PointsOnSharpBranch( nCircle, nPlane, nBranch, rmin, rmax, zmin, zmax, scale, fillet, 0., x, n, scales );

   // construct test points and normals offset around circle by r so not hitting
      PointsOnSharpBranch( nCircleTest, nPlaneTest, nBranchTest, rmin, rmax, zmin, zmax, scale, 2*fillet, pi2*0.13, xt0, nt, scalesTest );

   // displace test points off surface along normal
   // xt1 = off*nt + xt0
      const double off=1.0;
      for( i=0; i<nTest; i++ ){ vec.muladd( xt1[i], off*scalesTest[i], nt[i], xt0[i] ); }


   // build implicit surface interpolation
      rbf_surf isurf( DIMS, npts, x, n, scales, accuracy, rbf );
      i = isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // write points to file
     {
      std::ofstream ofile;
      ofile.open( "surface.dat" );

      std::ofstream ofileTest;
      ofileTest.open( "surfaceTest.dat" );

      for( i=0; i<nTest; i++ )
     {
         ofileTest << xt0[i][0] << " " << xt0[i][1] << " " << xt0[i][2] << " " << isurf.F( xt0[i] ) << std::endl;
     }


      for( i=0; i<npts; i++ )
     {
         ofile << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << isurf.F( x[i] ) << std::endl;
     }

      ofile.close();
      ofileTest.close();
     }


   // build projections
      for( i=0; i<nTest; i++ ){ std::cout << off*scalesTest[i] << ", "; isurf.project( xt1[i], xt2[i] ); }

   // print projection errors
      std::cout << "|error|, F(x), error vector " << std::endl;
      for( i=0; i<nTest; i++ )
     {
         vec.sub( v, xt0[i], xt2[i] );
         dr=vec.length( v );// / ( off*scalesTest[i] );
         r=isurf.F( xt2[i] );

         std::cout << dr << ", " << r << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }


   // deallocate pointers
     {
      for( i=0; i<npts;  i++ ){ delete[] x[  i]; x[  i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] n[  i]; n[  i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] nt[ i]; nt[ i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt0[i]; xt0[i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt1[i]; xt1[i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt2[i]; xt2[i]=NULL; }

      delete[]  u;    u=NULL;
      delete[]  v;    v=NULL;
      delete[]  w;    w=NULL;
      delete[]  x;    x=NULL;
      delete[]  n;    n=NULL;
      delete[] nt;   nt=NULL;
      delete[] xt0; xt0=NULL;
      delete[] xt1; xt1=NULL;
      delete[] xt2; xt2=NULL;
      delete[] scales; scales=NULL;
      delete[] scalesTest; scalesTest=NULL;

      delete  rbf;  rbf=NULL;
     }

      return 0;
  }

