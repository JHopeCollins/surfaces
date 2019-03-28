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
      double costheta, sintheta, theta;
      double r,xp,yp,zp, dr, dz;
      int i, j;

      i=0;
      for( int k=0; k<nCircle; k++ )
     {
         theta=k*dtheta+phase;
         costheta=cos(theta);
         sintheta=sin(theta);

      // construct points on plane
         zp=zmax;
         for( j=0; j<nPlane; j++ )
        {
            dr = j*(rmax-rmin)/(nPlane-1);

//          dr=rmax-rmin;
//          dr=dr*cosineSpacing( j, nPlane-1 );

            r  = rmin + fillet + dr;
            if( j==0 ){ r+=0.1*dr; }
            xp = r*costheta;
            yp = r*sintheta;

            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;

            n[i][0]=0.;
            n[i][1]=0.;
            n[i][2]=1.;

            vec.unit( n[i] );

            scales[i]=std::min( scale,     r-rmin );
            scales[i]=std::min( scales[i], fillet );

            i++;
        }

      // construct points on cylinder
         r = rmin;
         xp=r*costheta;
         yp=r*sintheta;
         for( j=0; j<nBranch; j++ )
        {
            dz = j*(zmax-zmin)/(nBranch-1);

//          dz=zmax-zmin;
//          dz=dz*cosineSpacing( j, nBranch-1 );

            zp = zmax - fillet - j*dz;
            if( j==0 ){ zp-=0.1*dz; }

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

/*
   void PointsOnFilletBranch( int nCircle, int nPlane, int nBranch, int nFillet, double rmin, double rmax, double zmin, double zmax, double scale, double rfillet, double phase, double **x, double **n, double *scales )
  {
      int i=nCircle*( nBranch + nPlane ); // number of points generated by PointsOnSharpBranch()
      printf( "iFillet: %d\n", i );

      int j,k;

      double dtheta= pi2/nCircle;      // theta is angle around branch
      double costheta, sintheta, theta;

      double dpsi= 0.5*pi/(nFillet-1);     // psi is angle around fillet (0 at branch, pi/2 at plane)
      double cospsi, sinpsi, psi;

      double rtorus=rmin+rfillet;      // radius of circle at centre of fillet torus
      double ractual;                  // radius of current point

      double ztorus=zmax-rfillet;       // height of torus centreline
      double xp,yp,zp;                 // coordinates of point

      vec_t  vec(DIMS);

      PointsOnSharpBranch( nCircle, nPlane, nBranch, rmin, rmax, zmin, zmax, scale, rfillet, phase, x, n, scales );

      for( k=0; k<nCircle; k++ )
     {
         theta=k*dtheta+phase;
         costheta=cos(theta);
         sintheta=sin(theta);

      //points on fillet
         for( j=0; j<nFillet; j++ )
        {
            psi=j*dpsi;
            printf( "iteration %d, psi=%f\n", j, psi );
            sinpsi=sin(psi);
            cospsi=cos(psi);

            ractual= rtorus - rfillet*cospsi;

            xp = ractual*costheta;
            yp = ractual*sintheta;
            zp = ztorus + rfillet*sinpsi;

            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;

            n[i][0]= -cospsi*costheta;
            n[i][1]= -cospsi*sintheta;
            n[i][2]=  sinpsi;

            vec.unit( n[i] );

            scales[i]=0.2*rfillet;

            i++;
        }
     }
  }
*/

   void PointsOnFilletBranch( int nCircle, int nPlane, int nBranch, int nFillet, double rmin, double rmax, double zmin, double zmax, double scale, double rfillet, double phase, double **x, double **n, double *scales )
  {

      int i,j,k;

      double dtheta= pi2/nCircle;      // theta is angle around branch
      double costheta, sintheta, theta;

      double dpsi= 0.5*pi/(nFillet-1);     // psi is angle around fillet (0 at branch, pi/2 at plane)
      double cospsi, sinpsi, psi;

      double rtorus=rmin+rfillet;      // radius of circle at centre of fillet torus
      double ractual;                  // radius of current point

      double ztorus=zmax-rfillet;      // height of torus centreline
      double xp,yp,zp;                 // coordinates of point
      double r, dz, dr;

      vec_t  vec(DIMS);

      i=0;
      for( k=0; k<nCircle; k++ )
     {
         theta=k*dtheta+phase;
         costheta=cos(theta);
         sintheta=sin(theta);

      // construct points on cylinder
         r = rmin;
         xp=r*costheta;
         yp=r*sintheta;
         dz = (ztorus-zmin)/(nBranch);
         for( j=1; j<nBranch+1; j++ )
        {
            zp = ztorus - j*dz;
//          if( j==0 ){ zp-=0.1*dz; }

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

      //points on fillet
         for( j=0; j<nFillet; j++ )
        {
            psi=j*dpsi;
            sinpsi=sin(psi);
            cospsi=cos(psi);

            ractual= rtorus - rfillet*cospsi;

            xp = ractual*costheta;
            yp = ractual*sintheta;
            zp = ztorus + rfillet*sinpsi;

            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;

            n[i][0]= -cospsi*costheta;
            n[i][1]= -cospsi*sintheta;
            n[i][2]=  sinpsi;

            vec.unit( n[i] );

            scales[i]=0.25*rfillet;

            i++;
        }

      // construct points on plane
         zp=zmax;
         dr = (rmax-rmin)/(nPlane);
         for( j=1; j<nPlane+1; j++ )
        {

            r  = rmin + rfillet + j*dr;
//          if( j==0 ){ r+=0.1*dr; }
            xp = r*costheta;
            yp = r*sintheta;

            x[i][0]=xp;
            x[i][1]=yp;
            x[i][2]=zp;

            n[i][0]=0.;
            n[i][1]=0.;
            n[i][2]=1.;

            vec.unit( n[i] );

            scales[i]=std::min( scale,      r-rmin );
            scales[i]=std::min( scales[i], rfillet );

            i++;
        }
     }
  }

   int main()
  {
      const int      nCircle = 8;
      const int      nBranch = 4;
      const int      nPlane  = 4;
      const int      nFillet = 6;
//    const int      npts = nCircle*( nBranch + nPlane );
      const int      npts  = 2*nCircle*( nBranch + nPlane + nFillet );
      const int      npts2 =   nCircle*( nBranch + nPlane + nFillet );

      const int      nCircleTest = 4;
      const int      nBranchTest = 3;
      const int      nPlaneTest  = 3;
      const int      nFilletTest = 3;
//    const int      nTest = nCircleTest*( nBranchTest + nPlaneTest );
      const int      nTest  = 2*nCircleTest*( nBranchTest + nPlaneTest + nFilletTest );
      const int      nTest2 =   nCircleTest*( nBranchTest + nPlaneTest + nFilletTest );

      const double   scale=0.25;

      const double   rmin=1., rmax= 2.;
      const double   zmin=0.1, zmax= 1.;
      const double   fillet=0.50;

      int					  i;
      double       **x=NULL;  // interpolated points
      double      **x1=NULL;  // displaced interpolated points
      double      **x2=NULL;  // projection of displaced points

      double       **n=NULL;
      double      **nt=NULL;
      double     **xt0=NULL;
      double     **xt1=NULL;
      double     **xt2=NULL;
      double   *scales=NULL;
      double   *scalesTest=NULL;

      double		r,dr;
      double      accuracy=10e-10;
      double     *u=NULL, *v=NULL, *w=NULL;

      vec_t         vec(DIMS);
      rbf_f         *rbf=NULL;

//    rbf= new rbf_biharmonic();
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);


      std::cout << std::scientific;


   // allocate arrays
     {
      u   = new double [DIMS];
      v   = new double [DIMS];
      w   = new double [DIMS];
      x   = new double*[npts ];
      x1  = new double*[npts ];
      x2  = new double*[npts ];
      n   = new double*[npts ];
      nt  = new double*[nTest];
      xt0 = new double*[nTest];
      xt1 = new double*[nTest];
      xt2 = new double*[nTest];
      scales=new double[npts];
      scalesTest=new double[npts];

      for( i=0; i<npts;  i++ ){ x[  i]=NULL; x[  i] = new double[DIMS]; }
      for( i=0; i<npts;  i++ ){ x1[ i]=NULL; x1[ i] = new double[DIMS]; }
      for( i=0; i<npts;  i++ ){ x2[ i]=NULL; x2[ i] = new double[DIMS]; }
      for( i=0; i<npts;  i++ ){ n[  i]=NULL; n[  i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ nt[ i]=NULL; nt[ i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt0[i]=NULL; xt0[i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt1[i]=NULL; xt1[i] = new double[DIMS]; }
      for( i=0; i<nTest; i++ ){ xt2[i]=NULL; xt2[i] = new double[DIMS]; }
     }


   // construct points and normals
//    PointsOnSharpBranch(  nCircle, nPlane, nBranch,          rmin, rmax, zmin, zmax, scale, fillet, 0., x, n, scales );
      PointsOnFilletBranch( nCircle, nPlane, nBranch, nFillet, rmin, rmax, zmin, zmax, scale, fillet, 0., x, n, scales );

   // construct test points and normals offset around circle by r*2*pi  so not hitting interpolation points
      r=0.13;
//    PointsOnSharpBranch(  nCircleTest, nPlaneTest, nBranchTest,              rmin, rmax, zmin, zmax, scale, fillet, r*pi2, xt0, nt, scalesTest );
      PointsOnFilletBranch( nCircleTest, nPlaneTest, nBranchTest, nFilletTest, rmin, rmax, zmin, zmax, scale, fillet, r*pi2, xt0, nt, scalesTest );

   // mirror all points over xy plane
      for( i=0; i<nTest2; i++ )
     {
         vec.eq( xt0[i+nTest2], xt0[i] );
         vec.eq( nt[ i+nTest2], nt[ i] );

         xt0[i+nTest2][2]*=-1;
         nt[ i+nTest2][2]*=-1;

         scalesTest[i+nTest2]=scalesTest[i];
     }

      for( i=0; i<npts2;  i++ )
     {
         vec.eq( x[i+npts2], x[i] );
         vec.eq( n[i+npts2], n[i] );

         x[i+npts2][2]*=-1;
         n[i+npts2][2]*=-1;

         scales[i+npts2]=scales[i];
     }

   // displace test points off surface along normal
   // xt1 = off*nt + xt0
      const double off=3.0;
      for( i=0; i<nTest; i++ ){ vec.muladd( xt1[i], off*scalesTest[i], nt[i], xt0[i] ); }
      for( i=0; i<npts;  i++ ){ vec.muladd(  x1[i], off*scales[i],      n[i],   x[i] ); }


   // build implicit surface interpolation
      rbf_surf isurf( DIMS, npts, x, n, scales, accuracy, rbf );
      i = isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // build projections
      for( i=0; i<nTest; i++ ){ isurf.project( xt1[i], xt2[i] ); }
      for( i=0; i<npts;  i++ ){ isurf.project( x1[i],  x2[i]  ); }


   // write points to file
     {
      std::ofstream ofile;
      ofile.open( "surface.dat" );

      std::ofstream ofileTest;
      ofileTest.open( "surfaceTest.dat" );

      for( i=0; i<npts; i++ )
     {
         ofileTest << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << isurf.F( x[i] ) << std::endl;
     }

      for( i=0; i<npts; i++ )
     {
         ofile << x1[i][0] << " " << x1[i][1] << " " << x1[i][2] << " " << isurf.F( x1[i] ) << std::endl;
     }

      ofile.close();
      ofileTest.close();
     }


   // print projection errors
      std::cout << std::endl;
      std::cout << "Test points: exact, not training points" << std::endl;
      std::cout << "|error|, F(x), error vector " << std::endl;
      for( i=0; i<nTest; i++ )
     {
         vec.sub( v, xt0[i], xt2[i] );
         dr=vec.length( v );// / ( off*scalesTest[i] );
         r=isurf.F( xt2[i] );

         std::cout << dr << ", " << r << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }

      std::cout << std::endl;
      std::cout << "Test points: exact, training points" << std::endl;
      std::cout << "|error|, F(x), error vector " << std::endl;
      for( i=0; i<npts; i++ )
     {
         vec.sub( v, x[i], x2[i] );
         dr=vec.length( v );
         r=isurf.F( x2[i] );

//       std::cout << dr << ", " << r << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }


   // test projection of point on implicit surface back onto implicit surface
      for( i=0; i<nTest; i++ )
     {
         vec.eq( xt0[i], xt2[i] );

         isurf.normal( xt0[i], nt[i] );

         vec.muladd( xt1[i], off*scalesTest[i], nt[i], xt0[i] );

         isurf.project( xt1[i], xt2[i] );
     }

      for( i=0; i<npts; i++ )
     {
         vec.eq( x[i], x2[i] );

         isurf.normal( x[i], n[i] );

         vec.muladd( x1[i], off*scales[i], n[i], x[i] );

         isurf.project( x1[i], x2[i] );
     }


      std::cout << std::endl;
      std::cout << "Test points: implicit, not training points" << std::endl;
      std::cout << "|error|, F(x), error vector " << std::endl;
      for( i=0; i<nTest; i++ )
     {
         vec.sub( v, xt0[i], xt2[i] );
         dr=vec.length( v );
         r=isurf.F( xt2[i] );

         std::cout << dr << ", " << r << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }

      std::cout << std::endl;
      std::cout << "Test points: implicit, training points" << std::endl;
      std::cout << "|error|, F(x), error vector " << std::endl;
      for( i=0; i<npts; i++ )
     {
         vec.sub( v, x[i], x2[i] );
         dr=vec.length( v );
         r=isurf.F( x2[i] );

//       std::cout << dr << ", " << r << ", (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
     }

      r=dr;

   // deallocate pointers
     {
      for( i=0; i<npts;  i++ ){ delete[] x[  i]; x[  i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] x1[ i]; x1[ i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] x2[ i]; x2[ i]=NULL; }
      for( i=0; i<npts;  i++ ){ delete[] n[  i]; n[  i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] nt[ i]; nt[ i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt0[i]; xt0[i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt1[i]; xt1[i]=NULL; }
      for( i=0; i<nTest; i++ ){ delete[] xt2[i]; xt2[i]=NULL; }

      delete[]  u;    u=NULL;
      delete[]  v;    v=NULL;
      delete[]  w;    w=NULL;
      delete[]  x;    x=NULL;
      delete[] x1;   x1=NULL;
      delete[] x2;   x2=NULL;
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

