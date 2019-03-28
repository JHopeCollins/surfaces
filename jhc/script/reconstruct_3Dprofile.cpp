#include <iostream>

#include <constants.h>
#include <rbf.h>
#include <rbf_surf.h>

#define DIMS 3


   int main()
  {
      int           i,j;

//    double      xmin=-3.5, xmax=0.;
//    double      ymin=-2.5, ymax=1.5;
//    double      zmin= 0.0, zmax=1.0;

//    double      dx,dy,dz, d;
//    double      x[DIMS];

//    int         nplot=51;

      vec_t vec(DIMS);
      rbf_f *rbf=NULL;

//    rbf= new rbf_biharmonic();
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);

      rbf_surf isurf;


   // build surface
      isurf.accuracy=EPS;
      isurf.set_RadialBasisFunction( rbf );
      isurf.fread( (char*)"data/pointclouds/profile3D" );

      i=isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // write surface points
      FILE *f=fopen( "data/pointclouds/profile3Drbfpoints", "w" );
      for( i=0; i<isurf.n; i++ )
     {
         fprintf( f, "%lf %lf %lf\n", isurf.pt[i][0], isurf.pt[i][1], isurf.pt[i][2] );
     }
      fclose( f );

/*
      f=fopen( "data/pointclouds/profile3Ddistance", "w" );
      dx=(xmax-xmin)/(nplot-1);
      dy=(ymax-ymin)/(nplot-1);
      dz=(zmax-zmin)/(nplot-1);
      for( i=0; i<nplot; i++ )
     {
         x[0]=xmin+i*dx;
         for( j=0; j<nplot; j++ )
        {
            x[1]=ymin+j*dy;
            for( k=0; k<nplot; k++ )
           {
               x[2]=zmin+k*dz;
               d=isurf.F( x );
               fprintf( f, "%lf, %lf %lf %lf\n", x[0], x[1], x[2], d );
           }
        }
     }
      fclose( f );
*/


   // test projections
      double **xt=NULL;
      int ntest=isurf.m;
      xt=new double*[ntest];
      for( i=0; i<ntest; i++ ){ xt[i]=NULL; xt[i]=new double[DIMS]; }

      double r=1.5;
      double   *yt=new double[DIMS];
   // generate points outside blade profile and project points onto blade profile
      for( i=0; i<ntest; i++ )
     {
         j=3*i;
         yt[0]=isurf.pt[j][0] + r*isurf.scales[i]*isurf.norm[i][0];
         yt[1]=isurf.pt[j][1] + r*isurf.scales[i]*isurf.norm[i][1];
         yt[2]=isurf.pt[j][2] + r*isurf.scales[i]*isurf.norm[i][2];

         xt[i][0]=yt[0];
         xt[i][1]=yt[1];
         xt[i][2]=yt[2];

         //std::cout << "(" << xt[i][0] << ", " << xt[i][1] << ", " << xt[i][2] << ") ";

         isurf.project( xt[i], xt[i] );

         //std::cout << "(" <<       xt[i][0] << ", " <<       xt[i][1] <<  ", " <<       xt[i][2] << ") " << std::endl;
         //std::cout << "(" << isurf.pt[j][0] << ", " << isurf.pt[j][1] <<  ", " << isurf.pt[j][2] << ") " << vec.radius( xt[i], yt )/isurf.scales[i] << std::endl;
     }


   // write projected points
      f=fopen( "data/pointclouds/profile3Dproject", "w" );
      for( i=0; i<ntest; i++ ){ fprintf( f, "%lf %lf %lf\n", xt[i][0], xt[i][1], xt[i][2] ); }
      fclose( f );


      for( i=0;i<ntest; i++ ){ delete[] xt[i]; xt[i]=NULL; }
      delete[] xt;  xt=NULL;
      delete[] yt;  yt=NULL;
      delete  rbf; rbf=NULL;

      return 0;
  }
