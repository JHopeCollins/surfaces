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

      rbf_surf<3, vtx_t, double> isurf0, isurf;


   // build surface
      isurf0.accuracy=EPS;
      isurf0.freadPointCloud( (char*)"data/pointclouds/profile3D" );

      i=isurf0.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;

      FILE *f=fopen( "data/profile3Dsurf.bin", "w" );
      isurf0.fwrite( f );
      fclose( f );

      f=fopen( "data/profile3Dsurf.bin", "r" );
      isurf.fread( f );
      fclose( f );

   // write surface points
      f=fopen( "data/pointclouds/profile3Drbfpoints", "w" );
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
      int ntest=isurf0.m;
      vtx_t *xt=new vtx_t[ntest], yt;
      double r=1.5;
   // generate points outside blade profile and project points onto blade profile
      for( i=0; i<ntest; i++ )
     {
         j=3*i;

         yt = isurf.pt[j] + r*isurf0.scales[i]*isurf0.norm[i];
         xt[i]=yt;

         //std::cout << "(" << xt[i][0] << ", " << xt[i][1] << ", " << xt[i][2] << ") ";

         isurf.project( xt[i], xt[i] );

         //std::cout << "(" <<       xt[i][0] << ", " <<       xt[i][1] <<  ", " <<       xt[i][2] << ") " << std::endl;
         //std::cout << "(" << isurf.pt[j][0] << ", " << isurf.pt[j][1] <<  ", " << isurf.pt[j][2] << ") " << length( xt[i] - yt )/isurf0.scales[i] << std::endl;
     }


   // write projected points
      f=fopen( "data/pointclouds/profile3Dproject", "w" );
      for( i=0; i<ntest; i++ ){ fprintf( f, "%lf %lf %lf\n", xt[i][0], xt[i][1], xt[i][2] ); }
      fclose( f );


      delete[] xt;  xt=NULL;

      return 0;
  }
