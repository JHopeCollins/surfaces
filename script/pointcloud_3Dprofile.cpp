#include <cstdio>
#include <cmath>

#include <ssurf.h>

   int main()
  {
      double       height=1.0;   // blade height
      double       dz;

      double      *scales;
      vtx_t       *normal; // normal/tangent vectors at point on blade
      vtx_t       *points;             // points on blade
      vtx_t        v,m;

      int         i,j;
      int         nspan=10;           // number of points across blade span
      int         nprofile;
      int         npoints;


      FILE *fr=fopen( "data/pointclouds/profile2D", "r" );

      fscanf( fr, "%d\n", &i );         // dimension of pointcloud space
      fscanf( fr, "%d\n", &nprofile );   // number of points in cloud


      points=new vtx_t [nprofile];
      normal=new vtx_t [nprofile];
      scales=new double[nprofile];

   // read 2D profile
      for( i=0; i<nprofile; i++ )
     {
         normal[i][2]=0.;
         points[i][2]=0.;
         fscanf( fr, "%lf %lf %lf %lf %lf\n", &points[i][0], &points[i][1], &normal[i][0], &normal[i][1], &scales[i] );
     }
      fclose( fr );

      npoints=(nspan+6)*nprofile;

      FILE *fw=fopen( "data/pointclouds/profile3D", "w" );
      fprintf( fw, "%d\n", 3 );
      fprintf( fw, "%d\n", npoints );

      dz= height    /(nspan-1);
      dz=(height-dz)/(nspan-1);

   // write 3D point cloud
      m[0]=0.;
      m[1]=0.;
      for( i=0; i<nprofile; i++ )
     {
         m[0]=0.;
         m[1]=0.;

      // end caps
         v=points[i] - 1.5*scales[i]*normal[i];

         v[2]= 0.;
         m[2]=-1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

         v[2]=height;
         m[2]= 1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

         v=0.9*points[i];

         v[2]= 0.;
         m[2]=-1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

         v[2]=height;
         m[2]= 1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

         v=0.6*points[i];

         v[2]= 0.;
         m[2]=-1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

         v[2]=height;
         m[2]= 1.;
         fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );

      // spanwise points
         for( j=0; j<nspan; j++ )
        {
            v=points[i];
            m=normal[i];
            v[2]=(j+0.5)*dz;
            fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], m[0], m[1], m[2], scales[i] );
        }
     }

      fclose( fw );

      delete[] points;
      delete[] normal;
      delete[] scales;

      return 0;
  }
