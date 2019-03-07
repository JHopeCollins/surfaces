#include <cstdio>
#include <cmath>

#include <ssurf.h>

   int main()
  {
      double       height=1.0;   // blade height
      double       scale=0.01;   // scale for point cloud
      double       dz;

      vtx_t       normal; // normal/tangent vectors at point on blade
      vtx_t       p;             // points on blade
      vtx_t       te;

      int         i,j;
      int         nspan=10;           // number of points across blade span
      int         nprofile;
      int         npoints;


      FILE *fr=fopen( "data/pointclouds/profile2Dunique", "r" );
      FILE *fw=fopen( "data/pointclouds/profile3Dunique", "w" );

      fscanf( fr, "%d\n", &i );         // dimension of pointcloud space
      fscanf( fr, "%d\n", &nprofile );   // number of points in cloud

      npoints=nspan*nprofile;

      fprintf( fw, "%d\n", 3 );
      fprintf( fw, "%d\n", npoints );

      dz=height/(nprofile-1);
      normal[2]=0.;
      for( i=0; i<nprofile; i++ )
     {
         fscanf( fr, "%lf %lf %lf %lf %lf\n", &p[0], &p[1], &normal[0], &normal[1], &scale );

         for( j=0; j<nspan; j++ )
        {
            p[2]=j*dz;

            fprintf( fw, "%lf %lf %lf %lf %lf %lf %lf\n", p[0], p[1], p[2], normal[0], normal[1], normal[2], scale );
        }
     }

      fclose( fr );
      fclose( fw );

      return 0;
  }
