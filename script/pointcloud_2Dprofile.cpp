#include <cstdio>

#include <constants.h>
#include <ssurf.h>

   int main()
  {
      stat_t        stat;
      ssurf_t      blade;
      tran_t   transform;

      int                   i;
      double       height=1.0;   // blade height
      double        theta=0.0;   // rotation angle
      double     z[2]={0.,0.};   // parametric coordinates on blade
      double       scale=0.01;   // scale for point cloud
      double         distance;   // distance between two points

      vtx_t       tangent[2]; // normal/tangent vectors at point on blade
      vtx_t      *points=NULL, *normal=NULL;             // points on blade
      vtx_t       v;

      transform.x0 = vtx_t( 0.,0.,0. ); // origin
      transform.c  = cos( theta );
      transform.s  = sin( theta );
      transform.d  = 0.1;               // size scaling

      blade.read( (char*)"data/profile", transform, height );
      int npoints=blade.y.n/2;

      points = new vtx_t[npoints];
      normal = new vtx_t[npoints];

   // first point
      int n=0;
      blade.nearest( blade.y[0], z[0], z[1] );
      blade.pos( z[0], z[1], points[n], tangent );

      normal[n] = cross( tangent[0], tangent[1] );
      normal[n]/= -length( normal[n] );
      n++;

      for( i=1; i<npoints; i++ )
     {
         blade.nearest( blade.y[i], z[0], z[1] );
         blade.pos( z[0], z[1], points[n], tangent );

         distance=length( points[n]-points[n-1] );

         if( distance > EPS )
        {
            normal[n] = cross( tangent[0], tangent[1] );
            normal[n]/= -length( normal[n] );

            n++;
        }
     }

      FILE *f=fopen( "data/pointclouds/profile2D", "w" );
      fprintf( f, "%d\n", 2 );         // dimension of pointcloud space
      fprintf( f, "%d\n", n );   // number of points in cloud
      for( i=0; i<n; i++ )
     {
         fprintf( f, "%lf %lf %lf %lf %lf\n", points[i][0], points[i][1], normal[i][0], normal[i][1], scale );
     }
      fclose( f );

      delete[] points; points=NULL;

      return 0;
  }
