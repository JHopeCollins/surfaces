#include <cstdio>
#include <cmath>

#include <ssurf.h>

   int main()
  {
      stat_t        stat;
      ssurf_t      blade;
      tran_t   transform;

      double       height=1.0;   // blade height
      double        theta=0.0;   // rotation angle
      double     z[2]={0.,0.};   // parametric coordinates on blade
      double       scale=0.01;   // scale for point cloud

      vtx_t       normal, tangent[2]; // normal/tangent vectors at point on blade
      vtx_t      *p=NULL;             // points on blade
      vtx_t       te;

      transform.x0 = vtx_t( 0.,0.,0. ); // origin
      transform.c  = cos( theta );
      transform.s  = sin( theta );
      transform.d  = 0.1;               // size scaling

      blade.read( (char*)"data/profile", transform, height );
      int npoints=blade.y.n/2;

      FILE *f=fopen( "data/pointclouds/profile2D", "w" );
      fprintf( f, "%d\n", 2 );         // dimension of pointcloud space
      fprintf( f, "%d\n", npoints );   // number of points in cloud

      p = new vtx_t[npoints];

      for( int i=0; i<npoints; i++ )
     {
         blade.nearest( blade.y[i], z[0], z[1] );
         blade.pos( z[0], z[1], p[i], tangent );

         normal = cross( tangent[0], tangent[1] );
         normal/= -length( normal );

         fprintf( f, "%lf %lf %lf %lf %lf\n", p[i][0], p[i][1], normal[0], normal[1], scale );
     }

/*
   // construct trailing edge point
      te = 0.5*( p[0] + p[npoints-1] );
      normal[0]= 1.;
      normal[1]=-1.;
      normal[2]=-0.;
      normal/= length( normal );


      fprintf( f, "%lf %lf %lf %lf %lf\n", te[0], te[1], normal[0], normal[1], scale );
*/

      fclose( f );

      delete[] p; p=NULL;

      return 0;
  }
