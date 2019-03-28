# include <iostream>

# include <constants.h>
# include <rbf.h>
# include <rbf_surf.h>

# define DIMS 2

   int main()
  {
      int           i,j;

      double      xmin=-1.5, xmax=1.25;
      double      ymin=-2.25, ymax=1.0;
      double      dx,dy, z;
      vtx_t       x;

      int         nplot=81;
      int         ntest=36;

      rbf_surf<2, vtx_t, double> isurf;


   // build surface
      isurf.accuracy=EPS;
      isurf.fread( (char*)"data/pointclouds/profile2D" );

      i=isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;


   // write surface points
      FILE *f=fopen( "data/pointclouds/profile2Drbfpoints", "w" );
      for( i=0; i<isurf.n; i++ ){ fprintf( f, "%lf, %lf\n", isurf.pt[i][0], isurf.pt[i][1] ); }
      fclose( f );

      f=fopen( "data/pointclouds/profile2Ddistance", "w" );
      dx=(xmax-xmin)/(nplot-1);
      dy=(ymax-ymin)/(nplot-1);
      x[2]=0.;
      for( i=0; i<nplot; i++ )
     {
         x[0]=xmin+i*dx;
         for( j=0; j<nplot; j++ )
        {
            x[1]=ymin+j*dy;
            z=isurf.F( x );
            fprintf( f, "%lf, %lf %lf\n", x[0], x[1], z );
        }
     }
      fclose( f );


   // test projections
      ntest=isurf.m;
      vtx_t *xt=new vtx_t [ntest], yt;
      double r=1.5;
   // generate points outside blade profile and project points onto blade profile
      for( i=0; i<ntest; i++ )
     {
         j=3*i;

         yt    = isurf.pt[j] + r*isurf.scales[i]*isurf.norm[i];
         xt[i] = yt;

         std::cout << "(" << xt[i][0] << ", " << xt[i][1] << ") ";

         isurf.project( xt[i], xt[i] );

         std::cout << "(" <<       xt[i][0] << ", " <<       xt[i][1] << ") " << std::endl;
         std::cout << "(" << isurf.pt[j][0] << ", " << isurf.pt[j][1] << ") " << length( xt[i] - yt )/isurf.scales[i] << std::endl;
     }


   // write projected points
      f=fopen( "data/pointclouds/profile2Dproject", "w" );
      for( i=0; i<ntest; i++ ){ fprintf( f, "%lf, %lf\n", xt[i][0], xt[i][1] ); }
      fclose( f );


      delete[] xt;  xt=NULL;

      return 0;
  }
