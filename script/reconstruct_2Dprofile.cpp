# include <iostream>

# include <constants.h>
# include <rbf.h>
# include <rbf_surf.h>

# define DIMS 2

   int main()
  {
      int           i,j;

      double      xmin=-3.5, xmax=0.;
      double      ymin=-2.5, ymax=1.5;
      double      dx,dy, z;
      double      x[2];

      int         nplot=51;

      vec_t vec(DIMS);
      rbf_f *rbf=NULL;

//    rbf= new rbf_biharmonic();
      rbf= new rbf_triharmonic();
//    rbf= new rbf_multiquadratic(scale);
//    rbf= new rbf_invmultiquadratic(scale);
//    rbf= new rbf_thinplate(scale);
//    rbf= new rbf_gaussian(scale);

      rbf_surf isurf;

      isurf.accuracy=EPS;
      isurf.set_RadialBasisFunction( rbf );
      isurf.fread( (char*)"data/pointclouds/profile2D" );

      i=isurf.build_weights();
      std::cout << "lapack info: " << i << std::endl << std::endl;

   // print isurf info
/*
     {
      std::cout << "isurf.dims: " << isurf.dims << std::endl;
      std::cout << "isurf.n: " << isurf.n << std::endl;
      std::cout << "isurf.m: " << isurf.m << std::endl;

      std::cout << "isurf.pt:" << std::endl;
      for( i=0; i<isurf.n; i++ )
     {
         std::cout << isurf.pt[i][0] << ", " << isurf.pt[i][1] << std::endl;
     }
     }
*/

      FILE *f=fopen( "data/pointclouds/profile2Drbfpoints", "w" );
      for( i=0; i<isurf.n; i++ )
     {
         fprintf( f, "%lf, %lf\n", isurf.pt[i][0], isurf.pt[i][1] );
     }
      fclose( f );

      f=fopen( "data/pointclouds/profile2Ddistance", "w" );
      dx=(xmax-xmin)/(nplot-1);
      dy=(ymax-ymin)/(nplot-1);
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

      delete rbf; rbf=NULL;

      return 0;
  }
