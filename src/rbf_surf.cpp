# include <cstdio>

# include <rbf_surf.h>
# include <vec.h>

// constructor - creates artificial points off the surface for distance function before interpolating
   rbf_surf::rbf_surf( int d, int nn, double **p, double **norms, double sc,  double acc, rbf_f *func )
  {
      scalearray=false;
      set_scalars( d, nn, acc );
      set_RadialBasisFunction( func );

      malloc();

      set_scale( sc );
      accuracy*=scale;

      set_PointsNormals( p, norms );
  }

   rbf_surf::rbf_surf( int d, int nn, double **p, double **norms, double *sc, double acc, rbf_f *func )
  {
      scalearray=true;
      set_scalars( d, nn, acc );
      set_RadialBasisFunction( func );

      malloc();

      set_scale( sc );
      accuracy*=scale;

      set_PointsNormals( p, norms );
  }

   void rbf_surf::malloc()
  {
      rbf_interp::malloc();
      if( scalearray ){ scales=new double[m]; }
      norm=new double*[m];
      for( int i=0; i<m; i++ ){ norm[i]=NULL; norm[i]=new double[dims]; }
  }

   void rbf_surf::free()
  {
      rbf_interp::free();
      if( scalearray ){ delete[] scales; scales=NULL; }
      for( int i=0; i<m; i++ ){ delete[] norm[i]; norm[i]=NULL; }
      delete[] norm; norm=NULL;
  }

   void rbf_surf::set_scalars( int d, int nn, double acc ){ dims=d; n=3*nn; m=nn; accuracy=acc; }

   void rbf_surf::set_RadialBasisFunction( rbf_f *func ){ rbf=func; }

   void rbf_surf::set_scale( double  sc ){ scale=sc; }

   void rbf_surf::set_scale( double *sc )
  {
      scale=0;
      for( int i=0; i<m; i++ )
     {
         scales[i]=sc[i];
         scale   +=sc[i];
     }
      scale/=m;
  }

   void rbf_surf::set_PointsNormals( double **p, double **norms )
  {
      on_surface_points( p, norms );
      off_surface_points();
  }

   void rbf_surf::fread( char *name )
  {
      scalearray=true;

      int i,j,k;

      FILE *f=fopen( name, "r" );
      assert( f );

   // read dimensions and number of points in cloud
      fscanf( f, "%d", &dims );
      fscanf( f, "%d", &m );
      n=3*m;

      malloc();

   // read points, normals and scales
      scale=0.;
      for( i=0; i<m; i++ )
     {
         k=3*i;
         for( j=0; j<dims; j++ ){ fscanf( f, "%lf", &pt[   k][j] ); }
         for( j=0; j<dims; j++ ){ fscanf( f, "%lf", &norm[ i][j] ); }
         fscanf( f, "%lf", &scales[i] );
         scale+=scales[i];
     }
      fclose( f );
   
      scale/=m;
      accuracy*=scale;

      off_surface_points();
  }

// construct on surface points with zero distance function to interpolate
   void rbf_surf::on_surface_points( double **p, double **norms )
  {
      int i,k;

   // points on surface
      for( i=0; i<m; i++ )
     {
         k=3*i;
         val[k]=0.0;

         eq(     pt[k],     p[i] );
         eq(   norm[i], norms[i] );
         unit( norm[i] );
     }
  }

// construct off surface points with non-zero distance function to interpolate
   void rbf_surf::off_surface_points()
  {
      int   i,j,k;

   // distance function increases along surface normal
      if( !scalearray )
     {
         for( i=0; i<m; i++ )
        {
            k=3*i;
            val[k+1] =  scale;
            val[k+2] = -scale;

            for( j=0; j<dims; j++ )
           {
               pt[k+1][j] = pt[k][j] + scale*norm[i][j];
               pt[k+2][j] = pt[k][j] - scale*norm[i][j];
           }
        }
     }
      else
     {
         for( i=0; i<m; i++ )
        {
            k=3*i;
            val[k+1] =  scales[i];
            val[k+2] = -scales[i];

            for( j=0; j<dims; j++ )
           {
               pt[k+1][j] = pt[k][j] + scales[i]*norm[i][j];
               pt[k+2][j] = pt[k][j] - scales[i]*norm[i][j];
           }
        }
     }
  }

// q is projection of point p onto the surface to within err % of surface scale
   void rbf_surf::project_bydistanceonly( double *p, double *q )
  {
      double *d, *x;
      double dist, err;
      double rlx=0.5;
      int    i;

      d=new double[dims];
      x=new double[dims];

      err = accuracy*scale;

      for( i=0; i<dims; i++ ){ x[i]=p[i]; }

      dist = distance( x );

      while( dist > err )
     {
      // displace point towards surface by distance*relaxation-factor

         derivative( x, d );
         unit(d);

         dist*=rlx;

         for( i=0; i<dims; i++ ){ x[i]=x[i]-dist*d[i]; }

         dist = distance( x );
     }

      eq( q, x );

      delete[] d;
      delete[] x;
  }
