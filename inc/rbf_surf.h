# ifndef RBFSURF_H
# define RBFSURF_H

# include <iostream>
# include <math.h>
# include <assert.h>

# include <rbf.h>
# include <vec.h>
# include <surf.h>
# include <rbf_interp.h>

# define  EPS  1.0e-8

// struct for creating an implicit representation of a surface from a point/normal cloud using an rbf interpolation
   template <int DIM_, typename POS_, typename VAL_>
   struct rbf_surf : public surfimplicit_t<DIM_, POS_, VAL_>,
                                rbf_interp<DIM_, POS_, VAL_>
  {
      int m;               // number of on-surface points
      POS_        *norm;   // surface normals at on-surface points
      REAL_       scale;   // distance to displace artificial points off surface
      REAL_     *scales;   // distance to displace artificial points off surface

      rbf_surf(){ norm=NULL; scales=NULL; }
     ~rbf_surf(){ free(); }

   // constructor - creates artificial points off the surface for distance function before interpolating
      rbf_surf( int nn, POS_ *p, POS_ *norms, REAL_ *sc, REAL_ acc )
  {
      this->n=3*nn;
      m=  nn;
      this->accuracy=acc;

      malloc();

      set_scale( sc );
      this->accuracy*=scale;

      set_PointsNormals( p, norms );
  }

   // memory management for arrays
      void malloc();
      void free();

   // setters
      void set_scale( REAL_ *sc );

      void set_PointsNormals( POS_ *p, POS_ *norms );

      void freadPointCloud( char *name );

   // construct on surface points with zero distance function to interpolate
      void on_surface_points( POS_ *p, POS_ *norms );

   // construct off surface points with non-zero distance function to interpolate
      void off_surface_points();

   // interpolated distance pseudo-function of point p0 from surface
      inline VAL_ F( POS_ p0 ){ return this->interpolate( p0 ); }

   // first gradient of interpolated pseudo-distance function at p0
      inline void dF( POS_ p0, VAL_ *d ){ this->derivative( p0, d ); }

   // second gradient (hessian) of interpolated pseudo-distance function at p0
      inline void d2F( POS_ p0, VAL_ *dd ){ this->hessian( p0, dd ); }

   // q is projection of point p onto the surface to within err % of surface scale, calculated by walking down distance function gradient
      void project_bydistanceonly( POS_ p, POS_ &q );

      void fread ( FILE *f );
      void fwrite( FILE *f );
  };


// constructor - creates artificial points off the surface for distance function before interpolating
/*
   template <int DIM_, typename POS_, typename VAL_>
   rbf_surf<DIM_, POS_, VAL_>::rbf_surf( int nn, POS_ *p, POS_ *norms, REAL_ *sc, REAL_ acc  )
  {
      n=3*nn;
      m=  nn;
      accuracy=acc;

      malloc();

      set_scale( sc );
      accuracy*=scale;

      set_PointsNormals( p, norms );
  }
*/

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::malloc()
  {
      rbf_interp<DIM_, POS_, VAL_>::malloc();
      scales=new double[m];
      norm = new POS_  [m];
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::free()
  {
      delete[] scales; scales=NULL;
      delete[] norm;   norm = NULL;
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::set_scale( REAL_ *sc )
  {
      scale=0;
      for( int i=0; i<m; i++ )
     {
         scales[i]=sc[i];
         scale   +=sc[i];
     }
      scale/=m;
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::set_PointsNormals( POS_ *p, POS_ *norms )
  {
      on_surface_points( p, norms );
      off_surface_points();
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::freadPointCloud( char *name )
  {
      int i,j,k;

      FILE *f=fopen( name, "r" );
      assert( f );

   // read dimensions and number of points in cloud
      int dims;
      fscanf( f, "%d", &dims );
      assert( dims==DIM_ );

      fscanf( f, "%d", &m );
      this->n=3*m;

      malloc();

   // read points, normals and scales
      scale=0.;
      for( i=0; i<m; i++ )
     {
         k=3*i;
         for( j=0; j<DIM_; j++ ){ fscanf( f, "%lf", &(this->pt)[k][j] ); }
         for( j=0; j<DIM_; j++ ){ fscanf( f, "%lf", &norm[i][j] ); }
         fscanf( f, "%lf", &scales[i] );
         scale+=scales[i];
     }
      fclose( f );
   
      scale/=m;
      this->accuracy*=scale;

      off_surface_points();
  }

// construct on surface points with zero distance function to interpolate
   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::on_surface_points( POS_ *p, POS_ *norms )
  {
      int i,k;
      REAL_ l;

   // points on surface
      for( i=0; i<m; i++ )
     {
         k=3*i;
         this->val[k]=0.0;

         this->pt[k]=p[i];

         norm[i]=norms[i];
         l=length( norm[i] );
         norm[i]/=l;
     }
  }

// construct off surface points with non-zero distance function to interpolate
   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::off_surface_points()
  {
      int   i,k;

   // distance function increases along surface normal
      for( i=0; i<m; i++ )
     {
         k=3*i;
         this->val[k+1] = VAL_( scales[i] );
         this->val[k+2] = VAL_(-scales[i] );

         this->pt[k+1] = this->pt[k] + scales[i]*norm[i];
         this->pt[k+2] = this->pt[k] - scales[i]*norm[i];
     }
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::fwrite( FILE *f )
  {
      rbf_interp<DIM_,POS_,VAL_>::fwrite( f );
    ::fwrite( &(this->accuracy), 1, sizeof(double), f );
  }

   template <int DIM_, typename POS_, typename VAL_>
   void rbf_surf<DIM_, POS_, VAL_>::fread( FILE *f )
  {
      rbf_interp<DIM_,POS_,VAL_>::fread( f );
    ::fread( &(this->accuracy), 1, sizeof(double), f );
  }


# endif
