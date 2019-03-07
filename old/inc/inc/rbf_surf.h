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
   struct rbf_surf : public rbf_interp, surfimplicit_t
  {
      int m;              // number of on-surface points
      double     scale;   // distance to displace artificial points off surface
      double    *scales;  // distance to displace artificial points off surface
      double    **norm;   // surface normals at on-surface points
      bool      scalearray=false;

      rbf_surf(){ norm=NULL; scales=NULL; }
     ~rbf_surf(){ if( arrays ){ free(); } }

   // constructor - creates artificial points off the surface for distance function before interpolating
      rbf_surf( int d, int nn, double **p, double **norms, double  sc, double acc, rbf_f *func );
      rbf_surf( int d, int nn, double **p, double **norms, double *sc, double acc, rbf_f *func );

   // memory management for arrays
      void malloc();
      void free();

   // setters
      void set_scalars( int d, int nn, double acc );

      void set_RadialBasisFunction( rbf_f *func );

      void set_scale( double  sc );

      void set_scale( double *sc );

      void set_PointsNormals( double **p, double **norms );

      void fread( char *name );

   // construct on surface points with zero distance function to interpolate
      void on_surface_points( double **p, double **norms );

   // construct off surface points with non-zero distance function to interpolate
      void off_surface_points();

   // interpolated distance pseudo-function of point p0 from surface
      inline double F( double *p0 ){ return interpolate( p0 ); }

   // first gradient of interpolated pseudo-distance function at p0
      inline void dF( double *p0, double *d ){ derivative( p0, d ); }

   // second gradient (hessian) of interpolated pseudo-distance function at p0
      inline void d2F( double *p0, double *dd ){ hessian( p0, dd ); }

   // q is projection of point p onto the surface to within err % of surface scale, calculated by walking down distance function gradient
      void project_bydistanceonly( double *p, double *q );
  };

# endif
