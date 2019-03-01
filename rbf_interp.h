# ifndef RBFINTERP_H
# define RBFINTERP_H

# include <iostream>
# include <math.h>
# include <assert.h>

# include "rbf.h"
# include "vec.h"
# include "surf.h"

# define  EPS  1.0e-10

/*
 * structs for interpolating with radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */

// struct for interpolating n points in dim dimensions using radial basis function rbf
   struct rbf_interp : public virtual vec_t
  {
      int           n;

      int    arrays=0;

      double     **pt;
      double     *val;
      double       *w;

      rbf_f      *rbf;

   // constructor
      rbf_interp( int d, int nn, double **pts, double *vals, rbf_f *func )
     {
         init( d, nn, pts, vals, func );
     }

      void init(  int d, int nn, double **pts, double *vals, rbf_f *func )
     {
         int i,j;

         n=nn;
         dims=d;
         rbf=func;

         malloc();

         for( i=0; i<n; i++ )
        {
            val[i] = vals[i];
            for( j=0; j<dims; j++ ){ pt[i][j]=pts[i][j]; }
        }
     }

      rbf_interp() { pt=NULL; val=NULL; w=NULL; rbf=NULL; }
   // destructor - deallocate arrays
      virtual ~rbf_interp(){ if( arrays ){ free(); } }

      void malloc()
     {
         int i,j;

         w   = new double [n];
         val = new double [n];
         pt  = new double*[n];

         for( i=0; i<n; i++ ){ pt[i] = new double[dims]; }
         arrays=1;
     }

      void free()
     {
         for( int i=0; i<n; i++ ){ delete[] pt[i]; }
         delete[] pt;
         delete[] val;
         delete[] w;
         arrays=0;
     }

   // solve for the centre weights given the radial basis function definition
      int build_weights()
     {
         double      *rbf_mat;    // y( |r_i - r_j| )
         double      *rhs;
         double      v,r;
         int         i,j;

         int  *ipiv;
         int  info;

         rbf_mat=new double[n*n];
         rhs =new      double[n];
         ipiv=new         int[n];

      // diagonals are y( |r_i - r_i| ) = y( 0 )
         v = (*rbf)( 0. );
         for( i=0; i<n; i++ )
        {
            j=index(i,i,n);

            rbf_mat[j] = v;

            rhs[i] = val[i];

            ipiv[i]=0;
        }

      // build matrix in row-major form for fortran based lapack
      // only lower triangular needed for symmetric lapack routine
      // std::cout << "matrix" << std::endl;
         for( i=0; i<n; i++ )
        {
            for(j=0; j<n; j++ )
           {
               r = radius( pt[i], pt[j] );
               v = (*rbf)( r );
               rbf_mat[ index(i,j,n) ] = v;
           }
        }

      // solve for weights
         i=1;
      //          N NRHS   A     LDA  IPIV   B  LDB,  info )
         dgesv_( &n, &i, rbf_mat, &n, ipiv, rhs, &n, &info );

         for( i=0; i<n; i++ ){ w[i] = rhs[i]; }

         delete[] rbf_mat;
         delete[] rhs;
         delete[] ipiv;

         return info;
     }

   // interpolated value at p0
      inline double interpolate( double *p0  )
     {
         double v=0, r=0;
         for( int i=0; i<n; i++ )
        {
            r = radius( p0, pt[i] );
            v+= w[i]* (*rbf)( r );
        }
         return v;
     }

   // first derivative of the interpolated value at p0
      inline void derivative( double *p0, double *dy )
     {
         int i,j;
         double *v, *q, d, r;

         v=new double[dims];
         q=new double[dims];

         eq( v, 0. );

         for( i=0; i<n; i++ )
        {
            sub( q, p0, pt[i] );

            r = length( q );
            if( r < EPS ){ continue; }
            d = w[i]*(*rbf).d( r )/r;

            mul( q, d );
            add( v, q );
        }

         for( j=0; j<dims; j++ ){ dy[j]=v[j]; }

         delete[] v;
         delete[] q;

         return;
     }

   // second derivative of the interpolated value at p0
      inline void hessian( double *p0, double *d2y )
     {
         int      i=0,j=0,k=0,l=0;
         double   *h=NULL;
         double   *x0=NULL, *xi=NULL;
         double   r, wr, df, d2f, val0, val1, val;

         x0=new double [dims];
         xi=new double [dims];
         h =new double [dims*dims];

         eq( x0, p0 );
         for( k=0; k<dims*dims; k++ ){ h[k]=0.; }

//       printf( "i, w, r, df, d2f\n" );

      // build hessian
         for( i=0; i<n; i++ )
        {
            sub( xi, x0, pt[i] );
            r = length(  xi );
            if( r < EPS ){ continue; }
            wr = w[i]/r;
            df  = (*rbf).d(  r );
            d2f = (*rbf).d2( r );

//          printf( "%d, %f, %f, %f, %f\n", i, w[i], r, df, d2f );

            val0 = wr*( d2f - df/r )/r;

            for( k=0; k<dims; k++ )
           {
               l = index( k,k, dims );
               h[l] += wr*df;

               val1 = val0*xi[k];

               for( j=0; j<dims; j++ )
              {
                  val = val1*xi[j];

                  l = index( k,j, dims );
                  h[l] += val;
              }
           }
        }

         for( k=0; k<dims*dims; k++ ){ d2y[k]=h[k]; }

         delete[] x0;
         delete[] xi;
         delete[] h;
     }

      inline double operator()(        double  *p0 ){     return interpolate(    p0 ); }
      inline void   operator()( int n, double **p0, double *y ){ interpolate( n, p0, y ); }

      inline void interpolate( int np, double **p0, double *y )
     {
         for( int i=0; i<np; i++ ){ y[i]=interpolate( p0[i] ); }
     }
      inline void derivative(  int np, double **p0, double **dy )
     {
         for( int i=0; i<np; i++ ){ derivative( p0[i], dy[i] ); }
     }
      inline void hessian(     int np, double **p0, double **d2y )
     {
         for( int i=0; i<np; i++ ){ hessian( p0[i], d2y[i] ); }
     }
  };


// struct for creating an implicit representation of a surface from a point/normal cloud using an rbf interpolation
   struct rbf_surf : public rbf_interp, surfimplicit_t
  {
      int m;              // number of on-surface points
      double     scale;   // distance to displace artificial points off surface
      double    **norm;   // surface normals at on-surface points

      rbf_surf(){ norm=NULL; }
     ~rbf_surf(){ if( arrays ){ free(); } }

   // constructor - creates artificial points off the surface for distance function before interpolating
      rbf_surf( int d, int nn, double **p, double **norms, double sc, double acc, rbf_f *func )
     {
         init( d, nn, p, norms, sc, acc, func );
     }

      void init( int d, int nn, double **p, double **norms, double sc, double acc, rbf_f *func )
     {
         int i,j,k;

         dims=d;
         n=3*nn;
         m=nn;
         rbf=func;
         scale=sc;
         accuracy=acc*scale;

         malloc();

      // points on surface
         for( i=0; i<m; i++ )
        {
            k=3*i;
            val[k]=0.0;

            for( j=0; j<dims; j++ )
           {
               pt[  k][j] =     p[i][j];
               norm[i][j] = norms[i][j];

           }
            unit( norm[i] );
        }

      // create new points with non-zero distance function
         off_surface_points();
     }

      void malloc()
     {
         rbf_interp::malloc();
         norm = new double*[m];
         for( int i=0; i<m; i++ ){ norm[i]=NULL; norm[i]=new double[dims]; }
     }

      void free()
     {
         rbf_interp::free();
         for( int i=0; i<m; i++ ){ delete[] norm[i]; norm[i]=NULL; }
         delete[] norm; norm=NULL;
     }

   // construct off surface points with non-zero distance function to interpolate
      void off_surface_points()
     {
         int   i,j,k;

      // distance function increases along surface normal
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

   // interpolated distance pseudo-function of point p0 from surface
      inline double F( double *p0 ){ return interpolate( p0 ); }

   // first gradient of interpolated pseudo-distance function at p0
      inline void dF( double *p0, double *d ){ derivative( p0, d ); }

   // second gradient (hessian) of interpolated pseudo-distance function at p0
      inline void d2F( double *p0, double *dd ){ hessian( p0, dd ); }

   // q is projection of point p onto the surface to within err % of surface scale
      void project_bydistanceonly( double *p, double *q )
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
  };

# endif
