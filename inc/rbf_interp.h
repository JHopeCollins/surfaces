# ifndef RBFINTERP_H
# define RBFINTERP_H

# include <iostream>
# include <math.h>
# include <assert.h>

# include <rbf.h>
# include <vec.h>
# include <surf.h>

# define  EPS  1.0e-8

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
      rbf_interp() { pt=NULL; val=NULL; w=NULL; rbf=NULL; }

      rbf_interp( int d, int nn, double **pts, double *vals, rbf_f *func );

   // destructor - deallocate arrays
      virtual ~rbf_interp(){ if( arrays ){ free(); } }

   // initialisation of struct
      void init(  int d, int nn, double **pts, double *vals, rbf_f *func );

   // memory management for arrays
      void malloc();
      void free();

   // solve for the centre weights given the radial basis function definition
      int build_weights();

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
               l = indx( k,k, dims );
               h[l] += wr*df;

               val1 = val0*xi[k];

               for( j=0; j<dims; j++ )
              {
                  val = val1*xi[j];

                  l = indx( k,j, dims );
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


# endif
