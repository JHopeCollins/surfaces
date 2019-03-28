# ifndef RBFINTERP_H
# define RBFINTERP_H

# include <cstdio>
# include <iostream>
# include <math.h>
# include <assert.h>

# include <rbf.h>
//# include <vec.h>
# include <surf.h>

# define  EPS  1.0e-8

/*
 * structs for interpolating with radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */

#  define RBF rbf_triharmonic

// struct for interpolating n points in dim dimensions using radial basis function rbf
   template <int DIM_, typename POS_, typename VAL_>
   struct rbf_interp
  {
      int           n;

      POS_        *pt;
      VAL_       *val;
      VAL_         *w;
      RBF         rbf;

   // constructor
      rbf_interp() { pt=NULL; val=NULL; w=NULL; }

      rbf_interp( int nn, POS_ *pts, VAL_ *vals )
  {
      init( nn, pts, vals );
  }

   // destructor - deallocate arrays
      virtual ~rbf_interp(){ free(); }

   // initialisation of struct
      void init( int nn, POS_ *pts, VAL_ *vals )
  {
      int i,j;

      n=nn;

      malloc();

      for( i=0; i<n; i++ )
     {
         val[i] = vals[i];
         for( j=0; j<DIM_; j++ ){ pt[i][j]=pts[i][j]; }
     }
  }

   // memory management for arrays
      void malloc()
  {
      w   = new VAL_ [n];
      val = new VAL_ [n];
      pt  = new POS_ [n];
  }

      void free()
  {
      delete[] pt;
      delete[] val;
      delete[] w;
  }

   // solve for the centre weights given the radial basis function definition
      int build_weights()
  {
      double      *rbf_mat;    // y( |r_i - r_j| )
      double      *rhs;
      double      v,r;
      int         i,j;
      POS_          q;

      int  *ipiv;
      int  info;

      rbf_mat=new double[n*n];
      rhs    =new double[n];
      ipiv   =new    int[n];

   // diagonals are y( |r_i - r_i| ) = y( 0 )
      v = rbf( 0. );
      for( i=0; i<n; i++ )
     {
         j=indx(i,i,n);

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
            q = pt[i]-pt[j];
            r = length( q );
            
            v = rbf( r );
            rbf_mat[ indx(i,j,n) ] = v;
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
      inline VAL_ interpolate( POS_ p0 )
     {
         double r=0;
         VAL_   v(0.);

         for( int i=0; i<n; i++ )
        {
            r = length( p0-pt[i] );
            v+= w[i]*rbf( r );
        }
         return v;
     }

      inline void interpolate(POS_ p0, VAL_ &y ){ y=interpolate( p0 ); }

   // first derivative of the interpolated value at p0
      inline void derivative( POS_ p0, VAL_ *dy )
     {
         int i,j;

         VAL_    s;
         POS_    q;
         REAL_   r;

         for( i=0; i<DIM_; i++ ){ dy[i]=VAL_(0.); }

         for( i=0; i<n; i++ )
        {
            q = p0 - pt[i];
            r = length( q );

            s = w[i]*rbf.d( r )/fmax( r,EPS );

            for( j=0; j<DIM_; j++ ){ dy[j]+= s*q[j]; }
        }

         return;
     }

   // second derivative of the interpolated value at p0

      inline void hessian( POS_ p0, VAL_ *d2y )
     {
         int      i=0,j=0,k=0,l=0;

         POS_     x0;
         POS_     xi;

         REAL_    r,r1;
         VAL_     wr, df, d2f, val0, val1, val;

         for( k=0; k<DIM_*DIM_; k++ ){ d2y[k]=VAL_(0.); }

      // build hessian
         for( i=0; i<n; i++ )
        {
            xi= p0- pt[i];
            r = length(  xi );
            r1= 1./fmax( r,EPS );
            wr  = w[i]*r1;
            df  = rbf.d(  r );
            d2f = rbf.d2( r );

            val0 = wr*( d2f - df*r1 )*r1;

            for( k=0; k<DIM_; k++ )
           {
               l = indx( k,k, DIM_ );
               d2y[l] += wr*df;

               val1 = val0*xi[k];

               for( j=0; j<DIM_; j++ )
              {
                  val = val1*xi[j];

                  l = indx( k,j, DIM_ );
                  d2y[l] += val;
              }
           }
        }

     }

      inline VAL_ operator()(        POS_  p0 ){   return interpolate(    p0 ); }
      inline void operator()( int n, POS_ *p0, VAL_ *y ){ interpolate( n, p0, y ); }

      inline void interpolate( int np, POS_ *p0, VAL_  *y )
     {
         for( int i=0; i<np; i++ ){ y[i]=interpolate( p0[i] ); }
     }
      inline void derivative(  int np, POS_ *p0, VAL_ **dy )
     {
         for( int i=0; i<np; i++ ){ derivative( p0[i], dy[i] ); }
     }
      inline void hessian(     int np, POS_ *p0, VAL_ **d2y )
     {
         for( int i=0; i<np; i++ ){ hessian( p0[i], d2y[i] ); }
     }

      virtual void fread( FILE *f )
     {
       ::fread( &n, 1, sizeof(int ), f );
         malloc();
         delete[] val; val=nullptr;
       ::fread( pt, n, sizeof(POS_), f );
       ::fread(  w, n, sizeof(VAL_), f );
         return;
     }

      virtual void fwrite( FILE *f )
     {
       ::fwrite( &n, 1, sizeof(int ), f );
       ::fwrite( pt, n, sizeof(POS_), f );
       ::fwrite(  w, n, sizeof(VAL_), f );
         return;
     }
  };

# endif
