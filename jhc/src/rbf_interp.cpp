# include <rbf_interp.h>
# include <vec.h>

/*
 * structs for interpolating with radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */


   rbf_interp::rbf_interp( int d, int nn, double **pts, double *vals, rbf_f *func )
  {
      init( d, nn, pts, vals, func );
  }

   void rbf_interp::init(  int d, int nn, double **pts, double *vals, rbf_f *func )
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

   void rbf_interp::malloc()
  {
      w   = new double [n];
      val = new double [n];
      pt  = new double*[n];

      for( int i=0; i<n; i++ ){ pt[i] = new double[dims]; }
      arrays=1;
  }

   void rbf_interp::free()
  {
      for( int i=0; i<n; i++ ){ delete[] pt[i]; }
      delete[] pt;
      delete[] val;
      delete[] w;
      arrays=0;
  }

// solve for the centre weights given the radial basis function definition
   int rbf_interp::build_weights()
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
            r = radius( pt[i], pt[j] );
            v = (*rbf)( r );
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

