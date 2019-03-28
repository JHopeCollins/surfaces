# include <rbf_interp.h>
# include <vec.h>

/*
 * structs for interpolating with radial basis functions.
 * adapted from code given in Numerical Recipes (C++) Third Edition section 3.7
 */


   rbf_interp<int DIM_, typename POS_, typename VAL_>::rbf_interp( int nn, POS_ *pts, VAL_ *vals )
  {
      init( d, nn, pts, vals );
  }

   void rbf_interp<int DIM_, typename POS_, typename VAL_>::init( int nn, POS_ *pts, VAL_ *vals )
  {
      int i,j;

      n=nn;

      malloc();

      for( i=0; i<n; i++ )
     {
         val[i] = vals[i];
         for( j=0; j<DIMS_; j++ ){ pt[i][j]=pts[i][j]; }
     }
  }

   void rbf_interp<int DIM_, typename POS_, typename VAL_>::malloc()
  {
      w   = new VAL_ [n];
      val = new VAL_ [n];
      pt  = new POS_ [n];
      arrays=1;
  }

   void rbf_interp<int DIM_, typename POS_, typename VAL_>::free()
  {
      delete[] pt;
      delete[] val;
      delete[] w;
      arrays=0;
  }

// solve for the centre weights given the radial basis function definition
   int rbf_interp<int DIM_, typename POS_, typename VAL_>::build_weights()
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

