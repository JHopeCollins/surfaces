#include <stdio.h>
#include <surf.h>


// class for surface defined by implicit function
   surfimplicit_t::surfimplicit_t( int d, double a )
     {
         dims=d;
         accuracy=a;
     }

// normal vector points along gradient of implicit function local to surface
   void surfimplicit_t::normal( double *p0, double *n )
  {
      double   *x, *d;
      x=new double[dims];
      d=new double[dims];

      project( p0, x );
      dF(      x,  d );
      unit(    d,  n );

      delete[] x;
      delete[] d;
  }

// t are surface tangent vectors at point on surface closest to x
   void surfimplicit_t::tangent( double *p, double *t )   // only for 2D currently
  {
      double *n, *x;

      n=new double[dims];
      x=new double[dims];

      project( p, x );
      normal(  x, n );
      t[0]=  n[1];
      t[1]= -n[0];

      delete[] n;
      delete[] x;
  }

// hessian matrix at point on surface closest to p0
   void surfimplicit_t::curvature( double *p0, double *c )
  {
      double   *x;
      x=new double[dims];

      project( p0, x );
      d2F(      x, c );

      delete[] x;
  }

// project point onto surface, subject to constraints of F(x)=0 and original point orthogonal to projected point
   void surfimplicit_t::project( double *p, double *q )
  {
      int           i=0;
      int           d=dims+1;
      int           nrhs=1;
      int           *ipiv=NULL, info=0;
      vec_t         vecd(d);
      double        *nu=NULL;                   // solution vector: {x,y,z,t}
      double        *x0=NULL;                   // starting point
      double        *df=NULL;                   // derivative of distance field
      double        res=10000000.;             // newton iteration residual
      double        *rhs=NULL, *jac=NULL;            // right hand side and jacobian matrix for newton iterations
      double         rlx=0.3;               // relaxation factor
      double        err=accuracy;    // convergence criteria

      if( F( p ) < err ){ eq(q,p); return; }

      ipiv = new    int[d];
      nu   = new double[d];
      rhs  = new double[d];
      jac  = new double[d*d];
      x0   = new double[dims];
      df   = new double[dims];

   // initialise initial point x0 and solution vector nu

   // tried offsetting initial guess down the distance function gradient slightly from p, didn't seem to make a difference, and often increased number of Newton iterations. Could experiment with it more later.

      eq(  nu, p );

      for( i=0; i<3; i++ )
     {
         dF(  nu, x0  );
         mul( x0, 0.4*F(nu) );
         add( nu, x0 );
     }

      eq(  x0, p );

      nu[dims] = radius( nu, x0 );

      int j=0;
   // start newton iterations
      while( res>err )
     {
//       printf( "projection iteration: %d\n", j );
      // construct rhs and jacobian
         newton_rhs( x0, nu, rhs, df );
         newton_jac( x0, nu, jac, df );

         for( i=0; i<d; i++ ){ ipiv[i]=0; }

      //          N   NRHS   A  LDA  IPIV   B  LDB,  info )
         dgesv_( &d, &nrhs, jac, &d, ipiv, rhs, &d, &info );
         assert( info==0 );

      // update solution vector
         res = vecd.length( rhs );

         vecd.mul( rhs, rlx );
         vecd.add( nu,  rhs );

         rlx*=1.3;
         rlx=std::min(rlx,1.0);

         j++;
     }

   // assign to q
      eq( q, nu );

      double l=radius( p, q );
      printf( "iteration %d, res= %f, dx= %f\n", j, res, l );

      delete[] ipiv;
      delete[] nu;
      delete[] rhs;
      delete[] jac;
      delete[] x0;
      delete[] df;
  }

