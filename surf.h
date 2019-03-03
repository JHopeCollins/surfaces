# ifndef SURF_H
# define SURF_H

#include <assert.h>
#include "vec.h"

// basic surface classes

// virtual base surface class defining interface
   struct surf_t : public virtual vec_t
  {

      surf_t(){}
      surf_t( int d )
     {
         dims=d;
     }
     ~surf_t(){}

   // calculate normal, gradient and curvature of surface at point closest to p0
      virtual void tangent(    double *p0, double *d ) = 0;
      virtual void normal(     double *p0, double *n ) = 0;
      virtual void curvature(  double *p0, double *c ) = 0;

   // q is projection of point p onto surface.
      virtual void project( double *p0, double *q0 ) = 0;
      void project( double *p ){ project( p, p ); }

   // returns distance from point p to closest point on surface
      inline double distance( double *p0 )
     {
         double   *q=NULL, l=0;
         q=new double[dims];
         eq(       q, 0.);
         project( p0, q );
         sub(     q, p0 );
         l = length(  q );
         delete[] q;
         return l;
     }
  };


// class for surface defined by parametric coordinates
   struct surfparametric_t : public surf_t
  {
      surfparametric_t() {}
     ~surfparametric_t() {}
      surfparametric_t( int d ){ dims=d; }

   // convert from 1D parameterisation to real space
      virtual void x( double s, double *x0 ) = 0;

   // s is arc length parameter of point of surface closest to x
      virtual double s( double *x0 ) = 0;

  };


// class for surface defined by implicit function
   struct surfimplicit_t : public surf_t
  {
      double   accuracy;   // uncertainty margin for projections onto surface

      surfimplicit_t(){}
     ~surfimplicit_t(){}
      surfimplicit_t( int d, double a )
     {
         dims=d;
         accuracy=a;
     }

   // implicit function (pseudo-distance function)
      virtual double F( double *p0 ) = 0;

   // gradient of implicit function
      virtual void dF( double *p0, double *d ) = 0;

   // second gradient (hessian) of implicit function
      virtual void d2F( double *p0, double *h ) = 0;

   // normal vector points along gradient of implicit function local to surface
      void normal( double *p0, double *n )
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
      void tangent( double *p, double *t )   // only for 2D currently
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
      void curvature( double *p0, double *c )
     {
         double   *x;
         x=new double[dims];

         project( p0, x );
         d2F(      x, c );

         delete[] x;
     }

   // right hand side for projection newton iterations
      inline void newton_rhs( double *x0, double *nu, double *rhs, double *df )
     {
         int i;

         dF(  nu, df );
         for( i=0; i<dims; i++ )
        {
            rhs[i] = -( nu[dims]*df[i] - ( x0[i] - nu[i] ) );
        }
         rhs[dims] = -F( nu );
     }

   // jacobian matrix for projection newton iterations
      inline void newton_jac( double *x0, double *nu, double *jac, double *df )
     {

         int      i,j,k,l, *t;
         double   *h;   // hessian

         h=new double[dims*dims];
         t=new    int[(dims+1)*(dims+1)];

         for(i=0;i<(dims+1)*(dims+1);i++){t[i]=0; jac[i]=0.;}

         d2F( nu, h );

/*
         printf("hessian\n");
         for(i=0;i<dims;i++)
        {
            for(j=0;j<dims;j++)
           {
               k=index(i,j,dims);
               printf("%f,  ", h[k]);
           }
            printf("\n");
        }
*/

         for( i=0; i<dims; i++ )
        {
            k = index( dims, i,    dims+1 );
            l = index( i,    dims, dims+1 );
            jac[ k ] = df[i]; t[k]+=1;
            jac[ l ] = df[i]; t[l]+=1;

            for( j=0; j<dims; j++ )
           {
               k = index( i,  j, dims+1 );
               l = index( i,  j, dims   );
               jac[ k ]+= nu[dims]*h[ l ]; t[k]+=1;
           }

            k= index( i,i, dims+1 );
            jac[ k ]+= 1.; t[k]+=1;
        }

         jac[(dims+1)*(dims+1)-1]=0.; t[(dims+1)*(dims+1)-1]+=1;

/*
         for(i=0;i<(dims+1)*(dims+1);i++){printf("%d\n",t[i]);}

         printf("jacobian\n");
         for(i=0;i<dims+1;i++)
        {
            for(j=0;j<dims+1;j++)
           {
               k=index(i,j,dims+1);
               printf("%f,  ", jac[k]);
           }
            printf("\n");
        }
*/

         delete[] h;
         delete[] t;
     }

   // project point onto surface, subject to constraints of F(x)=0 and original point orthogonal to projected point
      void project( double *p, double *q )
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
         const double  rlx=1.0;               // relaxation factor
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
/*
         for( i=0; i<3; i++ )
        {
            dF(  nu, x0  );
            mul( x0, 0.4*F(nu) );
            add( nu, x0 );
        }
*/
         eq(  x0, p );

         nu[dims] = radius( nu, x0 );

         int j=0;
      // start newton iterations
         while( res>err )
        {
//          printf( "projection iteration: %d\n", j );
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
            j++;
        }

      // assign to q
         eq( q, nu );

         double l=vecd.radius( p, q );
         printf( "iteration %d, res= %f, dx= %f\n", j, res, l );

         delete[] ipiv;
         delete[] nu;
         delete[] rhs;
         delete[] jac;
         delete[] x0;
         delete[] df;
     }
  };


# endif
