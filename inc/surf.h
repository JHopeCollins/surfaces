# ifndef SURF_H
# define SURF_H

#include <assert.h>
#include "vec.h"
#include "vtx.h"
// basic surface classes

// virtual base surface class defining interface

   template <int DIM_, typename POS_>
   struct surf_t
  {

      surf_t(){}
     ~surf_t(){}

   // calculate normal, gradient and curvature of surface at point closest to p0
      virtual void tangent(   POS_ p0,  POS_ *d ) = 0;
      virtual void normal(    POS_ p0,  POS_ &n ) = 0;
      virtual void curvature( POS_ p0,  POS_ *c ) = 0;

   // q is projection of point p onto surface.
      virtual void project( POS_ p0, POS_ &q0 ) = 0;
      void project( POS_ &p ){ project( p, p ); }

   // returns distance from point p to closest point on surface
      inline double distance( POS_ p0 )
     {
         POS_      q;
         project( p0, q );
         q = q-p0;
         return length( q );
     }
  };


// class for surface defined by parametric coordinates
   template <int DIM_, typename POS_>
   struct surfparametric_t : public surf_t<DIM_, POS_>
  {
      surfparametric_t() {}
     ~surfparametric_t() {}

   // convert from 1D parameterisation to real space
      virtual void x( double s, POS_ &x0 ) = 0;

   // s is arc length parameter of point of surface closest to x
      virtual double s( POS_ x0 ) = 0;
  };


// class for surface defined by implicit function
   template <int DIM_, typename POS_, typename IVAL_>
   struct surfimplicit_t : public surf_t<DIM_, POS_>
  {
      double   accuracy;   // uncertainty margin for projections onto surface

      surfimplicit_t(){}
     ~surfimplicit_t(){}
      surfimplicit_t( double a ) : accuracy(a) {}

   // calculate normal, gradient and curvature of surface at point closest to p0
      void tangent( POS_ p0, POS_ *t )
  {
      POS_ n, x;

      project( p0, x );
      normal(   x, n );
      t[0]=  n[1];
      t[1]= -n[0];
  }

      void normal( POS_ p0, POS_ &n )
  {
      int          i;
      POS_         x;
      double     l=0;
      IVAL_  d[DIM_];

      project( p0, x );
      dF(      x,  d );

      for( i=0; i<DIM_; i++ ){ l+=d[i]*d[i]; }
      l=sqrt(l);

      for( i=0; i<DIM_; i++ ){ n[i]=d[i]/l; }
  }

      void curvature( POS_ p0, POS_ *c )
  {
      POS_   x;
      IVAL_  h[DIM_];

      project( p0, x );
      d2F(      x, h );

      for( int i=0; i<DIM_*DIM_; i++ ){ c[i]=h[i]; }
  }

   // q is projection of point p onto surface.
      void project( POS_ p0, POS_ &q0 )
  {
      if( F( p0 ) < accuracy ){ q0=p0; return; }

      int           i=0;
      int           d=DIM_+1;
      int           nrhs=1;
      int           info=0;
      double        res=10000000.;             // newton iteration residual
      double        rlx=0.3;               // relaxation factor
      double        err=accuracy;    // convergence criteria

      int           ipiv[d];
      double       jac[d*d];
      double   nu[d],rhs[d];

      POS_        x0;
      IVAL_ df[DIM_];

   // initialise initial point x0 and solution vector nu
      x0=p0;
      for( i=0; i<DIM_; i++ ){ nu[i]=p0[i]; }
      nu[DIM_] = 0.;

      int j=0;
   // newton iterations
      while( res>err )
     {
      // construct rhs and jacobian
         newton_rhs( x0, nu, rhs, df );
         newton_jac( x0, nu, jac, df );

         for( i=0; i<d; i++ ){ ipiv[i]=0; }

      //          N   NRHS   A  LDA  IPIV   B  LDB,  info )
         dgesv_( &d, &nrhs, jac, &d, ipiv, rhs, &d, &info );
         assert( info==0 );

      // update solution vector
         for( i=0; i<d; i++ ){ nu[i]+=rlx*rhs[i]; }

      // calculate residual
         res = 0.;
         for( i=0; i<d; i++ ){ res+=rhs[i]*rhs[i]; }
         res=sqrt(res);

         rlx*=1.3;
         rlx=std::min(rlx,1.0);

         j++;
     }

   // assign to q
      for( i=0; i<DIM_; i++ ){ q0[i]=nu[i]; }

      return;
  }

   // implicit function (pseudo-distance function)
      virtual IVAL_ F( POS_ p0 ) = 0;

   // gradient of implicit function
      virtual void dF( POS_ p0, IVAL_ *d ) = 0;

   // second gradient (hessian) of implicit function
      virtual void d2F( POS_ p0, IVAL_ *h ) = 0;

   // right hand side for projection newton iterations
      inline void newton_rhs( POS_ x0, double *nu, double *rhs, IVAL_ *df )
     {
         int i;

         POS_  x; // current position

         for( i=0; i<DIM_; i++ ){ x[i]=nu[i]; }
         dF(  x, df );

         for( i=0; i<DIM_; i++ )
        {
            rhs[i] = -( nu[DIM_]*df[i] - ( x0[i] - nu[i] ) );
        }
         rhs[DIM_] = -F( x );
     }

   // jacobian matrix for projection newton iterations
      inline void newton_jac( POS_ x0, double *nu, double *jac, IVAL_ *df )
     {

         POS_              x;
         int         i,j,k,l;
         IVAL_  h[DIM_*DIM_];   // hessian

         for( i=0; i<DIM_; i++ ){ x[i]=nu[i]; }

         for( i=0; i<(DIM_+1)*(DIM_+1); i++ ){ jac[i]=0.; }

         d2F( x, h );

/*
         printf("hessian\n");
         for(i=0;i<DIM_;i++)
        {
            for(j=0;j<DIM_;j++)
           {
               k=indx(i,j,DIM_);
               printf("%f,  ", h[k]);
           }
            printf("\n");
        }
*/

         for( i=0; i<DIM_; i++ )
        {
            k = indx( DIM_, i,    DIM_+1 );
            l = indx( i,    DIM_, DIM_+1 );
            jac[ k ] = df[i];
            jac[ l ] = df[i];

            for( j=0; j<DIM_; j++ )
           {
               k = indx( i,  j, DIM_+1 );
               l = indx( i,  j, DIM_   );
               jac[ k ]+= nu[DIM_]*h[ l ];
           }

            k= indx( i,i, DIM_+1 );
            jac[ k ]+= 1.;
        }

         jac[(DIM_+1)*(DIM_+1)-1]=0.;

/*
         for(i=0;i<(DIM_+1)*(DIM_+1);i++){printf("%d\n",t[i]);}

         printf("jacobian\n");
         for(i=0;i<DIM_+1;i++)
        {
            for(j=0;j<DIM_+1;j++)
           {
               k=indx(i,j,DIM_+1);
               printf("%f,  ", jac[k]);
           }
            printf("\n");
        }
*/
     }
  };


# endif
