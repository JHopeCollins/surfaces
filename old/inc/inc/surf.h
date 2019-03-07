# ifndef SURF_H
# define SURF_H

#include <assert.h>
#include "vec.h"
// basic surface classes

// virtual base surface class defining interface
   struct surf_t : public virtual vec_t
  {

      surf_t(){}
      surf_t( int d ){ dims=d; }
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
      surfimplicit_t( int d, double a );

   // calculate normal, gradient and curvature of surface at point closest to p0
      void tangent(    double *p0, double *d );
      void normal(     double *p0, double *n );
      void curvature(  double *p0, double *c );

   // q is projection of point p onto surface.
      void project( double *p0, double *q0 );

   // implicit function (pseudo-distance function)
      virtual double F( double *p0 ) = 0;

   // gradient of implicit function
      virtual void dF( double *p0, double *d ) = 0;

   // second gradient (hessian) of implicit function
      virtual void d2F( double *p0, double *h ) = 0;

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
               k=indx(i,j,dims);
               printf("%f,  ", h[k]);
           }
            printf("\n");
        }
*/

         for( i=0; i<dims; i++ )
        {
            k = indx( dims, i,    dims+1 );
            l = indx( i,    dims, dims+1 );
            jac[ k ] = df[i]; t[k]+=1;
            jac[ l ] = df[i]; t[l]+=1;

            for( j=0; j<dims; j++ )
           {
               k = indx( i,  j, dims+1 );
               l = indx( i,  j, dims   );
               jac[ k ]+= nu[dims]*h[ l ]; t[k]+=1;
           }

            k= indx( i,i, dims+1 );
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
               k=indx(i,j,dims+1);
               printf("%f,  ", jac[k]);
           }
            printf("\n");
        }
*/

         delete[] h;
         delete[] t;
     }
  };


# endif
