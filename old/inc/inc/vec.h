# ifndef VEC_H
# define VEC_H

# include <math.h>
# include <algorithm>
# include <stdlib.h>

   extern "C" void dgesv_(         int *, int *, double *, int *, int *, double *,                  int *, int * );
   extern "C" void dsysv_( char *, int *, int *, double *, int *, int *, double *, int *, double *, int *, int * );

// index of 1D array aliasing column-major 2D array a[i][j] -> b[ index(i,j,imax) ]
   inline int indx( int i, int j, int n ){ return i*n+j; }

// base for any struct which needs to hold a dimension
   struct hasdimension_t
  {
      public:
      int   dims;

      hasdimension_t() {}
      hasdimension_t( int d ) : dims(d) {};
     ~hasdimension_t() {}
  };

// base struct for classes requiring vector operations

   struct vec_t : public virtual hasdimension_t
  {
      double   ***e;
      int      arrays=0;

      vec_t() {e=NULL;}
      vec_t( int d ){ dims=d; construct_levi_civita(); }

      virtual ~vec_t() { if( arrays ){ free(); } }

   // memory management for arrays
      virtual void malloc();

      virtual void free();

      void construct_levi_civita();

      inline void eq( double *u, double *v )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]; }
     }

      inline void eq( double *u, double  a )
     {
         for( int i=0; i<dims; i++ ){ u[i]=a; }
     }

   // multiply by factor
      inline void mul( double *u, double a, double *v )
     {
         for( int i=0; i<dims; i++ ){ u[i]=a*v[i]; }
     }

      inline void mul( double *u, double *v, double *w )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]*w[i]; }
     }

      inline void mul( double *u, double  a ){ mul( u, a, u ); }
      inline void mul( double *u, double *v ){ mul( u, u, v ); }

   // divide by factor
      inline void div( double *u, double a, double *v )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]/a; }
     }

      inline void div( double *u, double *v, double *w )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]/w[i]; }
     }

      inline void div( double *u, double  a ){ div( u, a, u ); }
      inline void div( double *u, double *v ){ div( u, u, v ); }

   // add constant
      inline void add( double *u, double a, double *v )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]+a; }
     }

      inline void add( double *u, double *v, double *w )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]+w[i]; }
     }

      inline void add( double *u, double  a ){ add( u, a, u ); }
      inline void add( double *u, double *v ){ add( u, u, v ); }

   // u = a*v+w
      inline void muladd( double *u, double a, double *v, double *w )
     {
         double *t=NULL; t=new double[dims];
         mul( t, a, v );
         add( u, t, w );
         delete[] t; t=NULL;
     }

   // subtract constant
      inline void sub( double *u, double a, double *v )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]-a; }
     }

      inline void sub( double *u, double *v, double *w )
     {
         for( int i=0; i<dims; i++ ){ u[i]=v[i]-w[i]; }
     }

      inline void sub( double *u, double  a ){ sub( u, a, u ); }
      inline void sub( double *u, double *v ){ sub( u, u, v ); }

   // dot product of u and v
      inline double dot( double *u, double *v )
     {
         double val=0.;
         for( int i=0; i<dims; i++ ){ val+=u[i]*v[i]; }
         return val;
     }

   // euclidian 2-norm of u
      inline double length( double *u ){ return sqrt( dot( u, u ) ); }

   // euclidian inf-norm of u
      inline double infnorm( double *u )
     {
         double   val=abs(u[0]);
         for( int i=1; i<dims; i++ ){ val=std::max( val, abs(u[i]) ); }
         return val;
     }

   // n is unit vector in direction u
      inline void unit( double *u, double *n ){ div( n, length(u), u ); }
      inline void unit( double *u ){ unit(u,u); }

   // euclidian 2-norm of distance between points p0 and p1
      inline double radius( double *p0, double *p1 )
     {
         double *d=NULL, r=0;

         d=new double[dims];

         sub( d, p0, p1 );
         r=length( d );

         delete[] d; d=NULL;

         return r;
     }

      inline void cross(  double *u, double *v, double *w )
     {
         switch( dims )
        {
            case 2: cross2( u, v, w ); return;
            case 3: cross3( u, v, w ); return;
        }

         int       i,j,k;
         double   *n=NULL;
         n=new double[dims];

         eq(n, 0.);

         for( i=0; i<dims; i++ )
        {
            for( j=0; j<dims; j++ )
           {
               for( k=0; k<dims; k++ )
              {
                  n[k]+= e[i][j][k]*v[i]*w[j];
              }
           }
        }
         eq( u, n );

         delete[] n;
         n=NULL;
     }

      inline void cross2( double *u, double *v, double *w )
     {
         u[0] = v[0]*w[1]-v[1]*w[0];
         u[1] = 0.;
         u[2] = 0.;
     }

      inline void cross3( double *u, double *v, double *w )
     {
         u[0] = v[1]*w[2] - v[2]*w[1];
         u[1] = v[2]*w[0] - v[0]*w[2];
         u[2] = v[0]*w[1] - v[1]*w[0];
     }

      inline void cross(  double *u, double *v ){ cross( u, u, v ); }
 };

#endif
