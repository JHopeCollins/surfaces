#  ifndef _VTX_
#  define _VTX_

#  include <typd.h>
#  include <cstring>
#  include <cmath>

   struct vtx_t
  {
      REAL_3       x;
      REAL_& operator[]( INT_ i ){ return x[i]; };

      vtx_t(){ x[0]=0;x[1]=0;x[2]=0; };
      vtx_t( REAL_ a ){ x[0]=a;x[1]=a;x[2]=a; };
      vtx_t( REAL_ a, REAL_ b, REAL_ c ){ x[0]=a;x[1]=b;x[2]=c; };
      vtx_t( REAL_3 v ){ x[0]=v[0];x[1]=v[1];x[2]=v[2]; };

      inline vtx_t &operator+=( vtx_t b  )
     {
         x[0]+= b[0];
         x[1]+= b[1];
         x[2]+= b[2];
         return *this;
     };
      inline vtx_t &operator-=( vtx_t b  )
     {
         x[0]-= b[0];
         x[1]-= b[1];
         x[2]-= b[2];
         return *this;
     };
   
   
      inline vtx_t &operator*=( REAL_ a  )
     {
         x[0]*=a;
         x[1]*=a;
         x[2]*=a;
         return *this;
     };

      inline vtx_t &operator/=( REAL_ a  )
     {
         REAL_ b=1./a;
         x[0]*=b;
         x[1]*=b;
         x[2]*=b;
         return *this;
     };
      
  };


   struct bvtx_t  // : public vtx_t
  {
      INT_         id;   
      INT_         ref;   
      REAL_2       y;

      bvtx_t()
     {
         id=-1;
         ref=0;
         memset( y,-1,sizeof(REAL_2));
     };
  };

   inline vtx_t operator*( REAL_ a, vtx_t b )
  {
      vtx_t c;
      c[0]= a*b[0];
      c[1]= a*b[1];
      c[2]= a*b[2];
      return c;
  };

   inline vtx_t operator/( vtx_t b, REAL_ a )
  {
      vtx_t c;
      REAL_ d= 1./a;
      c[0]= d*b[0];
      c[1]= d*b[1];
      c[2]= d*b[2];
      return c;
  };

   inline vtx_t operator+( vtx_t a, vtx_t b )
  {
      vtx_t c;
      c[0]= a[0]+ b[0];
      c[1]= a[1]+ b[1];
      c[2]= a[2]+ b[2];
      return c;
  };

   inline vtx_t operator-( vtx_t a, vtx_t b )
  {
      vtx_t c;
      c[0]= a[0]- b[0];
      c[1]= a[1]- b[1];
      c[2]= a[2]- b[2];
      return c;
  };


   inline REAL_ operator*( vtx_t a, vtx_t b )
  {
      return a[0]*b[0]+a[1]*b[1]+ a[2]*b[2];
  };


   inline vtx_t operator*( vtx_t a, REAL_ b )
  {
      vtx_t val( a[0]*b, a[1]*b, a[2]*b );
      return val;
  };

/* inline vtx_t operator*( vtx_t a, INT_ b )
  {
      vtx_t val( a[0]*b, a[1]*b, a[2]*b );
      return val;
  };*/


   inline vtx_t cross( vtx_t a, vtx_t b )
  {
      vtx_t c;
      c[0]= a[1]*b[2]- a[2]*b[1];
      c[1]= a[2]*b[0]- a[0]*b[2];
      c[2]= a[0]*b[1]- a[1]*b[0];
      return c;
  };

   inline REAL_ norminf( vtx_t a )
  {
      REAL_ val=-1;
      val= fabs(a[0]);
      val= fmax( val,fabs(a[1]));
      val= fmax( val,fabs(a[2]));
      return val;
  };


   inline vtx_t plane_normal( vtx_t v0, vtx_t v1, vtx_t v2, vtx_t v3 )
  {
      vtx_t    d0, d1;
      vtx_t    n;

      INT_     i;
      REAL_3   buf;

      for( i=0; i<3; i++ ){ buf[i] = v3[i]-v0[i]; }
      d0 = vtx_t( buf );

      for( i=0; i<3; i++ ){ buf[i] = v2[i]-v1[i]; }
      d1 = vtx_t( buf );

      n = cross( d0, d1 );

      return n;
  };

   inline REAL_ length( vtx_t a )
  {
      REAL_ val;

      val  = a[0]*a[0];
      val += a[1]*a[1];
      val += a[2]*a[2];

      val = sqrt( val );

      return val;
  };


   inline void ldl2( REAL_ *a, REAL_ *x )
  {
      REAL_ w;

// factorize

      w= a[1];
      a[0]= 1./a[0];
      a[1]*=   a[0];
      a[2]-= w*a[1];
      a[2]= 1./a[2];
         
// solve

      x[1]-= a[1]*x[0];
      x[0]*= a[0];
      x[1]*= a[2];
      x[0]-= a[1]*x[1];

      return;
  }

   inline void orth( vtx_t p, vtx_t l, vtx_t q, vtx_t &r )
  {
      REAL_ d;
      vtx_t s;
      
      r= q-p;

      d= r*l;

      s= l;
      s*=d;

      r-= s;
      r+= p;
      
  };

   struct aprj_t
  {
      INT_        it;
      INT_       mit;
      REAL_      tol;
      REAL_      err;
      REAL_        d;
      REAL_        b[2];
      REAL_        a[3];
      vtx_t        x;
      bool      hint;
      bool       dbg;

      aprj_t()
     {
         tol=1.e-16;
         mit=10;
         d=0;
         b[0]=0;
         b[1]=0;
         a[0]=0;
         a[1]=0;
         a[2]=0;

         dbg=false;
     }
  };

   inline void aprj( vtx_t &x0, vtx_t &x1, vtx_t &x2, vtx_t &x3, vtx_t &p, REAL_ *y, aprj_t &stat )
  {
      REAL_   q0[4];
      REAL_   q1[2][4];
      REAL_   q2[4];

      REAL_   rlx;
      REAL_   err;

      REAL_    d;
      REAL_    b[2];
      REAL_    a[3];

      vtx_t   dx[3];
      vtx_t   d2x;

      vtx_t   x;
      vtx_t   r;

      INT_    j;

      rlx= 0.5;

      for( j=0;j<stat.mit;j++ )
     {

// Bilinear shape function

         q0[0]= (1-y[0])*(1-y[1]);
         q0[1]=    y[0] *(1-y[1]);
         q0[2]= (1-y[0])*   y[1] ;
         q0[3]=    y[0] *   y[1] ;

         q1[0][0]= -(1-y[1]); q1[1][0]= -(1-y[0]);
         q1[0][1]=  (1-y[1]); q1[1][1]=    -y[0] ;
         q1[0][2]=    -y[1] ; q1[1][2]=  (1-y[0]);
         q1[0][3]=     y[1] ; q1[1][3]=     y[0] ;

         q2[0]=  1;
         q2[1]= -1;
         q2[2]= -1;
         q2[3]=  1;

// Position
         x= q0[0]*x0+ q0[1]*x1+ q0[2]*x2+ q0[3]*x3;

// Tangent vectors
         dx[0]= q1[0][0]*x0+ q1[0][1]*x1+ q1[0][2]*x2+ q1[0][3]*x3;
         dx[1]= q1[1][0]*x0+ q1[1][1]*x1+ q1[1][2]*x2+ q1[1][3]*x3;

// Curvature tensor
         d2x=   q2[0]*   x0+ q2[1]*   x1+ q2[2]*   x2+ q2[3]*   x3;

         r= x-p;

// Distance
         d=  r*r;
         d*= 0.5;

// Derivative of distance wrt parametric coordinates

         b[0]= r*dx[0];
         b[1]= r*dx[1];

// Second derivatives

         a[0]=  dx[0]*dx[0];
         a[1]=  dx[0]*dx[1]+ d2x*r;
         a[2]=  dx[1]*dx[1];

// LDL' factorization and solution

         ldl2( a,b );

// Newton update

         y[0]-= rlx*b[0];
         y[1]-= rlx*b[1];

// Convergence criterion

         err= fabs(b[0])+ fabs(b[1]);
         if( err < stat.tol ){ break; };

         rlx*= 1.2;
         rlx= fmin( rlx,1. );

//       printf( "%2d % 9.3e % 9.3e % 9.3e % 9.3e % 9.3e\n", j, d, y[0],y[1], b[0],b[1] );

     }
      stat.a[0]= a[0];
      stat.a[1]= a[1];
      stat.a[2]= a[2];

      stat.x   = x;

      stat.b[0]= b[0];
      stat.b[1]= b[1];

      stat.d=    d;

      stat.err=  err;
      stat.it=   j;

      return;
  }

#  endif
