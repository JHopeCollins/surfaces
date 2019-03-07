#  ifndef _SPLINE_
#  define _SPLINE_

#  include <typd.h>
#  include <cmath>
#  include <cassert>

/** Spline fitting functions **/

/** Natural spline through n points with abscissa x and values y **/

   template <typename TYPE_> void spline( INT_ n, REAL_ *x, TYPE_ *y, TYPE_ *r )
  {
      INT_   i;
      REAL_  a,b,w;
      REAL_ *c;

      c= new REAL_[n];

      r[0]= TYPE_(0.);
      c[0]= 0.;

      for(i=1;i<n-1;i++) 
     {

// evaluate matrix coefficients

         a=     x[i  ]- x[i-1];
         b=     x[i+1]- x[i-1];
         b*=    2;
         w=     x[i+1]- x[i  ];

         r[i]=  ( y[i+1]-y[i  ] )/w- ( y[i  ]-y[i-1] )/a;
         r[i]*= 6;

// in-place inversion

         b-=    a*c[i-1];
         r[i]-= a*r[i-1];

         b= 1./b;
         r[i]*= b;
         c[i]=  w*b;
         
     }

      r[n-1]=  TYPE_(0.);

// backward sweep

      for(i=n-2;i>=0;i--)
     {
         r[i]-= c[i]*r[i+1];
     }

      delete[] c; c= NULL; 
      return;
  };

/** Constrained spline through n points with abscissa x, values y and slopes d0,d1 (if specified) **/

   template <typename TYPE_> void spline( INT_ n, REAL_ *x, TYPE_ *y, bool b0, TYPE_ d0, bool b1, TYPE_ d1, TYPE_ *r )
  {
      INT_   i;
      REAL_  a,b,w;
      REAL_ *c;

      c= new REAL_[n];

      if( !b0 )
     {
         r[0]= TYPE_(0.);
         c[0]= TYPE_(0.);
     }
      else
     {
         w= x[1]-x[0];
         
         r[0]=  (y[1]-y[0])/w;
         r[0]-= d0;
         r[0]*= 3;

         w= 1./w;
         r[0]*= w;
         c[0] = 0.5;
     }

      for(i=1;i<n-1;i++) 
     {

// evaluate matrix coefficients

         a=     x[i  ]- x[i-1];
         b=     x[i+1]- x[i-1];
         b*=    2;
         w=     x[i+1]- x[i  ];

         r[i]=  ( y[i+1]-y[i  ] )/w- ( y[i  ]-y[i-1] )/a;
         r[i]*= 6;

// in-place inversion

         b-=    a*c[i-1];
         r[i]-= a*r[i-1];

         b= 1./b;
         r[i]*= b;
         c[i]=  w*b;
         
     }
      if( !b1 )
     {
         r[n-1]=  TYPE_(0.);
         a= 0;
         b= 1;
     }
      else 
     {
         a= x[n-1]- x[n-2];
         b= a; 
         b*=2;

         r[n-1]= d1;
         r[n-1]-= (y[n-1]-y[n-2])/a;
         r[n-1]*= 6;
     }

      b-=      a*c[n-2];
      r[n-1]-= a*r[n-2];

      b= 1./b;
      r[n-1]*= b;

// backward sweep

      for(i=n-2;i>=0;i--)
     {
         r[i]-= c[i]*r[i+1];
     }

      delete[] c; c= NULL; 
      return;
  };


/** Periodic spline through n points with abscissa x, period t and values y **/

   template <typename TYPE_> void spline( INT_ n, REAL_ t, REAL_ *x, TYPE_ *y, TYPE_ *r )
  {
      INT_   i;
      REAL_  a,b,w;
      REAL_ *c;
      REAL_ *s;

      c= new REAL_[n];
      s= new REAL_[n];

      a=   t+x[0  ]- x[n-1];
      b=   t+x[1  ]- x[n-1];
      b*=    2;
      w=     x[1  ]- x[0  ];

      r[0]=  ( y[1  ]-y[0  ] )/w- ( y[0  ]-y[n-1] )/a;
      r[0]*= 6;
      s[0]=  a;

      b= 1./b;
      r[0]*= b;
      s[0]*= b;
      c[0]=  w*b;

      for(i=1;i<n-2;i++) 
     {

// evaluate matrix coefficients

         a=     x[i  ]- x[i-1];
         b=     x[i+1]- x[i-1];
         b*=    2;
         w=     x[i+1]- x[i  ];

         r[i]=  ( y[i+1]-y[i  ] )/w- ( y[i  ]-y[i-1] )/a;
         r[i]*= 6;
         s[i]=  0.;

// in-place inversion

         b-=    a*c[i-1];
         r[i]-= a*r[i-1];
         s[i]-= a*s[i-1];

         b= 1./b;
         r[i]*= b;
         s[i]*= b;
         c[i]=  w*b;
         
     }

// evaluate matrix coefficients
      i= n-2;

      a=     x[i  ]- x[i-1];
      b=     x[i+1]- x[i-1];
      b*=    2;
      w=     x[i+1]- x[i  ];

      r[i]=  ( y[i+1]-y[i  ] )/w- ( y[i  ]-y[i-1] )/a;
      r[i]*= 6;
      s[i]=  w;

// in-place inversion

      b-=    a*c[i-1];
      r[i]-= a*r[i-1];
      s[i]-= a*s[i-1];

      b= 1./b;
      r[i]*= b;
      s[i]*= b;

// backward sweep

      for(i=n-3;i>=0;i--)
     {
         r[i]-= c[i]*r[i+1];
         s[i]-= c[i]*s[i+1];
     }

// Schur complement section

      i= n-1;

      a=     x[n-1]- x[n-2];
      b=   t+x[0  ]- x[n-2];
      b*=    2;
      w=   t+x[0  ]- x[n-1];

      r[n-1]=  ( y[0  ]-y[n-1  ] )/w- ( y[n-1  ]-y[n-2] )/a;
      r[n-1]*= 6;

      r[n-1]-= r[0  ]*w;
      r[n-1]-= r[n-2]*a;

      b-= s[0  ]*w;
      b-= s[n-2]*w;
 
      b= 1./b;
      r[n-1]*= b;

      for( i=0;i<n-1;i++ ) 
     {
         r[i]-= r[n-1]*s[i];
     }

      delete[] c; c= NULL; 
      delete[] s; s= NULL; 

      return;
  };


   template <typename TYPE_> void splint0( REAL_ x0, REAL_ x1, TYPE_ y0, TYPE_ y1, TYPE_ r0, TYPE_ r1, REAL_ u, TYPE_ &v )
  {
      REAL_ a,b, c,d;
      REAL_ h;

      h=x1-x0;
      assert( h > 0 );

      a=(x1-u)/h;
      b=(u-x0)/h;

      c= a*a;
      c*= a;
      c-= a;

      d= b*b;
      d*= b;
      d-= b;

      h*= h;
      h/= 6;

      v=   a* y0+ b* y1;
      v+= (c* r0+ d* r1)*h;
  }

   template <typename TYPE_> void splint0( REAL_ x0, REAL_ x1, TYPE_ y0, TYPE_ y1, TYPE_ r0, TYPE_ r1, REAL_ u, TYPE_ &v, TYPE_ &dv )
  {
      REAL_ a,b, c,d, e,f;
      REAL_ h,w,k;

      h=x1-x0;
      assert( h > 0 );
      w= 1./h;
      k= h/6;

      a=(x1-u)*w;
      b=(u-x0)*w;

      c=  a*a;
      e=  3*c;
      c*= a;
      c-= a;
      e--;

      d= b*b;
      f= 3*d;
      d*= b;
      d-= b;
      f--;

      dv= y1- y0;
      dv*=w;
      dv+=(f*r1-e*r0)*k;

      h*= k;

      v=   a* y0+ b* y1;
      v+= (c*r0+ d*r1)*h;

  }

   template <typename TYPE_> void splint0( REAL_ x0, REAL_ x1, TYPE_ y0, TYPE_ y1, TYPE_ r0, TYPE_ r1, REAL_ u, TYPE_ &v, TYPE_ &dv, TYPE_ &d2v )
  {
      REAL_ a,b, c,d, e,f;
      REAL_ h,w,k;

      h=x1-x0;
      assert( h > 0 );
      w= 1./h;
      k= h/6.;

      a=(x1-u)*w;
      b=(u-x0)*w;

      c=  a*a;
      e=  3*c;
      c*= a;
      c-= a;
      e--;

      d= b*b;
      f= 3*d;
      d*= b;
      d-= b;
      f--;

      dv= y1- y0;
      dv*=w;
      dv+=(f*r1- e*r0)*k;

      d2v= b*r1+ a*r0;

      h*= k;

      v=   a* y0+ b* y1;
      v+= (c*r0+ d*r1)*h;
  }

/** Spline interpolation functions **/

/** Open splines **/

   template <typename TYPE_> void splint(INT_ n, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v )
  {

      INT_ k0,k1,k;
      TYPE_ dv;

      if( u < 0 )
     {
         splint0( x[0],x[1], y[0],y[1], r[0],r[1], 0.,v,dv );
         v+= u*dv;
     }
      else
     {
         if( u <= 1 )
        {
            k0=0;
            k1=n-1;
       
            while( k1-k0 > 1) 
           {
               k=(k1+k0) >> 1;
               if (x[k] > u){ k1=k; }else{ k0=k; }
           }
            splint0( x[k0],x[k1], y[k0],y[k1], r[k0],r[k1], u,v );
        }
         else
        {
            splint0( x[n-2],x[n-1], y[n-2],y[n-1], r[n-2],r[n-1], 1.,v,dv );
            v+= (u-1.)*dv;
        }
     }


      return;
  }

   template <typename TYPE_> void splint(INT_ n, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v, TYPE_ &dv )
  {

      INT_ k0,k1,k;
      if( u < 0 )
     {
         splint0( x[0],x[1], y[0],y[1], r[0],r[1], 0.,v,dv );
         v+= u*dv;
     }
      else
     {
         if( u < 1 )
        {
            k0=0;
            k1=n-1;
            while( k1-k0 > 1) 
           {
               k=(k1+k0) >> 1;
               if (x[k] > u){ k1=k; }else{ k0=k; }
           }
            splint0( x[k0],x[k1], y[k0],y[k1], r[k0],r[k1], u,v,dv );
        }
         else
        {
            splint0( x[n-2],x[n-1], y[n-2],y[n-1], r[n-2],r[n-1], 1.,v,dv );
            v+= (u-1.)*dv;
        }
     }


      return;
  }

   template <typename TYPE_> void splint(INT_ n, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v, TYPE_ &dv, TYPE_ &d2v )
  {

      INT_ k0,k1,k;

      if( u < 0 )
     {
         splint0( x[0],x[1], y[0],y[1], r[0],r[1], 0.,v,dv );
         v+= u*dv;
         d2v=TYPE_(0.);
     }
      else
     {
         if( u < 1 )
        {
            k0=0;
            k1=n-1;
            while( k1-k0 > 1) 
           {
               k=(k1+k0) >> 1;
               if (x[k] > u){ k1=k; }else{ k0=k; }
           }
            splint0( x[k0],x[k1], y[k0],y[k1], r[k0],r[k1], u,v,dv,d2v );
        }
         else
        {
            splint0( x[n-2],x[n-1], y[n-2],y[n-1], r[n-2],r[n-1], 1.,v,dv );
            v+= (u-1.)*dv;
            d2v=TYPE_(0.);
        }
     }

      return;
  }

/** Periodic splines **/

   template <typename TYPE_> void splint(INT_ n, REAL_ t, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v )
  {

      INT_ k0,k1,k;
      REAL_    x0,x1;
      TYPE_    y0,y1, r0,r1;

      while( u <     x[0] ){ u+= t; };
      while( u >=  t+x[0] ){ u-= t; };

      if( u > x[n-1] )
     {
         x0= x[n-1];
         x1= x[0  ]+t;
         y0= y[n-1];
         y1= y[0  ];
         r0= r[n-1];
         r1= r[0  ];
     } 
      else
     {
         k0=0;
         k1=n-1;
         while( k1-k0 > 1) 
        {
            k=(k1+k0) >> 1;
            if (x[k] > u){ k1=k; }else{ k0=k; }
        }
         x0= x[k0];
         x1= x[k1];
         y0= y[k0];
         y1= y[k1];
         r0= r[k0];
         r1= r[k1];
     }
 
      splint0( x0,x1, y0,y1, r0,r1, u,v );

      return;
  }

   template <typename TYPE_> void splint(INT_ n, REAL_ t, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v, TYPE_ &dv )
  {

      INT_ k0,k1,k;
      REAL_    x0,x1;
      TYPE_    y0,y1, r0,r1;

      while( u <     x[0] ){ u+= t; };
      while( u >=  t+x[0] ){ u-= t; };

      if( u > x[n-1] )
     {
         x0= x[n-1];
         x1= x[0  ]+t;
         y0= y[n-1];
         y1= y[0  ];
         r0= r[n-1];
         r1= r[0  ];
     } 
      else
     {
         k0=0;
         k1=n-1;
         while( k1-k0 > 1) 
        {
            k=(k1+k0) >> 1;
            if (x[k] > u){ k1=k; }else{ k0=k; }
        }
         x0= x[k0];
         x1= x[k1];
         y0= y[k0];
         y1= y[k1];
         r0= r[k0];
         r1= r[k1];
     }
 
      splint0( x0,x1, y0,y1, r0,r1, u,v,dv );

      return;
  }

   template <typename TYPE_> void splint(INT_ n, REAL_ t, REAL_ *x, TYPE_ *y, TYPE_ *r, REAL_ u, TYPE_ &v, TYPE_ &dv, TYPE_ &d2v )
  {

      INT_ k0,k1,k;
      REAL_    x0,x1;
      TYPE_    y0,y1, r0,r1;

      while( u <     x[0] ){ u+= t; };
      while( u >=  t+x[0] ){ u-= t; };

      if( u > x[n-1] )
     {
         x0= x[n-1];
         x1= x[0  ]+t;
         y0= y[n-1];
         y1= y[0  ];
         r0= r[n-1];
         r1= r[0  ];
     } 
      else
     {
         k0=0;
         k1=n-1;
         while( k1-k0 > 1) 
        {
            k=(k1+k0) >> 1;
            if (x[k] > u){ k1=k; }else{ k0=k; }
        }
         x0= x[k0];
         x1= x[k1];
         y0= y[k0];
         y1= y[k1];
         r0= r[k0];
         r1= r[k1];
     }
 
      splint0( x0,x1, y0,y1, r0,r1, u,v,dv,d2v );

      return;
  }
#  endif

