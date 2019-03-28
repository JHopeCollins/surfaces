#  ifndef _SPLINE2_
#  define _SPLINE2_

#  include <spline.h>

/** Open splines **/

   template <typename TYPE_> void spline( INT_ m, INT_ n, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2 )
  {
      INT_ j,k;

      k=0;
      for(j=0;j<n;j++)
     { 
         spline(m, x,z+k, z2+k); 
         k+=m;
     }
  }

   template <typename TYPE_> void splint( INT_ m, INT_ n, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v )
  {
      INT_   j,k;
      TYPE_ *w,*w2;

      w=  new TYPE_[n];
      w2= new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, x, z+k,z2+k, s, w[j]);
         k+= m;
     }
      spline( n, y, w,w2);
      splint( n, y, w,w2, t,v);

      delete[] w;
      delete[] w2;
  }

   template<typename TYPE_> void splint( INT_ m, INT_ n, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v, TYPE_ *dv )
  {
      INT_   j,k;
      TYPE_  *w2, *w;
      TYPE_ *dw2,*dw;

      w=  new TYPE_[n];
      w2= new TYPE_[n];

      dw= new TYPE_[n];
      dw2=new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, x, z+k,z2+k, s, w[j],dw[j] );
         k+= m;
     }
      spline( n, y,  w, w2);
      spline( n, y, dw,dw2);

      splint( n, y,  w, w2, t,v,dv[1]);
      splint( n, y, dw,dw2, t,  dv[0]);

      delete[]   w;
      delete[]  w2;
      delete[]  dw;
      delete[] dw2;
  }

   template <typename TYPE_> void splint( INT_ m, INT_ n, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v, TYPE_ *dv, TYPE_ *d2v )
  {
      INT_   j,k;

      TYPE_   *w2  ,*w;
      TYPE_  *dw2, *dw;
      TYPE_ *d2w2,*d2w;

      w=    new TYPE_[n];
      w2=   new TYPE_[n];

      dw=   new TYPE_[n];
      dw2=  new TYPE_[n];

      d2w=  new TYPE_[n];
      d2w2= new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, x, z+k,z2+k, s, w[j],dw[j],d2w[j]);
         k+= m;
     }

      spline( n, y,   w,  w2 );
      spline( n, y,  dw, dw2 );
      spline( n, y, d2w,d2w2 );

      splint( n, y,   w,  w2, t,   v,    dv[1],d2v[2]);
      splint( n, y,  dw, dw2, t,  dv[0],d2v[1]);
      splint( n, y, d2w,d2w2, t, d2v[0]);

      delete[]    w;
      delete[]   w2;
      delete[]   dw;
      delete[]  dw2;
      delete[]  d2w;
      delete[] d2w2;
  }

/** Open splines **/

   template <typename TYPE_> void spline( INT_ m, INT_ n, REAL_ p, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2 )
  {
      INT_ j,k;

      k=0;
      for(j=0;j<n;j++)
     { 
         spline(m, p,x,z+k, z2+k); 
         k+=m;
     }
  }

   template <typename TYPE_> void splint( INT_ m, INT_ n, REAL_ p, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v )
  {
      INT_   j,k;
      TYPE_ *w,*w2;

      w=  new TYPE_[n];
      w2= new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, p,x, z+k,z2+k, s, w[j]);
         k+= m;
     }
      spline( n, y, w,w2);
      splint( n, y, w,w2, t,v);

      delete[] w;
      delete[] w2;
  }

   template<typename TYPE_> void splint( INT_ m, INT_ n, REAL_ p, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v, TYPE_ *dv )
  {
      INT_   j,k;
      TYPE_  *w2, *w;
      TYPE_ *dw2,*dw;

      w=  new TYPE_[n];
      w2= new TYPE_[n];

      dw= new TYPE_[n];
      dw2=new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, p,x, z+k,z2+k, s, w[j],dw[j] );
         k+= m;
     }
      spline( n, y,  w, w2);
      spline( n, y, dw,dw2);

      splint( n, y,  w, w2, t,v,dv[1]);
      splint( n, y, dw,dw2, t,  dv[0]);

      delete[]   w;
      delete[]  w2;
      delete[]  dw;
      delete[] dw2;
  }

   template <typename TYPE_> void splint( INT_ m, INT_ n, REAL_ p, REAL_ *x, REAL_ *y, TYPE_ *z, TYPE_ *z2, REAL_ s, REAL_ t, TYPE_ &v, TYPE_ *dv, TYPE_ *d2v )
  {
      INT_   j,k;

      TYPE_   *w2  ,*w;
      TYPE_  *dw2, *dw;
      TYPE_ *d2w2,*d2w;

      w=    new TYPE_[n];
      w2=   new TYPE_[n];

      dw=   new TYPE_[n];
      dw2=  new TYPE_[n];

      d2w=  new TYPE_[n];
      d2w2= new TYPE_[n];

      k= 0;
      for(j=0;j<n;j++)
     {
         splint( m, p,x, z+k,z2+k, s, w[j],dw[j],d2w[j]);
         k+= m;
     }

      spline( n, y,   w,  w2 );
      spline( n, y,  dw, dw2 );
      spline( n, y, d2w,d2w2 );

      splint( n, y,   w,  w2, t,   v,    dv[1],d2v[2]);
      splint( n, y,  dw, dw2, t,  dv[0],d2v[1]);
      splint( n, y, d2w,d2w2, t, d2v[0]);

      delete[]    w;
      delete[]   w2;
      delete[]   dw;
      delete[]  dw2;
      delete[]  d2w;
      delete[] d2w2;
  }


#endif

