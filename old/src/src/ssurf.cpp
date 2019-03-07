#  include <ssurf.h>

   void ssurf_t::read( char *name, tran_t &r, REAL_ h )
  {
      INT_       i;
      INT_       n;
      REAL_      a,b,c;

      FILE      *f;

      f= fopen( name,"r" ); assert( f );
 
      fscanf( f, "%d", &n );
      
      x0.n= n; x0.resize(-1.);
      x1.n= 2; x1.resize(-1.);

      y.n=  2*n;  y.resize(vdum);
      y2.n= 2*n; y2.resize(vdum);

      for( INT_ i=0;i<n;i++ )
     {
         fscanf( f, "%lf %lf %lf", &a,&b,&c );
         y[i][0]= a;
         y[i][1]= b;
         y[i][2]= c;
         r.tran( y[i] ); 
         y[i+n][0]= y[i][0];
         y[i+n][1]= y[i][1];
         y[i+n][2]= y[i][2]+ h;

     }
      fclose( f );

      x0[0]=0;
      for( i=1;i<n;i++ )
     { 
         x0[i]= x0[i-1]+ length( y[i]-y[i-1] ); 
     }
      for( i=1;i<n-1;i++ ){ x0[i]/= x0[n-1]; };
      x0[n-1]= 1;

      x1[0]= 0;
      x1[1]= 1;

      spline( n,2, x0.data(),x1.data(), y.data(),y2.data() );
  }

   void ssurf_t::write( INT_ b, char *name )
  {
      INT_     i,j,k;

      printf( "%2d sur 5 %s %d %d ", b, name, x0.n,x1.n );
      for( i=0;i<x0.n;i++ ){ printf( "%12.5e ", x0[i] ); }
      for( i=0;i<x1.n;i++ ){ printf( "%12.5e ", x1[i] ); }

      k=0; 
      for( j=0;j<x1.n;j++ )
     { 
         for( i=0;i<x0.n;i++ )
        { 
            printf( "% 12.5e % 12.5e % 12.5e ", y[k][0], y[k][1], y[k][2] ); 
            k++;
        }
     }
     
      printf( "\n" );
  }

   void ssurf_t::test( char *name )
  {
      INT_       i,j;

      vtx_t      x;

      FILE      *f;

      f= fopen( name,"w" ); assert( f );

      for( j=0;j<64;j++ )
     {
         for( i=0;i<64;i++ )
        {
            REAL_ u= i; u/= 63.;
            REAL_ v= j; v/= 63.;
            pos( u,v,x );
            fprintf( f,"% 9.3e % 9.3e % 9.3e % 9.3e % 9.3e\n", u,v,x[0],x[1],x[2] );
        }
         fprintf( f,"\n" );
     }
      fclose( f );
  }

   void ssurf_t::pos( REAL_ u, REAL_ v, vtx_t &x )
  {
      splint( x0.n,x1.n, x0.data(),x1.data(), y.data(),y2.data(), u,v, x );
  }

   void ssurf_t::pos( REAL_ u, REAL_ v, vtx_t &x, vtx_t *dx )
  {
      splint( x0.n,x1.n, x0.data(),x1.data(), y.data(),y2.data(), u,v, x,dx );
  }

   void ssurf_t::pos( REAL_ u, REAL_ v, vtx_t &x, vtx_t *dx, vtx_t *d2x )
  {
      splint( x0.n,x1.n, x0.data(),x1.data(), y.data(),y2.data(), u,v, x,dx,d2x );
  }

   void ssurf_t::prj( vtx_t &x0, REAL_ *y, vtx_t &x, stat_t &stat )
  {

      REAL_   rlx;
      REAL_   err;

      REAL_    d;
      REAL_    b[2];
      REAL_    a[3];

      vtx_t   dx[3];
      vtx_t   d2x[3];

      vtx_t   r;

      INT_    j;

      rlx= 0.5;

      for( j=0;j<stat.mit;j++ )
     {

// Bilinear shape function
         pos( y[0],y[1],x,dx,d2x );

         r= x-x0;

// Distance
         d=  r*r;
         d*= 0.5;

// Derivative of distance wrt parametric coordinates

         b[0]= r*dx[0];
         b[1]= r*dx[1];

// Second derivatives

         a[0]=  dx[0]*dx[0]+ r*d2x[0];
         a[1]=  dx[0]*dx[1]+ r*d2x[2];
         a[2]=  dx[1]*dx[1]+ r*d2x[1];

// LDL' factorization and solution

         ldl2( a,b );

// Newton update

         y[0]-= rlx*b[0];
         y[1]-= rlx*b[1];

// Convergence criterion

         err= fabs(b[0])+ fabs(b[1]);

//       printf( "%9.3e % 12.5e % 12.5e % 12.5e % 12.5e % 12.5e % 12.5e\n", err, x0[0],x0[1],x0[2], x[0],x[1],x[2] );
         if( err < stat.tol ){ break; };

         rlx*= 1.2;
         rlx= fmin( rlx,1. );

     }

      return;
  }

   void ssurf_t::nearest( vtx_t x, REAL_ &a, REAL_ &b )
  {
      REAL_ dmin;
      REAL_ amin;
      vtx_t z;
      REAL_ d;

      dmin= 9999;
      b= 0.5;
      for( INT_ i=0;i<64;i++ )
     {
         a= i;
         a/= 63;
         pos( a,b,  z );
   
         d= length( x-z );       
         if( d < dmin  )
        {
            amin= a;
            dmin= d;
        }
          
     }
      a= amin;
  }
