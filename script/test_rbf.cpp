#include <iostream>
#include <assert.h>
#include <math.h>
#include <rbf.h>
#include <rbf_interp.h>
#include <surf2.h>

#define SCALE 0.1

   const double off=1.2342;

   double f( double x )
  {
//    return x+off;
//    return x*x;
      return log( x+off );
  }

   double df( double x )
  {
//    return 1;
//    return 2.*x;
      return 1./(x+off);
  }

   int main()
  {
      const int      dims =1;
      const int      npts =30;
      const int      ntest=10;
      const double   xmin=0.0;
      const double   xmax=5.0;


      vtx_t        *x;
      double       *y;

      vtx_t         *xtest;
      double        *ytest;
      double      **dytest;

      double      dx, x0, xt, y0, yt;

      int         i;


   // allocate arrays
      x      = new vtx_t  [npts];
      y      = new double [npts];
      xtest  = new vtx_t  [ntest];
      ytest  = new double [ntest];
      dytest = new double*[ntest];

      for( i=0; i<ntest; i++ ){ dytest[i] = new double[dims]; }

   // absicca and interpolation values
      dx= ( xmax - xmin )/(npts-1.);
      for( i=0; i<npts; i++ )
     {
         x0 =  xmin + i*dx;
         x[i][0] = x0;
         y[i] = f( x0 );
     }

      dx= ( xmax - xmin )/(ntest-1.);
      dx= ( xmax - xmin - dx ) / (ntest-1.);
      x0= xmin+0.5*dx;
      for( i=0; i<ntest; i++ ){ xtest[i][0] = x0+i*dx; }

   // build interpolation
      rbf_interp<1, vtx_t, double> interpolation( npts, x, y );

      i = interpolation.build_weights();
      std::cout << "info: " << i << std::endl;

      // std::cout << "w" << std::endl;
      // for( i=0; i<npts; i++ ){ std::cout << interpolation.w[i] << std::endl; }

      interpolation(            ntest, xtest,  ytest );
      interpolation.derivative( ntest, xtest, dytest );


   // print arrays
      // std::cout << "x,  y" << std::endl;
      // for( i=0; i<npts; i++ ){ std::cout << x[i][0] << " , " << y[i] << std::endl; }

      std::cout << "xt,         y,          yt          error (%)" << std::endl;
      for( i=0; i<ntest; i++ )
     {
         xt = xtest[i][0];
         yt = ytest[i];
         y0 = f( xt );
         std::cout << xt << " ,   " << y0 << " ,   " << yt << " ,   " << 0.01*int(10000.*abs((y0-yt)/y0)) << std::endl;
     }

      std::cout << "xt,        dy,         dyt          error (%)" << std::endl;
      for( i=0; i<ntest; i++ )
     {
         xt =  xtest[i][0];
         yt = dytest[i][0];
         y0 = df( xt );
         std::cout << xt << " ,   " << y0 << " ,   " << yt << " ,   " << 0.01*int(10000.*abs((y0-yt)/y0)) << std::endl;
     }

   // deallocate arrays
      for( i=0; i<ntest; i++ ){ delete[] dytest[i]; }

      delete[]  x;
      delete[]  y;
      delete[]  xtest;
      delete[]  ytest;
      delete[] dytest;

      return 0;
  }

