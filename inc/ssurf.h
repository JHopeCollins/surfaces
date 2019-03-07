#  ifndef _SSURF_
#  define _SSURF_

#  include <cstdio>
#  include <cstdlib>
#  include <cassert>
#  include <cmath>

#  include <vtx.h>
#  include <spline2.h>
#  include <vsize.h>

#  include <tran.h>

   struct stat_t
  {
      INT_   mit;
      REAL_  tol;

      stat_t()
     {
         mit=20;
         tol=1.e-24;
     }
  };

   struct ssurf_t
  {
      vsize_t<REAL_,32>            x0;
      vsize_t<REAL_,32>            x1;
      vsize_t<vtx_t,32>            y,y2;

      void read( char *name, tran_t &r, REAL_ h );
      void write( INT_, char *name );
      void test( char *name );
      void pos( REAL_ u, REAL_ v, vtx_t &x );
      void pos( REAL_ u, REAL_ v, vtx_t &x, vtx_t *dx );
      void pos( REAL_ u, REAL_ v, vtx_t &x, vtx_t *dx, vtx_t *d2x );

      void prj( vtx_t &x0, REAL_ *y, vtx_t &x, stat_t &stat );
      void nearest( vtx_t x, REAL_ &a, REAL_ &b );
  };

#  endif
