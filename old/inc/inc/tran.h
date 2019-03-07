#  ifndef _TRAN_
#  define _TRAN_

#  include <vtx.h>

   const vtx_t vdum;

   struct tran_t
  {
      vtx_t                        x0;
      REAL_                        c,s,d;

      tran_t()
     {
         x0= vtx_t(0.,0.,0.);
         c=  1;
         s=  0;
         d=  1;
     };

      inline void tran( vtx_t &v )
     {
         REAL_ w;    
         w= v[0];
         v[0]= c*w- s*v[1];
         v[1]= s*w+ c*v[1];
         v*= d;
         v+= x0;
         return;
     };
  };

#  endif
