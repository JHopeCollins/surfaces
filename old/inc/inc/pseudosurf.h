#include <math.h>
#include <vec.h>


// virtual base for surface in 3D space defined by x and y coordinates (non-overlapping z)
   struct pseudo3surf_t : public vec_t
  {
      pseudo3surf_t(){ dims=3; construct_levi_civita(); }

      virtual ~pseudo3surf_t(){}

      virtual double F(     double x, double y ) = 0;
      virtual double dFx(   double x, double y ) = 0;
      virtual double dFy(   double x, double y ) = 0;
      virtual double d2Fxx( double x, double y ) = 0;
      virtual double d2Fyy( double x, double y ) = 0;
      virtual double d2Fxy( double x, double y ) = 0;

      inline double F(     double *p ){ return F(     p[0], p[1] ); }
      inline double dFx(   double *p ){ return dFx(   p[0], p[1] ); }
      inline double dFy(   double *p ){ return dFy(   p[0], p[1] ); }
      inline double d2Fxx( double *p ){ return d2Fxx( p[0], p[1] ); }
      inline double d2Fyy( double *p ){ return d2Fyy( p[0], p[1] ); }
      inline double d2Fxy( double *p ){ return d2Fxy( p[0], p[1] ); }

      void tangent( double *p, double *t0, double *t1 );
      void normal(  double *p, double *n );
  };


// z = a + bx + cy + dxy + ex^2 + fy^2
   struct quadratic_t : public pseudo3surf_t
  {
      double a, b, c, d, e, f;

      quadratic_t(){}
     ~quadratic_t(){}

      quadratic_t( double aa, double bb, double cc, double dd, double ee, double ff )
            : a(aa), b(bb), c(cc), d(dd), e(ee), f(ff) {}

      double F(     double x, double y );
      double dFx(   double x, double y );
      double dFy(   double x, double y );
      double d2Fxx( double x, double y );
      double d2Fyy( double x, double y );
      double d2Fxy( double x, double y );

//    double F(     double *p ){ return F(     p[0], p[1] ); }
//    double dFx(   double *p ){ return dFx(   p[0], p[1] ); }
//    double dFy(   double *p ){ return dFy(   p[0], p[1] ); }
//    double d2Fxx( double *p ){ return d2Fxx( p[0], p[1] ); }
//    double d2Fyy( double *p ){ return d2Fyy( p[0], p[1] ); }
//    double d2Fxy( double *p ){ return d2Fxy( p[0], p[1] ); }
  };


// z = a*cos(2*pi*x*lx) * b*sin(2*pi*y*ly)
   struct sinusoid_t : public pseudo3surf_t
  {
      double  a,  b;
      double lx, ly;
      double pi;
      double pi2;
      double pi2lxa;
      double pi2lyb;
      double pi2lxa2;
      double pi2lyb2;

      sinusoid_t(){}
     ~sinusoid_t(){}

      sinusoid_t( double aa, double llx, double bb, double lly )
            : a(aa), b(bb), lx(llx), ly(lly)
     {
         pi= 4.0*atan(1.0);
         pi2=2.*pi;
         pi2lxa=pi2*lx*a;
         pi2lyb=pi2*ly*b;
         pi2lxa2=pi2lxa*pi2lxa;
         pi2lyb2=pi2lyb*pi2lyb;
     }

      double F(     double x, double y );
      double dFx(   double x, double y );
      double dFy(   double x, double y );
      double d2Fxx( double x, double y );
      double d2Fyy( double x, double y );
      double d2Fxy( double x, double y );

//    double F(     double *p ){ return F(     p[0], p[1] ); }
//    double dFx(   double *p ){ return dFx(   p[0], p[1] ); }
//    double dFy(   double *p ){ return dFy(   p[0], p[1] ); }
//    double d2Fxx( double *p ){ return d2Fxx( p[0], p[1] ); }
//    double d2Fyy( double *p ){ return d2Fyy( p[0], p[1] ); }
//    double d2Fxy( double *p ){ return d2Fxy( p[0], p[1] ); }

  };

