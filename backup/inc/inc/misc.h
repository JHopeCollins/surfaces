#  ifndef _MISC_
#  define _MISC_

#  include <typd.h>

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Swap a and b **/

   template <typename type> void swap( type &a, type &b )
  {
      type c;
      c= a;
      a= b;
      b= c;
  };

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Find maximum between a and b **/

   template <typename type> type max( type a, type b ){ return (a>b)?a:b; };
   template <typename type> type min( type a, type b ){ return (a<b)?a:b; };

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Return 1 if v is positive, -1 otherwise **/

   template <typename type> type sign( type v ){ return (v>0)?1:-1; };

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Find the smallest multiple of j larger than i **/

   inline INT_ pad( INT_ i,  INT_ j )
  {
      INT_ v;
      v= i/j;
      if( i%j > 0 ){ v++; }
      v*= j;
      return v;
  };

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Sort a short array with three entries with no more than three comparisons **/

   template <typename type> void sort3( type *v, INT_ *prm )
  {
      prm[0]=0;
      prm[1]=1;
      prm[2]=2;

      if( v[prm[0]] > v[prm[1]] )
     { 
         swap(prm[0],prm[1]); 
     }
      if( v[prm[1]] > v[prm[2]] )
     { 
         swap(prm[1],prm[2]); 
         if( v[prm[0]] > v[prm[1]] )
        { 
            swap(prm[1],prm[0]); 
        }
     }
  }

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Sort a short array with four entries with no more than three comparisons **/

   template <typename type> void sort4( type *v, INT_ *iprm )
  {
      iprm[0]=0;
      iprm[1]=1;
      iprm[2]=2;
      iprm[3]=3;

      if( v[iprm[0]]>v[iprm[1]] ){ swap( iprm[0],iprm[1] ); };
      if( v[iprm[2]]>v[iprm[3]] ){ swap( iprm[2],iprm[3] ); };

      if( v[iprm[1]]>v[iprm[2]] )
     { 
         swap( iprm[1],iprm[2] ); 
         if( v[iprm[0]]>v[iprm[1]] ){ swap( iprm[0],iprm[1] ); };
         if( v[iprm[2]]>v[iprm[3]] ){ swap( iprm[2],iprm[3] ); };
     };
      if( v[iprm[1]]>v[iprm[2]] ){ swap( iprm[1],iprm[2] ); }
      return;

  }

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         Tue  7 Aug 09:37:52 BST 2018
// Changes History -
// Next Change(s)  -

/** Binary search **/

/* template < typename type > INT_ bsearch( type var, INT_ n, type *data )
  {
      INT_ val=-1;
      INT_ i0,i1;
      if( n > 0 )
     {
        if( var <= data[0] )
       {
           val=0; 
       }
        else
       {
           if( var >= data[n-1] )
          {
              val= n-1;
          }
           else
          {
              i0=0;
              i1=n-1;
              val= i0+i1;
              val/= 2;
              while( true )
             {
                 if( var < data[val] )
                {
                    i1= val;
                }
                 else
                {
                    i0= val;
                }
                 if( i1 <= i0+1 )
                {
                    val= i1;
                    break;
                }
                 else
                {
                    val= i0+i1;
                    val/= 2;
                }
              }
           }
        }
     }
      return val;
  }*/

#  endif
