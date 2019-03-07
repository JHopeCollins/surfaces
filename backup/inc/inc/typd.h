#  ifndef _TYPD_
#  define _TYPD_

/*3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
         1         2         3         4         5         6         7         8         9         0         1         2

   Author          Luca di Mare <l.di.mare@ic.ac.uk>
   Created         Wed Jul 14 18:11:25 BST 2010
   Changes History -
   Next Change(s)  ( work in progress )
 */

#    define  VLEN       8
             
#    define  INT_       int
#    define  UINT_      unsigned int
#    define  REAL_      double


#    define  REAL_MPI   MPI_DOUBLE

     typedef INT_        INT_2[2];
     typedef INT_        INT_3[3];
     typedef INT_        INT_4[4];
     typedef INT_        INT_6[6];
     typedef INT_        INT_8[8];

     typedef REAL_      REAL_2[2];
     typedef REAL_      REAL_3[3];
     typedef REAL_      REAL_4[4];
     typedef REAL_      REAL_8[8];
     typedef REAL_      REAL_9[9];

     typedef REAL_      REAL_V[VLEN];

#  endif
