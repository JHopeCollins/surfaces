# ifndef _VSIZE_
# define _VSIZE_

# include <typd.h>
# include <misc.h>

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         
// Changes History
// Next Change(s)  -

/** A data structure to store variable size linear arrays with direct addressing **/

  template <typename type,size_t DSIZE> struct vsize_t
 {
     private:

        type         *v;
        INT_          m; 

     public: 

        INT_          n;
   
/** Constructor **/

        vsize_t()
       {
           n=0;
           m=0;
           v=NULL;
       };
   
/** Access to the storage **/

        type *data(){ return v; };
   
/** Destructor **/

       ~vsize_t()
       {
           n=0;
           m=0;
           free(v); v=NULL;
       };



/** Destructor **/

        void clear()
       {
           n=0;
           m=0;
           free(v); v=NULL;
       };

/** Copy constructor **/

        vsize_t &operator=( vsize_t data )
       {
           clear();
           n= data.n; resize( data[0]);
           assert( m == data.m );
           memcpy( v,data.v, m*sizeof(type) );
           return *this;
       }


/*      void copy( vsize_t data )
       {
//         clear();
           n= data.n; resize( data[0] );
           assert( m == data.m );
//         memcpy( v,data.v, m*sizeof(type) );
//         for( INT_ i=0;i<n;i++ ){ v[i]= data[i]; };
           return;
       };*/
   
/** Constructor **/

        vsize_t( INT_ l, type data)
       {
           n=l;
           m=pad(l,DSIZE);
           v=(type*)malloc(m*sizeof(type));
           for( INT_ i=0;i<m;i++ ){ v[i]= data; };
       };
   
/** Resize storage to fit current logical size and initialize values **/

        void resize( type data )
       {
           resize( n,data );
           return ;
       }
   
/** Resize storage to fit l size and initialize values **/

        void resize( INT_ l, type data )
       {
           if( l >= m )
          {
              INT_ o=m;
              m= pad(l,DSIZE); 
              v= (type*)realloc(v,m*sizeof(type));
              for( INT_ i=o;i<m;i++ ){ v[i]= data; };
          }
           return ;
       }

/** Overload operator [] so that the vsize can be accessed like any array **/
   
        type &operator[]( INT_ i )
       { 
           assert( i > -1 && i < n );
           return v[i]; 
       }; 

/** Add entry and return its position **/

        INT_ append( type val )
       {
           INT_ i= n++;
           resize( val );
           return i;
       }
     
/** Add entry at position k **/

        void  insert( INT_ k, type val )
       {
           if( k >= n ){ n=k+1; resize(val); };
       };

/** Test for entry **/
        inline INT_ has( type var )
       {
           INT_ val=-1;
           for( INT_ i=0;i<n;i++ )
          {
              if( v[i] == var )
             {
                 val= i;
                 break;
             }
          }
           return val;
       };

        inline void compare( INT_ ist, INT_ ien, vsize_t<type,DSIZE> &var, INT_ jst, INT_ jen, INT_ &h, INT_ &k )
       {
           INT_ i,j;
           h=-1;
           k=-1;
           for( i=ist;i<ien;i++ )
          {
              for( j=jst;j<jen;j++ )
             {
                 if( v[i] == var[j]  )
                {
                    h= i;
                    k= j;
                    break;
                }
             }
          }
       }


 };

//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//       1         2         3         4         5         6         7         8         9         0         1         2
//
// Author          Luca di Mare <l.di.mare@ic.ac.uk>
// Created         
// Changes History
// Next Change(s)  -

/** A data structure to store two-dimensinoal, variable size linear arrays with direct addressing **/

  template <typename type,INT_ V,size_t DSIZE> struct wsize_t
 {
        typedef type  vtype[V];
     private:
        vtype        *v;
        INT_          m; 

     public:
        INT_          n;

/** Constructor **/

        wsize_t()
       {
           n=0;
           m=0;
           v=NULL;
       };

/** Access storage **/   

       ~wsize_t()
       {
           n=0;
           m=0;
           free(v); v=NULL;
       };

/** Access storage **/   

        vtype *data(){ return (vtype *)v; };

/** Constructor **/   
   
        wsize_t( INT_ l, type data )
       {
           n=l;
           m=pad(l,DSIZE);
           v=(vtype*)malloc(m*sizeof(vtype));
           for( INT_ i=0;i<m;i++ ){ for( INT_ j=0;j<V;j++ ){ v[i][j]= data; }; }
       };

/** Resize to fit logical size **/   
   
        void resize( type data )
       {
           resize( n,data );
           return;
       }
   
/** Resize to fit l entries **/   

        void resize( INT_ l, type data )
       {
           if( l >= m )
          {
              INT_ o=m;
              m= pad(l,DSIZE); 
              v= (vtype*)realloc(v,m*sizeof(vtype));
              for( INT_ i=o;i<m;i++ ){ for( INT_ j=0;j<V;j++ ){ v[i][j]= data; } }
          }
           return ;
       }
   
/** Overload operator [] **/

        type* operator[]( INT_ i )
       { 
           assert( i > -1 && i < n );
           return v[i]; 
       }; 

/** Destructor **/

        void clear()
       {
           n=0;
           m=0;
           free(v); v=NULL;
       };

/*      wsize_t &operator=( wsize_t data )
       {
           clear();
           n= data.n; resize( data[0][0] );
           assert( m == data.m );
           memcpy( v,data.v, m*sizeof(vtype) );
           return *this;
       }*/

        void copy( wsize_t data )
       {
//         clear();
           n= data.n; resize( data[0][0] );
           assert( m == data.m );
//         memcpy( v,data.v, m*sizeof(vtype) );
//         for( INT_ i=0;i<n;i++ ){ for( INT_ j=0;j<V;j++ ){ v[i][j]= data[i][j]; }};
           return;
       };

/** Add entry and return its position **/

        INT_ append( type val )
       {
           INT_ i= n++;
           resize( val );
           return i;
       }

/** Add entry at position k **/

       void  insert( INT_ k, type val )
     {
         if( k >= n ){ n=k+1; resize(val); };
     };

/** Test for entry **/
        inline INT_ has( INT_ j, type var )
       {
           INT_ val=-1;
           for( INT_ i=0;i<n;i++ )
          {
              if( v[i][j] == var )
             {
                 val= i;
                 break;
             }
          }
           return val;
       };

        inline void compare( INT_ ist, INT_ ien, INT_ l, wsize_t<type,V,DSIZE> &var, INT_ jst, INT_ jen, INT_ m, INT_ &h, INT_ &k )
       {
           INT_ i,j;
           h=-1;
           k=-1;
           for( i=ist;i<ien;i++ )
          {
              for( j=jst;j<jen;j++ )
             {
                 if( v[i][l] == var[j][m]  )
                {
                    h= i; 
                    k= j; 
                    break;
                }
             }
          }
       }
     
 };


# endif
