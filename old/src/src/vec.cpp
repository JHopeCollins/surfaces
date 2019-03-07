# include <vec.h>

// base struct for classes requiring vector operations

      void vec_t::malloc()
     {
         int   i,j;

         e=NULL;
         e=new double**[dims];

         for( i=0; i<dims; i++ )
        {
            e[i]=NULL;
            e[i]=new double*[dims];

            for( j=0; j<dims; j++ )
           {
               e[i][j]=NULL;
               e[i][j]=new double[dims];
           }
        }
         arrays=1;
     }

      void vec_t::free()
     {
         int   i,j;

         for( i=0; i<dims; i++ )
        {
            for( j=0; j<dims; j++ )
           {
               delete[] e[i][j];
               e[i][j]=NULL;
           }
            delete[] e[i];
            e[i]=NULL;
        }
         delete[] e;
         e=NULL;
     }

      void vec_t::construct_levi_civita()
     {
         int   i,j,k;

         if( !arrays){ malloc(); }

         for( i=0; i<dims; i++ )
        {
            for( j=0; j<dims; j++ )
           {
               for( k=0; k<dims; k++ )
              {
                  if( i==j || i==k || j==k ){           e[i][j][k]= 0.; }
                  if( j==(i+1)%dims && k==(j+1)%dims ){ e[i][j][k]= 1.; }
                  if( j==(i-1)%dims && k==(j-1)%dims ){ e[i][j][k]=-1.; }
              }
           }
        }
     }

