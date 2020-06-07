#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* A very high precition FLAME clock to
   time Level 1 BLAS operations */
double FLA_Clock();

/* Generates a Random Matrix using a MT-19937 RNG.
   Signature : (m, n, A, ldA) */
void RandomMatrix(int, int, double *, int);

/* Double Precition Copy Operation. 
   Signature : (m, x, incx, y, incy) */
void dcopy_(int *, double *, int *, double *, int *);

/* Double Precition Scaling Operation.
   Signature : (m, alpha, x, incx) */
void dscal_(int *, double *, double *, int *);

/* Double Precition Axpy Operation.
   Signature : (m, alpha, x, incx, y, incy) */
void daxpy_(int *, double *, double *, int *, double *, int *);

/* Double Precition Dot Operation.
   Signature : (m, x, incx, y, incy) */
void ddot_(int *, double *, int *, double *, int *);

/* Double Precition L2-Norm Operation.
   Signature : (m, x, incx) */
void dnrm2_(int *, double *, int *);

/* Double Precition Asum (L1-norm) Operation.
   Signature : (m, x, incx) */
void dasum_(int *, double *, int *);


int main(int argc, char *argv[])
{
   if( argc != 6 ) {
      fprintf(stderr, "usage: ./driver.x <opname> <nrepeats> <first> <last> <inc>\n");
      return 1;
   }

   int
      size, i_one = 1,
      first, last, inc,
      nrep, m, n, k, i;

   double
      dtime, dtime_best,
      d_one = 1.0;

   double
      *x, *y, *yold;

   nrep  = atoi( argv[2] );
   first = atoi( argv[3] );
   last  = atoi( argv[4] );
   inc   = atoi( argv[5] );

   last  = (  last / inc ) * inc;
   first = ( first / inc ) * inc;
   first = ( first == 0 ? inc : first );

   // printf("\nfirst: %5d\tlast: %5d\tinc: %5d\n", first, last, inc);

   printf("size, time(ms)\n");

   if( strcmp( "copy", argv[1] ) == 0 ) {
      /* Benchmark dcopy operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(m, 1, x, m);

            dtime = FLA_Clock();

            dcopy_(&m, x, &i_one, y, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
      }
   }
   else if( strcmp( "scal", argv[1] ) == 0 ) {
      /* Benchmark dscal operation. */
      /* We don't repeat otherwise the value of x gets
         cached and we see unwanted results. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         RandomMatrix(m, 1, x, m);

         dtime = FLA_Clock();

         dscal_(&m, &d_one, x, &i_one);

         dtime = FLA_Clock() - dtime;

         fprintf(stdout, "%5d, %8.4le\n", size, dtime);

         free(x);
         free(y);
      }
   }
   else if( strcmp( "axpy", argv[1] ) == 0 ) {
      /* Benchmark daxpy operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(m, 1, x, m);

            dtime = FLA_Clock();

            daxpy_(&m, &d_one, x, &i_one, y, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
      }
   }
   else if( strcmp( "dot", argv[1] ) == 0 ) {
      /* Benchmark ddot operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(m, 1, x, m);

            dtime = FLA_Clock();

            ddot_(&m, x, &i_one, y, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
      }
   }
   else if( strcmp( "nrm2", argv[1] ) == 0 ) {
      /* Benchmark dnrm2 operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(m, 1, x, m);

            dtime = FLA_Clock();

            dnrm2_(&m, x, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
      }
   }
   else if( strcmp( "asum", argv[1] ) == 0 ) {
      /* Benchmark dasum operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = size;

         x = ( double * ) malloc( m * 1 * sizeof( double ) );
         y = ( double * ) malloc( m * 1 * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(m, 1, x, m);

            dtime = FLA_Clock();

            dasum_(&m, x, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
      }
   }
   else {
      fprintf(stderr, "Not implemented!\n");
      return 1;
   }

   return 0;
}
