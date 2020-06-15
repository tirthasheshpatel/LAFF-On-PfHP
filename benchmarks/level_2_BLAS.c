/**
 * The naming convention for Level-2 BLAS routines is given by
 *                     _XXYY
 * where _ = double precition, single precition, etc.
 *       XX = shape of the matrix ==> 'ge' for general,
 *                                    'sy' for symetric,
 *                                    'he' for hermitian,
 *                                    'tr' for triangular.
 *       YY = operation to be performed ==> mv = matrix vector multiply,
 *                                          sv = solve vector,
 *                                           r = rank-1 update,
 *                                          r2 = rank-2 update.
 * I have benchmarked the most used BLAS functions.
 *  - gemv : general matrix-vector multiplication
 *  - symv : symmetric matrix-vector multiplication
 *  - trmv : triangular matrix-vector multiplication
 *  - trsv : triangular solve vector
 *  - ger  : general rank-1 update
 *  - syr  : symmetric rank-1 update
 *  - syr2 : symmetric rank-2 update
 */

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

/* Generates a Random Lower Triangular Matrix using a MT-19937 RNG.
   Signature : (m, A, ldA) */
void RandomTril(int, double *, int);

/* Generates a Random Symmetric Matrix using a MT-19937 RNG.
   Signature : (m, A, ldA) */
void RandomSym(int, double *, int);

/* General Matrix Vector Multiply operation.
   Signature : (trans, m, n, alpha, A, ldA, x, incx, beta, y, incy) */
void dgemv_(char *,
            int *, int *,
            double *, double *, int *,
                      double *, int *,
            double *, double *, int *);

/* Symmetric Matrix Vector Multiply operation.
   Signature : (uplo, m, alpha, A, ldA, x, incx, beta, y, incy) */
void dsymv_(char *, int *,
            double *, double *, int *,
                      double *, int *,
            double *, double *, int *);

/* Triangular Matrix Vector Solve operation.
   Signature : (uplo, trans, diag, m, A, ldA, x, incx) */
void dtrsv_(char *, char *, char *,
            int *,
            double *, int *,
            double *, int *);

/* General Rank-1 Update operation.
   Signature : (m, n, alpha, x, incx, y, incy, A, ldA) */
void dger_(int *, int *,
           double *, double *, int *,
                     double *, int *,
           double *, int *);

/* Symmetric Rank-1 Update operation.
   Signature : (uplo, m, alpha, x, incx, A, ldA) */
void dsyr_(char *, int *,
           double *, double *, int *,
           double *, int *);

/* Symmetric Rank-2 Update operation.
   Signature : (uplo, m, alpha, x, incx, y, incy, A, ldA) */
void dsyr2_(char *, int *,
            double *, double *, int *,
                      double *, int *,
            double *, int *);


int main(int argc, char *argv[])
{
   if( argc != 6 ) {
      fprintf(stderr, "usage: ./driver.x <opname> <nrepeats> <first> <last> <inc>\n");
      return 1;
   }

   int
      size, i_one = 1,
      first, last, inc,
      ldA, nrep,
      m, n, k, i;

   char
      up = 'U', lo = 'L',
      t = 'T', not = 'N',
      diag = 'N';

   double
      dtime, dtime_best,
      d_one = 1.0;

   double
      *x, *y, *A;

   nrep  = atoi( argv[2] );
   first = atoi( argv[3] );
   last  = atoi( argv[4] );
   inc   = atoi( argv[5] );

   last  = (  last / inc ) * inc;
   first = ( first / inc ) * inc;
   first = ( first == 0 ? inc : first );

   printf("size, time(ms)\n");

   if( strcmp( "gemv", argv[1] ) == 0 ) {
      /* Benchmark dgemv operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         y = ( double * ) malloc(   m * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
            RandomMatrix(m, 1, y,   m);
            RandomMatrix(m, n, A, ldA);

            dtime = FLA_Clock();

            dgemv_("No transpose", &m, &n, &d_one, A, &ldA,
                   x, &i_one, &d_one, y, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
         free(A);
      }
   }
   else if( strcmp( "symv", argv[1] ) == 0 ) {
      /* Benchmark dsymv operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         y = ( double * ) malloc(   m * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
            RandomMatrix(m, 1, y,   m);
               RandomSym(m,    A, ldA);

            dtime = FLA_Clock();

            dsymv_(&lo, &m, &d_one, A, &ldA,
                   x, &i_one, &d_one, y, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
         free(A);
      }
   }
   else if( strcmp( "trsv", argv[1] ) == 0 ) {
      /* Benchmark dtrsv operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
              RandomTril(m,    A, ldA);

            dtime = FLA_Clock();

            dtrsv_(&lo, &not, &diag, &m, A, &ldA, x, &i_one);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(A);
      }
   }
   else if( strcmp( "ger", argv[1] ) == 0 ) {
      /* Benchmark dger operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         y = ( double * ) malloc(   m * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
            RandomMatrix(m, 1, y,   m);
            RandomMatrix(m, n, A, ldA);

            dtime = FLA_Clock();

            dger_(&m, &n, &d_one, x, &i_one, y, &i_one, A, &ldA);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
         free(A);
      }
   }
   else if( strcmp( "syr", argv[1] ) == 0 ) {
      /* Benchmark dsyr operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
               RandomSym(m,    A, ldA);

            dtime = FLA_Clock();

            dsyr_(&lo, &m, &d_one, x, &i_one, A, &ldA);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(A);
      }
   }
   else if( strcmp( "syr2", argv[1] ) == 0 ) {
      /* Benchmark dsyr2 operation. */
      for( size=last ; size>=first ; size-=inc ) {
         m = n = size;
         ldA = size;

         x = ( double * ) malloc(   n * 1 * sizeof( double ) );
         y = ( double * ) malloc(   m * 1 * sizeof( double ) );
         A = ( double * ) malloc( ldA * n * sizeof( double ) );

         for( i=0 ; i<nrep ; i++ ) {
            RandomMatrix(n, 1, x,   n);
            RandomMatrix(m, 1, y,   m);
               RandomSym(m,    A, ldA);

            dtime = FLA_Clock();

            dsyr2_(&lo, &m, &d_one, x, &i_one, y, &i_one, A, &ldA);

            dtime = FLA_Clock() - dtime;
            if( i == 0 )                  dtime_best = dtime;
            else if( dtime < dtime_best ) dtime_best = dtime;
         }

         fprintf(stdout, "%5d, %8.4le\n", size, dtime_best);

         free(x);
         free(y);
         free(A);
      }
   }

   return 0;
}
