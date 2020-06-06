#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* FLAME Clock to time our operations.
   Taken from intel's BLIS framework. */
double FLA_Clock();

/* MT-19937 random number generator to generate
   random matrices. It signature is : m, n, A, ldA */
void RandomMatrix(int, int, double *, int);


/* Calculates maximum absolute difference between
   two matrices. Used to compare our impementations
   with BLAS and LAPACK */
double MaxAbsDiff(int, int,
                  double *, int,
                  double *, int);

/* Out implementation of Rank-1 updates. */
void MyGer(int, int,                      // m, n
           double *, int,                 // x, incx
           double *, int,                 // y, incy
           double *, int);                // A. ldA


/* BLIS routine for General Rank-1 update.
   We compare our algorithm with this implementation. */
void dger_(int *, int *,                  // m, n
           double *, double *, int *,     // alpha, x, incx
                     double *, int *,     //        y, incy
           double *, int *);              // A, ldA


int main(int argc, char *argv[])
{
    if ( argc != 4 ) {
        fprintf(stderr, "usage: ./driver.x <start> <last> <inc>\n");
        return 1;
    }

    int
        size, m, n, i = 1,
        first, last, inc,
        i_one = 1, lda,
        AllCorrect = 1;

    double
        diff, maxdiff = 0.0,
        d_one = 1.0;

    double
        *x, *y, *A,
        *Aold, *Aref;

    first = atoi( argv[1] );
    last  = atoi( argv[2] );
    inc   = atoi( argv[3] );

    last  = ( last / inc ) * inc;
    first = ( first / inc ) * inc;
    first = ( first == 0 ? inc : first );

    fprintf(stdout, "\n");

    for( size=last ; size>=first ; size-=inc ) {
        m = n = size;
        lda = size;

        A    = ( double * ) malloc( lda * n * sizeof( double ) );
        Aold = ( double * ) malloc( lda * n * sizeof( double ) );
        Aref = ( double * ) malloc( lda * n * sizeof( double ) );
        x    = ( double * ) malloc(   m * 1 * sizeof( double ) );
        y    = ( double * ) malloc(   n * 1 * sizeof( double ) );

        RandomMatrix( m, n, Aold, lda );
        RandomMatrix( m, 1,    x,   m );
        RandomMatrix( n, 1,    y,   n );

        memcpy( Aref, Aold, lda * n * sizeof( double ) );

        /* Call the LAPACK routine. */
        dger_(&m, &n,
              &d_one, x, &i_one,
                      y, &i_one,
              Aref, &lda);

        memcpy( A, Aold, lda * n * sizeof( double ) );

        MyGer(m, n, x, 1, y, 1, A, lda);

        diff = MaxAbsDiff(m, n, A, lda, Aref, lda);
        maxdiff = ( diff > maxdiff ? diff : maxdiff );

        if ( maxdiff > 1e-10 ) {
            fprintf(stdout, "FAILED: test case %5d\tsize: %5d\tdiff: %8.4le\n", i, size, diff);
            AllCorrect = 0;
        }

        free(    A );
        free( Aold );
        free( Aref );
        free(    x );
        free(    y );

        i++;
    }

    if( AllCorrect ) fprintf(stdout, "LGTM!\n\n");
    else             fprintf(stdout, "some tests fail!\n\n");

    return 0;
}
