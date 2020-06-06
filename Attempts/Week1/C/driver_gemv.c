#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define dabs( x ) ( (x) < 0 ? -(x) : x )

double FLA_Clock();      // This is a routine for extracting elapsed
			 // time borrowed from the libflame library

/* MaxAbsDiff computes the maximum absolute difference over all
   corresponding elements of two matrices */
double MaxAbsDiff( int, int, double *, int, double *, int );

/* RandomMatrix overwrites a matrix with random values */
void RandomMatrix( int, int, double *, int );

/* BLAS implementation of general matrix-vector multiplication.
   Used as a referance to compare with our implementation. */ 
void dgemv_(char *,                         // transpose,
            int *, int *,                   // m, n
            double *, double *, int *,      // alpha, A, ldA
                      double *, int *,      //        x, incx
            double *, double *, int *);     //  beta, y, incy

void MyGemv(int, int,                       // m, n
            double *, int,                  // A, ldA
            double *, int,                  // x, incx
            double *, int);                 // y, incy

int main(int argc, char *argv[])
{
    if ( argc != 4 ) {
        fprintf(stderr, "usage: ./driver.x <first> <last> <inc>\n");
        return 1;
    }

    int
        m, n,
        size, first, last, inc,
        lda, i_one = 1,
        AllCorrect = 1;

    double
        d_one = 1,
        diff, maxdiff = 0.0;

    double
        *A, *x, *y,
        *yold, *yref;

    first = atoi( argv[1] );
    last  = atoi( argv[2] );
    inc   = atoi( argv[3] );

    last = ( last / inc ) * inc;
    first  = ( first / inc ) * inc;
    first = ( first == 0 ? inc : first );

    printf("\n");

    for( size=last ; size>=first ; size-=inc ) {
        m = n = size;
        lda = size;

        A    = ( double * ) malloc( lda * n * sizeof( double ) );
        x    = ( double * ) malloc(   n * 1 * sizeof( double ) );
        yold = ( double * ) malloc(   m * 1 * sizeof( double ) );
        yref = ( double * ) malloc(   m * 1 * sizeof( double ) );
        y    = ( double * ) malloc(   m * 1 * sizeof( double ) );

        RandomMatrix( lda, n,    A, lda );
        RandomMatrix(   n, 1,    x,   n );
        RandomMatrix(   m, 1, yold,   m );

        memcpy( yref, yold, m * 1 * sizeof( double ) );

        dgemv_("No transpose",
               &m, &n,
               &d_one,    A, &lda,
                          x, &i_one,
               &d_one, yref, &i_one);

        memcpy( y, yold, m * 1 * sizeof( double ) );

        MyGemv(m, n,
               A, lda,
               x, i_one,
               y, i_one);

        diff = MaxAbsDiff(m, 1, yref, m, y, m);
        maxdiff = ( diff > maxdiff ? diff : maxdiff );

        if( maxdiff > 1e-08 ) {
            AllCorrect = 0;
            fprintf(stderr, "FAILED: size=%5d, diff=%8.4le\n", size, diff);
        }

        free(    A );
        free(    x );
        free(    y );
        free( yold );
        free( yref );
    }

    if( AllCorrect ) fprintf(stdout, "LGTM!\n\n");
    else {
        fprintf(stderr, "some tests fail!\n\n");
        return 1;
    }

    return 0;
}
