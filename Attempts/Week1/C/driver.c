#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define dabs(x) ( (x) < 0 ? -(x) : x )

double FLA_Clock(); // FLAME Clock to time operations

/**
 * Calculates the maximum absolute difference between two
 * matrices. Used to assert the correctness of the algoirthm
 * with the referance routine.
 */
double MaxAbsDiff(int, int,
                  double *, int,
                  double *, int);

/**
 * Generates a psuedo-random matrix using MT-19937 RNG.
 */
void RandomMatrix(int, int, double *, int);

/**
 * General matrix matirx multiplication implementation
 * done in different forms. All the functions must have same
 * signature which is (m, n, k, ap, lda, bp, ldb, cp, ldc)
 */
void MyGemm(int, int, int,                      // m, n, k
            const double *, int,                // A, ldA
            const double *, int,                // B, ldB
                  double *, int);               // C, ldc

/**
 * BLAS/LAPACK implementation of the routine. This is the refernce
 * routine we are going to compare with.
 */
void dgemm_( char *, char *,                 // transA, transB
             int *, int *, int *,            // m, n, k
             double *, double *, int *,      // alpha, A, ldA
                       double *, int *,      //        B, ldB
             double *, double *, int * );    // beta,  C, ldC

// Tests for general matrix matrix multiplication using
// blis package with openmp.
int main(int argc, char **argv)
{
    if ( argc != 5 ) {
        fprintf(stderr, "usage: ./driver.x <nrepeats> <first> <last> <inc>\n");
        return 1;
    }

    int
        m, n, k,
        lda, ldb, ldc,
        size, first, last, inc,
        i, irep,
        nrepeats;

    double
        d_one = 1.0,
        dtime, dtime_best,
        diff, maxdiff = 0.0, gflops;

    double
        *A, *B, *C,
        *Cold, *Cref;

    nrepeats = atoi(argv[1]);
    first    = atoi(argv[2]);
    last     = atoi(argv[3]);
    inc      = atoi(argv[4]);

    last  = ( last / inc )  * inc;
    first = ( first / inc ) * inc;
    first = ( first == 0 ) ? inc : first;

    fprintf(stdout, "first: %d\tlast: %d\tinc: %d\n\n", first, last, inc);
    i = 1;

    for ( size=last ; size>=first; size-=inc ) {
        m = n = k = size;
        lda = ldb = ldc = size;

        gflops = ( 2. * m * n * k ) * 1e-09;

        A    = ( double * ) malloc( lda * k * sizeof( double ) );
        B    = ( double * ) malloc( ldb * n * sizeof( double ) );
        C    = ( double * ) malloc( ldc * n * sizeof( double ) );
        Cold = ( double * ) malloc( ldc * n * sizeof( double ) );
        Cref = ( double * ) malloc( ldc * n * sizeof( double ) );

        RandomMatrix(m, k,    A, lda);
        RandomMatrix(k, n,    B, ldb);
        RandomMatrix(m, n, Cold, ldc);

        /* Time reference implementation provided by the BLAS library
           routine dgemm (double precision general matrix-matrix
           multiplicationn */
        for ( irep=0 ; irep<nrepeats ; irep++ ) {
            memcpy( Cref, Cold, ldc * n * sizeof( double ) );

            dtime = FLA_Clock();

            dgemm_("No transpose", "No transpose",
                   &m, &n, &k,
                   &d_one, A, &lda,
                           B, &ldb,
                   &d_one, Cref, &ldc);

            dtime = FLA_Clock() - dtime;

            if      ( irep == 0 )           dtime_best = dtime;
            else if ( dtime < dtime_best )  dtime_best = dtime;
        }

        printf("LAPACK   ==> size: %5d\tTime (ms): %8.4le\tGFLOPS/sec: %8.4le\n",
               n, dtime_best, gflops/dtime_best);
        fflush(stdout);

        /* Time our implementation takes to compute the
           general matrix matrix multiplication */
        for ( irep=0 ; irep<nrepeats ; irep++ ) {
            memcpy( C, Cold, ldc * n * sizeof( double ) );

            dtime = FLA_Clock();

            MyGemm(m, n, k,
                   A, lda,
                   B, ldb,
                   C, ldc);

            dtime = FLA_Clock() - dtime;

            if      ( irep == 0 )          dtime_best = dtime;
            else if ( dtime < dtime_best ) dtime_best = dtime;
        }

        diff = MaxAbsDiff(m, n, C, ldc, Cref, ldc);
        maxdiff = ( diff > maxdiff ? diff : maxdiff ) ;

        printf("C NATIVE ==> size: %5d\tTime (ms): %8.4le\tGFLOPS/sec: %8.4le\tdiff: %8.4le\n",
               n, dtime_best, gflops/dtime_best, diff);
        fflush(stdout);

        free( A    );
        free( B    );
        free( C    );
        free( Cold );
        free( Cref );

        i++;
    }

    printf("\n%% Maximum difference between reference and your implementation: %le.\n",
           maxdiff);
    fflush(stdout);

    return 0;
}
