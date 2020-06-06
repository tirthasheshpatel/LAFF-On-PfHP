#define alpha( i,j ) A[ (j)*lda + j ]
#define     chi( i ) x[ (i)*incx ]
#define     psi( i ) y[ (i)*incy ]

void Dots(int,
          const double *, int,
          const double *, int,
          double *);

void MyGemv(int m, int n,
            double *A, int lda,
            double *x, int incx,
            double *y, int incy)
/**
 * General matrix vector multiplication with scalar updates.
 * Computes `y = Ax + y` using dot products with scalar updates.
 * Ax + y = | A_0 * x + y_0 |
 *          | A_1 * x + y_1 |
 *          |      ...      |
 *          | A_m * x + y_m |
 */
{
    int i;
    for( i=0 ; i<m ; i++ )
        Dots( n,
              &A[i], lda,
              x, incx,
              &psi( i ) );
}
