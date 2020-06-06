#define chi( i ) x[ (i)*incx ]
#define psi( i ) y[ (i)*incy ]

void Axpy(int m,
          double alpha, double *x, int incx,
          double *y, int incy)
/**
 * Axpy operation.
 * computes `y = y + alpha*x`
 */
{
    int i;
    for( i=0 ; i<m ; i++ )
        psi( i ) += alpha * chi( i );
}
