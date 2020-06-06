#define alpha( i,j ) A[ (j)*lda + (i) ]
#define     chi( i ) x[ (i)*incx ]
#define     psi( i ) y[ (i)*incy ]

void Axpy(int, double, double *, int, double *, int);

void MyGemv(int m, int n,
            double *A, int lda,
            double *x, int incx,
            double *y, int incy)
/**
 * General Matrix-Vector Multiplication using Axpy Operation.
 * Computes `y = Ax + y` where `x,y` are vectors and `A` is a matrix.
 * Ax + y = | A(   0,0 )*x( 0 ) | + | A(   0,1 )*x( 1 ) | + ... + | A(   0,n-1 )*x( n-1 ) | + |   y( 0 ) |
 *          | A(   1,0 )*x( 0 ) | + | A(   1,1 )*x( 1 ) | + ... + | A(   1,n-1 )*x( n-1 ) | + |   y( 1 ) |
 *          |        ...        | + |       ...         | + ... + |          ...          | + |    ...   |
 *          | A( m-1,0 )*x( 0 ) | + | A( m-1,1 )*x( 1 ) | + ... + | A( m-1,n-1 )*x( n-1 ) | + | y( m-1 ) |
 */
{
    int j;
    for( j=0 ; j<n ; j++ )
        Axpy(m, chi( j ), &alpha( 0,j ), 1, y, incy);
}
